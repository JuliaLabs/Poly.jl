using JuLoop
using ISL
import ISL.API
using Printf

export create_context_and_universe


"""
create an ISL basic set representing the domains in kernel
"""
function create_context_and_universe(kernel::LoopKernel)#::Tuple{ISL.API.isl_basic_set, ISL.API.isl_ctx}
    context = ISL.API.isl_ctx_alloc()
    println("GOT CONTEXT")
    # @show domains_isl_rep(kernel)
    space = ISL.API.isl_set_read_from_str(context, "{:}")
    println("GOT SPACE")
    @show instructions_isl_rep(kernel)
    instructions = ISL.API.isl_union_set_read_from_str(context, instructions_isl_rep(kernel))
    println("GOT INSTS")
    @show dependencies_isl_rep(kernel)
    dependencies = ISL.API.isl_union_map_read_from_str(context, dependencies_isl_rep(kernel))
    println("GOT DEPS")
    schedule_map = ISL.API.isl_union_map_intersect_domain(dependencies, instructions)
    println("GOT MAP")
    build = ISL.API.isl_ast_build_from_context(space)
    println("GOT BUILD")
    ast = ISL.API.isl_ast_build_node_from_schedule_map(build, schedule_map)
    println("GOT AST")
    print(ISL.API.isl_ast_node_dump(ast))
end


"""
helper to turn gensyms into strings
"""
function sym_to_str(s::Symbol)::String
   str = string(s)
   str = replace.(str, r"#" => "")
   str = replace.(str, r":" => "")
   return str
end


"""
modify loop domains to be in the ISL syntax form
ex: i = 1:10 becomes [] -> {[i]: 1 <= i <= 10}
ex: i = 1:2:n becomes  [n] -> {[i]: exists a: i = 2a and 1 <= i <= n}
"""
function domains_isl_rep(kernel::LoopKernel)::String
    keys = string([domain.iname for domain in kernel.domains])
    keys = replace.(keys, r":" => "")

    params = string(kernel.args)
    params = replace.(params, r":" => "")
    params = replace.(params, r"Symbol" => "")

    conditions = ""
    count = 1
    for domain in kernel.domains
        step = domain.recurrence.args[2]
        if count != 1
            conditions = string(conditions, " and ")
        end
        count += 1
        if step != 1
            conditions = string(conditions, @sprintf("exists (a: %s = %da and %s <= %s <= %s)", string(domain.iname), step, string(domain.lowerbound), string(domain.iname), string(domain.upperbound)))
        else
            conditions = string(conditions, @sprintf("%s <= %s <= %s", string(domain.lowerbound), string(domain.iname), string(domain.upperbound)))
        end
    end

    return @sprintf("%s -> {%s: %s}", params, keys, conditions)
end


"""
construct the instructions representation in ISL syntax form
ex: (id = :mult) C[i, j] += A[i, k] * B[k, j] in matrix multiplication becomes
{mult[i, j, k] : 0 <= i <= n and 0 <= j <= m and 0 <= k <= r}
"""
function instructions_isl_rep(kernel::LoopKernel)::String
    params = string(kernel.args)
    params = replace.(params, r":" => "")
    params = replace.(params, r"Symbol" => "")

    insn_map = "{"
    icount = 1
    for instruction in kernel.instructions
        count = 1
        name = string(sym_to_str(instruction.iname), "[")
        conditions = ""
        for iname in instruction.dependencies
            for domain in kernel.domains
                if domain.iname == iname
                    step = domain.recurrence.args[2]
                    if count != 1
                        conditions = string(conditions, " and ")
                        name = string(name, ", ")
                    end
                    name = string(name, domain.iname)
                    lb = domain.lowerbound
                    ub = domain.upperbound
                    count += 1
                    if step != 1
                        conditions = string(conditions, @sprintf("exists (a: %s = %da and %s <= %s <= %s)", string(domain.iname), step, string(lb), string(domain.iname), string(ub)))
                    else
                        conditions = string(conditions, @sprintf("%s <= %s <= %s", string(lb), string(domain.iname), string(ub)))
                    end
                end
            end
        end
        name = string(name, "] : ")
        if conditions != ""
            if icount != 1
                insn_map = string(insn_map, "; ")
            end
            icount += 1
            insn_map = string(insn_map, name, conditions)
        end
    end
    insn_map = string(insn_map, "}")
    if insn_map == "{}"
        insn_map = "{:}"
    end
    return return @sprintf("%s -> %s", params, insn_map)
end


"""
helpers for extracting accesses from rhs of an expression expr
"""
function rhs_dependencies(expr::Expr)::Vector{Union{Expr, Symbol}}
    deps = []
    if expr.head == :call
        # RHS is operation on multiple terms, add all parts
        for arg in expr.args[2:end]
            append!(deps, rhs_dependencies(arg))
        end
    elseif expr.head == :macrocall
        # RHS is macro call on multiple terms, add all parts except macro and line number node
        for arg in expr.args[3:end]
            append!(deps, rhs_dependencies(arg))
        end
    else
        push!(deps, expr)
    end
    return deps
end

function rhs_dependencies(sym::Symbol)::Vector{Union{Expr, Symbol}}
    return [sym]
end

function rhs_dependencies(num::Number)::Vector{Union{Expr, Symbol}}
    return []
end

"""
construct the access relations and dependencies map in ISL syntax form
ex: (id = :mult) C[i, j] += A[i, k] * B[k, j] becomes
{mult[i, j, k] -> A[k, i];
mult[i, j, k] -> B[j, k]}
"""
function dependencies_isl_rep(kernel::LoopKernel)::String
    params = string(kernel.args)
    params = replace.(params, r":" => "")
    params = replace.(params, r"Symbol" => "")

    deps_map = "{"
    count = 1
    for instruction in kernel.instructions
        dcount = 1
        ds = "["
        for iname in instruction.dependencies
            for domain in kernel.domains
                if domain.iname == iname
                    if dcount != 1
                        ds = string(ds, ", ")
                    end
                    dcount += 1
                    ds = string(ds, domain.iname)
                end
            end
        end
        ds = string(ds, "] ")
        dcount = 1
        for iname in instruction.dependencies
            for instruction2 in kernel.instructions
                if instruction2.iname == iname
                    dcount = 1
                    ds2 = "["
                    for iname2 in instruction2.dependencies
                        for domain in kernel.domains
                            if domain.iname == iname2
                                if dcount != 1
                                    ds2 = string(ds2, ", ")
                                end
                                dcount += 1
                                ds2 = string(ds2, domain.iname)
                            end
                        end
                    end
                    ds2 = string(ds2, "]")
                    if count != 1
                        deps_map = string(deps_map, "; ")
                    end
                    count += 1
                    deps_map = string(deps_map, sym_to_str(instruction.iname), ds, " -> ")
                    # get LHS of instruction
                    deps_map = string(deps_map, sym_to_str(instruction2.iname), ds2)
                end
            end
        end
        rhs = instruction.body.args[2]
        rhs_parts = rhs_dependencies(rhs)
        for part in rhs_parts
            if count != 1
                deps_map = string(deps_map, "; ")
            end
            count += 1
            deps_map = string(deps_map, sym_to_str(instruction.iname), ds, " -> ")
            # get LHS of instruction
            deps_map = string(deps_map, part)
        end
    end

    deps_map = string(deps_map, "}")
    if deps_map == "{}"
        deps_map = "{:}"
    end
    return @sprintf("%s -> %s", params, deps_map)
end
