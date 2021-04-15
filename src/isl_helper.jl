using JuLoop
using ISL
import ISL.API
using Printf

export run_polyhedral_model


"""
Run polyhedral model on a kernel. Does the following:
creates a context and space for ISL
create an ISL union set representing the domains and instructions in kernel
create an ISL union map representing dependencies and access relations in kernel
create a schedule map by intersecting the above
create an AST from the scheule map
save C code to out.txt representing the AST
"""
function run_polyhedral_model(kernel::LoopKernel)#::Tuple{ISL.API.isl_basic_set, ISL.API.isl_ctx}
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
    ISL.API.isl_ast_node_dump(ast)
    println("DUMPED NODE")
    c = ISL.API.isl_ast_node_get_ctx(ast)
    # p = ISL.API.isl_printer_to_str(c)
    file = ccall((:fopen,), Ptr{Libc.FILE}, (Ptr{Cchar}, Ptr{Cchar}), "out.txt", "w+")
    # file = fopen("/tmp/test.txt", "w+")
    p = ISL.API.isl_printer_to_file(c, file)
    p = ISL.API.isl_printer_set_output_format(p, 4) # 4 = C code
    q = ISL.API.isl_printer_print_ast_node(p, ast)
    p = ISL.API.isl_printer_flush(p)
    # s = ISL.API.isl_printer_get_str(q)
    # println(Base.unsafe_load(s))
    println("PRINTED CODE TO out.txt")
end


"""
helper to turn gensyms into strings without special characters
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
construct the instructions' domains representation in ISL syntax form
ex: (id = :mult) C[i, j] += A[i, k] * B[k, j] in matrix multiplication becomes
[n, m, r, A, B, C] -> {mult[i, j, k] : 0 <= i <= n and 0 <= j <= m and 0 <= k <= r}
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
    return @sprintf("%s -> %s", params, insn_map)
end


"""
helpers for extracting accesses from an expression expr
"""
function expr_accesses(expr::Expr)::Vector{Union{Expr, Symbol}}
    deps = []
    if expr.head == :call
        # RHS is operation on multiple terms, add all parts
        for arg in expr.args[2:end]
            append!(deps, expr_accesses(arg))
        end
    elseif expr.head == :macrocall
        # RHS is macro call on multiple terms, add all parts except macro and line number node
        for arg in expr.args[3:end]
            append!(deps, expr_accesses(arg))
        end
    else
        push!(deps, expr)
    end
    return deps
end

function expr_accesses(sym::Symbol)::Vector{Union{Expr, Symbol}}
    return [sym]
end

function expr_accesses(num::Number)::Vector{Union{Expr, Symbol}}
    return [:(num[])]
end

"""
construct the access relations and dependencies map in ISL syntax form

ex: (id = :mult) C[i, j] += A[i, k] * B[k, j] becomes
[n, A, B, C] -> {mult[i, j, k] -> A[k, i] : 0 <= i <= n and 0 <= j <= n and 0 <= k <= n;
mult[i, j, k] -> B[j, k] : 0 <= i <= n and 0 <= j <= n and 0 <= k <= n;
mult[i, j, k] -> C[i, j] : 0 <= i <= n and 0 <= j <= n and 0 <= k <= n}

Any dependencies between instructions are also added, such as:
{mult[i, j, k] -> earlierinst[i, j] : 0 <= i <= n and 0 <= j <= n and 0 <= k <= n}
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
                    deps_map = string(deps_map, sym_to_str(instruction.iname), ds, " -> ", sym_to_str(instruction2.iname), ds2)
                    # get domains
                    conditions = ""
                    ccount = 1
                    for domain in kernel.domains
                        if domain.iname in instruction.dependencies || domain.iname in instruction2.dependencies
                            step = domain.recurrence.args[2]
                            if ccount != 1
                                conditions = string(conditions, " and ")
                            end
                            ccount += 1
                            if step != 1
                                conditions = string(conditions, @sprintf("exists (a: %s = %da and %s <= %s <= %s)", string(domain.iname), step, string(domain.lowerbound), string(domain.iname), string(domain.upperbound)))
                            else
                                conditions = string(conditions, @sprintf("%s <= %s <= %s", string(domain.lowerbound), string(domain.iname), string(domain.upperbound)))
                            end
                        end
                    end
                    if conditions != ""
                        deps_map = string(deps_map, " : ", conditions)
                    end
                end
            end
        end
        rhs_parts = expr_accesses(instruction.body.args[2])
        lhs_parts = expr_accesses(instruction.body.args[1])
        for part in append!(rhs_parts, lhs_parts)
            if count != 1
                deps_map = string(deps_map, "; ")
            end
            count += 1
            deps_map = string(deps_map, sym_to_str(instruction.iname), ds, " -> ", part)
            # get domains
            conditions = ""
            ccount = 1
            for domain in kernel.domains
                if domain.iname in instruction.dependencies
                    step = domain.recurrence.args[2]
                    if ccount != 1
                        conditions = string(conditions, " and ")
                    end
                    ccount += 1
                    if step != 1
                        conditions = string(conditions, @sprintf("exists (a: %s = %da and %s <= %s <= %s)", string(domain.iname), step, string(domain.lowerbound), string(domain.iname), string(domain.upperbound)))
                    else
                        conditions = string(conditions, @sprintf("%s <= %s <= %s", string(domain.lowerbound), string(domain.iname), string(domain.upperbound)))
                    end
                end
            end
            if conditions != ""
                deps_map = string(deps_map, " : ", conditions)
            end
        end
    end

    deps_map = string(deps_map, "}")
    if deps_map == "{}"
        deps_map = "{:}"
    end
    return @sprintf("%s -> %s", params, deps_map)
end
