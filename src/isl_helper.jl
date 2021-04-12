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
    @show domains_isl_rep(kernel)
    space = ISL.API.isl_set_read_from_str(context, domains_isl_rep(kernel))
    println("GOT SPACE")
    @show instructions_isl_rep(kernel)
    instructions = ISL.API.isl_set_read_from_str(context, instructions_isl_rep(kernel))
    println("GOT INSTS")
    @show dependencies_isl_rep(kernel)
    dependencies = ISL.API.isl_map_read_from_str(context, dependencies_isl_rep(kernel))
    println("GOT DEPS")
    schedule_map = ISL.API.isl_map_intersect_domain(dependencies, instructions)
    println("GOT MAP")
    schedule_map = ISL.API.isl_union_map_from_map(schedule_map)
    println("GOT MAP UNION")
    build = ISL.API.isl_ast_build_from_context(space)
    println("GOT BUILD")
    ast = ISL.API.isl_ast_build_node_from_schedule_map(build, schedule_map)
    println("GOT AST")
    print(ISL.API.isl_ast_node_dump(ast))
end


"""
modify loop domains to be in the ISL syntax form
ex: i = 1:10 becomes {[i]: 1 <= i <= 10}
ex: i = 1:2:10 becomes  {[i]: exists a: i = 2a and 1 <= i <= 10}
"""
function domains_isl_rep(kernel::LoopKernel)::String
    keys = string([domain.iname for domain in kernel.domains])
    keys = replace.(keys, r":" => "")

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

    return @sprintf("{%s: %s}", keys, conditions)
end


"""
construct the instructions representation in ISL syntax form
ex: C[i, j] = A[i, k] * B[k, j] in matrix multiplication becomes
{C[i, j] : 0 <= i <= n and 0 <= j <= m;
A[i, k]: 0 <= i <= n and 0 <= k <= r;
B[k, j]: 0 <= k <= r and 0 <= j <= m}
"""
function instructions_isl_rep(kernel::LoopKernel)::String
    insn_map = "{"
    icount = 1
    for instruction in kernel.instructions
        if icount != 1
            insn_map = string(insn_map, ", ")
        end
        icount += 1
        count = 1
        insn_map = string(insn_map, instruction.body.args[1], ": ")
        conditions = ""
        in_domains = []
        for iname in instruction.dependencies
            for domain in kernel.domains
                if domain.iname == iname
                    if inexpr(instruction.body.args[1], domain.iname)
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
                    push!(in_domains, domain)
                end
            end
        end
        count = 1
        insn_map = string(insn_map, conditions)
        if icount != 1
            insn_map = string(insn_map, ", ")
        end
        icount += 1
        exp = instruction.body.args[2]
        rhs_parts = []
        if typeof(exp) == Expr && typeof(exp.args[2]) == Expr && sizeof(exp.args[2].args, 1) > 1
            # RHS is operation on multiple terms
            append!(rhs_parts, exp.args[2].args)
        elseif typeof(exp) == Expr || typeof(exp) == Symbol
            push!(rhs_parts, exp)
        end
        for part in rhs_parts
            insn_map = string(insn_map, part, ": ")
            conditions = ""
            for domain in in_domains
                if inexpr(part, domain.iname)
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
            end
            insn_map = string(insn_map, conditions)
        end
    end

    insn_map = string(insn_map, "}")
    if insn_map == "{}"
        insn_map = "{:}"
    end
    return insn_map
end


"""
construct the dependencies map in ISL syntax form
ex: :mult = C[i, j] = A[i, k] * B[k, j] becomes
{C[i, j] -> A[i, k];
C[i, j] -> B[k, j];
C[i, j] -> C[i, j]}
"""
function dependencies_isl_rep(kernel::LoopKernel)::String
    deps_map = "{"
    count = 1
    for instruction in kernel.instructions
        for iname in instruction.dependencies
            for instruction2 in kernel.instructions
                if instruction2.iname == iname
                    if count != 1
                        deps_map = string(deps_map, ", ")
                    end
                    count += 1
                    deps_map = string(deps_map, instruction.body.args[1], "-> ")
                    # get LHS of instruction
                    deps_map = string(deps_map, instruction2.body.args[1])
                end
            end
        end
        exp = instruction.body.args[2]
        rhs_parts = []
        if typeof(exp) == Expr && typeof(exp.args[2]) == Expr && sizeof(exp.args[2].args, 1) > 1
            # RHS is operation on multiple terms
            append!(rhs_parts, exp.args[2].args)
        elseif typeof(exp) == Expr || typeof(exp) == Symbol
            push!(rhs_parts, exp)
        end
        if instruction.body.head != :(=)
            push!(rhs_parts, instruction.body.args[1])
        end
        for part in rhs_parts
            if count != 1
                deps_map = string(deps_map, ", ")
            end
            count += 1
            deps_map = string(deps_map, instruction.body.args[1], "-> ")
            # get LHS of instruction
            deps_map = string(deps_map, part)
        end
    end

    deps_map = string(deps_map, "}")
    if deps_map == "{}"
        deps_map = "{:}"
    end
    return deps_map
end
