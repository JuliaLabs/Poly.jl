using JuLoop
using ISL
import ISL.API
using Printf

export run_polyhedral_model


"""
Run polyhedral model on a kernel. Does the following:
creates a context and space for ISL
create an ISL union set representing the domains and instructions in kernel
create an ISL union maps representing access relations in kernel
create an ISL schedule from intial implementation
create dependencies maps by analyzing the accesses and schedule
add loop optimal ordering to dependencies
create a new schedule from all dependencies
create an AST from the schedule
parse AST to julia code
return expression

if debug is set, ISL errors are printed to stderror
if verbose is
    0: no printing (unless debug is on, where ISL errors are printed)
    1: final Julia code is printed
    2: initial C code, loop orderings, new schedule, final C code, and final Julia code printed
    3: all of 2 plus domain, access relations, dependence relations, original schedule, and new schedule constraints printed
"""
function run_polyhedral_model(kernel::LoopKernel; debug=false, verbose=0)::Expr
    """
    initialize context and set error printing options to warn/silence
    """
    context = ISL.API.isl_ctx_alloc()
    if debug
        ISL.API.isl_options_set_on_error(context, 0) # warn errors
    else
        ISL.API.isl_options_set_on_error(context, 1) # silence errors
    end

    """
    get the domain (instructions) of the kernel
    """
    instructions = ISL.API.isl_union_set_read_from_str(context, instructions_isl_rep(kernel))
    if verbose >= 3
        println("===DOMAIN===")
        ISL.API.isl_union_set_dump(instructions)
    end

    """
    construct the access relations
    """
    may_read, may_write, must_write = accesses_isl_rep(kernel)
    if verbose >= 3
        println("===ACCESS RELATIONS===")
        @show may_read
        @show may_write
        @show must_write
    end
    may_read = ISL.API.isl_union_map_read_from_str(context, may_read)
    may_write = ISL.API.isl_union_map_read_from_str(context, may_write)
    must_write = ISL.API.isl_union_map_read_from_str(context, must_write)

    """
    construct the original schedule
    """
    schedule = schedule_tree_isl_rep(context, kernel.domains[1], kernel)
    if verbose >= 3
        println("===ORIGINAL SCHEDULE===")
        ISL.API.isl_schedule_dump(schedule)
    end
    if verbose >= 2
        # show original code in C representation
        space = ISL.API.isl_set_read_from_str(context, "{:}")
        build = ISL.API.isl_ast_build_from_context(space)
        ast = ISL.API.isl_ast_build_node_from_schedule(build, ISL.API.isl_schedule_copy(schedule))
        c = ISL.API.isl_ast_node_get_ctx(ast)
        p = ISL.API.isl_printer_to_str(c)
        p = ISL.API.isl_printer_set_output_format(p, 4) # 4 = C code
        q = ISL.API.isl_printer_print_ast_node(p, ast)
        str = ISL.API.isl_printer_get_str(q)
        str = Base.unsafe_convert(Ptr{Cchar}, str)
        println("===ORIGINAL C CODE===")
        println(Base.unsafe_string(str))
    end

    """
    dependency analysis
    """
    # read after write deps
    access = ISL.API.isl_union_access_info_from_sink(ISL.API.isl_union_map_copy(may_read))
    access = ISL.API.isl_union_access_info_set_must_source(access, ISL.API.isl_union_map_copy(must_write))
    access = ISL.API.isl_union_access_info_set_schedule(access, ISL.API.isl_schedule_copy(schedule))
    flow = ISL.API.isl_union_access_info_compute_flow(access)
    raw_deps = ISL.API.isl_union_flow_get_must_dependence(flow)
    # write after read and write after write deps
    access = ISL.API.isl_union_access_info_from_sink(ISL.API.isl_union_map_copy(must_write))
    access = ISL.API.isl_union_access_info_set_must_source(access, must_write)
    access = ISL.API.isl_union_access_info_set_may_source(access, may_read)
    access = ISL.API.isl_union_access_info_set_schedule(access, schedule)
    flow = ISL.API.isl_union_access_info_compute_flow(access)
    waw_deps = ISL.API.isl_union_flow_get_must_dependence(flow)
    war_deps = ISL.API.isl_union_flow_get_may_dependence(flow)
    if verbose >= 3
        println("===DEPENDENCIES (RAW, WAW, WAR)===")
        ISL.API.isl_union_map_dump(raw_deps)
        ISL.API.isl_union_map_dump(waw_deps)
        ISL.API.isl_union_map_dump(war_deps)
    end

    """
    add loop ordering dependencies based on striding analysis
    """
    loop_ordering_deps = get_loop_ordering_deps_isl(kernel)
    if verbose >= 2
        println("===LOOP ORDERINGS===")
        @show loop_ordering_deps
    end
    loop_ordering_deps = ISL.API.isl_union_map_read_from_str(context, loop_ordering_deps)

    """
    construct new schedule constraints from dependencies
    """
    # dependence constraints
    all_deps = ISL.API.isl_union_map_union(waw_deps, war_deps)
    all_deps = ISL.API.isl_union_map_union(all_deps, raw_deps)
    loop_ordering_deps = ISL.API.isl_union_map_union(ISL.API.isl_union_map_copy(all_deps), loop_ordering_deps)
    schedule_constraints = ISL.API.isl_schedule_constraints_on_domain(ISL.API.isl_union_set_copy(instructions))
    schedule_constraints = ISL.API.isl_schedule_constraints_set_validity(schedule_constraints, ISL.API.isl_union_map_copy(loop_ordering_deps))

    # proximity constraints (keeps loops nested based on dependencies)
    schedule_constraints = ISL.API.isl_schedule_constraints_set_proximity(schedule_constraints, ISL.API.isl_union_map_copy(all_deps))

    # coincidence constraints (things to be scheduled together- doesn't usually do much)
    # schedule_constraints = ISL.API.isl_schedule_constraints_set_coincidence(schedule_constraints, ISL.API.isl_union_map_copy(all_deps))

    if verbose >= 3
        println("===SCHEDULE CONSTRAINTS===")
        ISL.API.isl_schedule_constraints_dump(schedule_constraints)
    end

    """
    scheduling options
    """
    # maximize band size- better for tiling
    ISL.API.isl_options_set_schedule_maximize_band_depth(context, 1)
    # schedule weakly-connected components together
    ISL.API.isl_options_set_schedule_whole_component(context, 1)
    ISL.API.isl_options_set_schedule_treat_coalescing(context, 1)

    """
    compute the new schedule
    """
    schedule = ISL.API.isl_schedule_constraints_compute_schedule(schedule_constraints)

    # handle scheduling error (usually when loop orderings are illegal, try without)
    if ISL.API.isl_ctx_last_error_line(context) != -1
        # reset context error
        ISL.API.isl_ctx_reset_error(context)
        # try schedule without loop reordering
        schedule_constraints = ISL.API.isl_schedule_constraints_on_domain(instructions)
        schedule_constraints = ISL.API.isl_schedule_constraints_set_validity(schedule_constraints, ISL.API.isl_union_map_copy(all_deps))
        schedule_constraints = ISL.API.isl_schedule_constraints_set_proximity(schedule_constraints, ISL.API.isl_union_map_copy(all_deps))
        schedule = ISL.API.isl_schedule_constraints_compute_schedule(schedule_constraints)
    end

    # other modifications to schedule tree

    if verbose >= 2
        println("===NEW SCHEDULE===")
        ISL.API.isl_schedule_dump(schedule)
    end

    """
    construct AST from new schedule
    """
    space = ISL.API.isl_set_read_from_str(context, "{:}")
    build = ISL.API.isl_ast_build_from_context(space)
    ast = ISL.API.isl_ast_build_node_from_schedule(build, schedule)

    if verbose >= 2
        # print C code
        c = ISL.API.isl_ast_node_get_ctx(ast)
        p = ISL.API.isl_printer_to_str(c)
        p = ISL.API.isl_printer_set_output_format(p, 4) # 4 = C code
        q = ISL.API.isl_printer_print_ast_node(p, ast)
        str = ISL.API.isl_printer_get_str(q)
        str = Base.unsafe_convert(Ptr{Cchar}, str)
        println("===FINAL C CODE===")
        println(Base.unsafe_string(str))
    end

    """
    parse back to Julia code
    """
    expr = parse_ast(ast, kernel)
    if verbose >= 1
        println("===FINAL JULIA CODE===")
        @show expr
    end
    return expr
end


"""
POLYHEDRAL MODEL GENERATION FUNCTIONS
Julia -> ISL
"""

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
helper to return a unique identifier for a kernel existential quantifier
"""
function unique_identifier()::String
    id = gensym("id")
    return sym_to_str(id)
end


"""
helper to get string rep of params in kernel
"""
function get_params_str(kernel::LoopKernel)
    params = string(kernel.consts)
    params = replace.(params, r":" => "")
    params = replace.(params, r"Symbol" => "")
    return params
end


"""
format a condition on an instruction properly
"""
function format_instruction_cond(cond::Expr)::String
    stringcond = string(cond)
    stringcond = replace.(stringcond, r"%" => " mod ")
    stringcond = replace.(stringcond, r"==" => " = ")
    stringcond = replace.(stringcond, r"&&" => " and ")
    stringcond = replace.(stringcond, r"&" => " and ")
    stringcond = replace.(stringcond, "||" => " or ")
    stringcond = replace.(stringcond, "|" => " or ")
    return stringcond
end


"""
helper to get the domains related to instruction(s) in the ISL format
ex: 0 <= i <= n and 0 <= j <= n
"""
function get_instructions_domains(instructions::Vector{Instruction}, kernel::LoopKernel)::String

    conditions = ""
    count = 1
    for domain in kernel.domains
        for instruction in instructions
            if domain.iname in instruction.dependencies
                if count != 1
                    conditions = string(conditions, " and ")
                end
                count += 1
                instructioncond = ""
                if instruction.cond != :()
                    # instruction has separate domain condition
                    instructioncond = string(" and ", format_instruction_cond(instruction.cond))
                end
                if domain.step != 1
                    identifier = unique_identifier()
                    conditions = string(conditions, @sprintf("exists %s: %s = %s*%s + %s and %s <= %s <= %s", identifier, string(domain.iname), string(domain.step), identifier, string(domain.lowerbound), string(domain.lowerbound), string(domain.iname), string(domain.upperbound)), instructioncond)
                else
                    conditions = string(conditions, @sprintf("%s <= %s <= %s", string(domain.lowerbound), string(domain.iname), string(domain.upperbound)), instructioncond)
                end
            end
        end
    end
    return conditions
end


"""
get the related iteration spaces for an instruction
ex: [i, j]
if primes is true, returns [i', j']
"""
function get_related_domains(instruction::Instruction, kernel::LoopKernel, primes=false)
    count = 1
    domains = "["
    for domain in kernel.domains
        for iname in instruction.dependencies
            if domain.iname == iname
                if count != 1
                    domains = string(domains, ", ")
                end
                count += 1
                domains = string(domains, domain.iname)
                if primes
                    domains = string(domains, "\'")
                end
            end
        end
    end
    domains = string(domains, "]")
    return domains
end


"""
construct the instructions' domains representation in ISL syntax form
ex: (id = :mult) C[i, j] += A[i, k] * B[k, j] in matrix multiplication becomes
[n, m, r, A, B, C] -> {mult[i, j, k] : 0 <= i <= n and 0 <= j <= m and 0 <= k <= r}
"""
function instructions_isl_rep(kernel::LoopKernel)::String
    params = get_params_str(kernel)

    insn_map = "{"
    icount = 1
    for instruction in kernel.instructions
        count = 1
        name = string(sym_to_str(instruction.iname), get_related_domains(instruction, kernel))
        conditions = get_instructions_domains([instruction], kernel)

        if icount != 1
            insn_map = string(insn_map, "; ")
        end
        icount += 1

        if conditions != ""
            insn_map = string(insn_map, name, " : ", conditions)
        else
            insn_map = string(insn_map, name)
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
function expr_accesses(expr::Expr, kernel::LoopKernel)::Vector{Union{Expr, Symbol}}
    deps = []
    if expr.head == :call
        # RHS is operation on multiple terms, add all parts
        for arg in expr.args[2:end]
            append!(deps, expr_accesses(arg, kernel))
        end
    elseif expr.head == :macrocall
        # RHS is macro call on multiple terms, add all parts except macro and line number node
        for arg in expr.args[3:end]
            append!(deps, expr_accesses(arg, kernel))
        end
    else
        push!(deps, expr)
    end
    return deps
end

function expr_accesses(sym::Symbol, kernel::LoopKernel)::Vector{Union{Expr, Symbol}}
    if sym in kernel.consts
        return []
    end
    if sym in [d.iname for d in kernel.domains]
        return []
    end
    return [:($sym[])]
end

function expr_accesses(num::Number, kernel::LoopKernel)::Vector{Union{Expr, Symbol}}
    return []
end

function expr_accesses(l::LineNumberNode, kernel::LoopKernel)::Vector{Union{Expr, Symbol}}
    return []
end

"""
construct the access relations map in ISL syntax form

ex: (id = :mult) C[i, j] += A[i, k] * B[k, j] becomes
may_read:
[n, A, B, C] -> {
mult[i, j, k] -> A[k, i] : 0 <= i <= n and 0 <= j <= n and 0 <= k <= n;
mult[i, j, k] -> B[j, k] : 0 <= i <= n and 0 <= j <= n and 0 <= k <= n}
may_write:
[n, A, B, C] -> {mult[i, j, k] -> C[i, j] : 0 <= i <= n and 0 <= j <= n and 0 <= k <= n}
must_write:
[n, A, B, C] -> {mult[i, j, k] -> C[i, j] : 0 <= i <= n and 0 <= j <= n and 0 <= k <= n}
"""
function accesses_isl_rep(kernel::LoopKernel)::Tuple{String, String, String}
    params = get_params_str(kernel)

    may_read = string(params, " -> {")
    may_write = string(params, " -> {")
    must_write = string(params, " -> {")
    rcount = 1
    wcount = 1
    mcount = 1
    for instruction in kernel.instructions
        ds = get_related_domains(instruction, kernel)
        reads = expr_accesses(instruction.body.args[2], kernel)
        writes = expr_accesses(instruction.body.args[1], kernel)
        head = instruction.body.head
        if head != :(=)
            # read from LHS also
            append!(reads, writes)
        end
        # must writes in LHS
        for part in writes
            if mcount != 1
                must_write = string(must_write, "; ")
            end
            mcount += 1
            must_write = string(must_write, sym_to_str(instruction.iname), ds, " -> ", part)
            # get domains
            conditions = get_instructions_domains([instruction], kernel)
            if conditions != ""
                must_write = string(must_write, " : ", conditions)
            end
        end
        # may reads in RHS and LHS
        for part in reads
            if rcount != 1
                may_read = string(may_read, "; ")
            end
            rcount += 1
            may_read = string(may_read, sym_to_str(instruction.iname), ds, " -> ", part)
            # get domains
            conditions = get_instructions_domains([instruction], kernel)
            if conditions != ""
                may_read = string(may_read, " : ", conditions)
            end
        end
    end

    may_read = string(may_read, "}")
    if may_read == "{}"
        may_read = "{:}"
    end
    must_write = string(must_write, "}")
    if must_write == "{}"
        must_write = "{:}"
    end
    may_write = must_write
    return may_read, may_write, must_write
end


"""
Get the loop ordering dependencies for a kernel
For example, to specify the order j, k, i in a loop nest:
    "[n, r, m] -> {out253[i, j, k] -> out253[i', j', k']: j, k, i << j', k', i'}"
"""
function get_loop_ordering_deps_isl(kernel::LoopKernel)::String
    params = get_params_str(kernel)

    ordering_deps = "{"
    count = 0

    loop_orderings = get_best_nesting_orders(kernel)

    for order in loop_orderings
        for domain in kernel.domains
            if domain.iname in order
                for instruction in domain.instructions
                    if typeof(instruction) == Instruction
                        # found instruction to use
                        if count != 0
                            ordering_deps = string(ordering_deps, "; ")
                        end
                        count += 1
                        # sym[i, j] -> sym[i', j']
                        ordering_deps = string(ordering_deps, sym_to_str(instruction.iname), get_related_domains(instruction, kernel), " -> ", sym_to_str(instruction.iname), get_related_domains(instruction, kernel, true), ": ")
                        order_string = ""
                        order_string_prime = ""
                        ocount = 0
                        for iname in order
                            if iname in instruction.dependencies
                                if ocount != 0
                                    order_string = string(order_string, ", ")
                                    order_string_prime = string(order_string_prime, ", ")
                                end
                                ocount += 1
                                order_string = string(order_string, iname)
                                order_string_prime = string(order_string_prime, iname, "\'")
                            end
                        end
                        # i, j << i', j'
                        ordering_deps = string(ordering_deps, order_string, " << ", order_string_prime)
                        break
                    end
                end
            end
        end
    end

    ordering_deps = string(ordering_deps, "}")
    return @sprintf("%s -> %s", params, ordering_deps)
end


"""
construct a partial schedule for a grouping of instructions in domains
"""
function partial_schedule_isl(domains::Vector{Domain}, instructions::Vector{Instruction}, kernel::LoopKernel)::String

    mupa = "["

    dcount = 0
    for domain in domains
        set = ""

        icount = 0
        for instruction in instructions
            if domain.iname in instruction.dependencies
                if icount != 0
                    set = string(set, "; ")
                end
                icount += 1
                set = string(set, sym_to_str(instruction.iname), get_related_domains(instruction, kernel), "-> [(", domain.iname, ")]")
            end
        end

        if set != ""
            if dcount != 0
                mupa = string(mupa, ", ")
            end
            dcount += 1
            mupa = string(mupa, "{", set, "}")
        end
    end

    mupa = string(mupa, "]")

    return mupa
end


"""
construct an original schedule for dependency analysis
"""
function schedule_tree_isl_rep(context::Ptr{ISL.API.isl_ctx}, domain::Domain, kernel::LoopKernel, parent_doms=[])::Ptr{ISL.API.isl_schedule}

    new_parent_doms = [domain]
    append!(new_parent_doms, parent_doms)
    schedule = :nothing
    # navigate through original scheduling
    for instruction in domain.instructions
        if typeof(instruction) == Instruction
            # construct instruction schedule
            temp_kern = LoopKernel([instruction], kernel.domains, kernel.args, kernel.argtypes, kernel.consts)
            inst_dom = instructions_isl_rep(temp_kern)
            inst_dom = ISL.API.isl_union_set_read_from_str(context, inst_dom)
            inst_schedule = ISL.API.isl_schedule_from_domain(inst_dom)

            # put in sequence
            if schedule == :nothing
                schedule = inst_schedule
            else
                schedule = ISL.API.isl_schedule_sequence(schedule, inst_schedule)
            end
        else # Domain
            dom_schedule = schedule_tree_isl_rep(context, instruction, kernel, new_parent_doms)
            if schedule == :nothing
                schedule = dom_schedule
            else
                schedule = ISL.API.isl_schedule_sequence(schedule, dom_schedule)
            end
        end
    end

    # insert partial schedule
    insts = Instruction[]
    for instruction in kernel.instructions
        if domain.iname in instruction.dependencies
            push!(insts, instruction)
        end
    end
    partial = partial_schedule_isl([domain], insts, kernel)
    if partial != "[]"
        partial = ISL.API.isl_multi_union_pw_aff_read_from_str(context, partial)
        schedule = ISL.API.isl_schedule_insert_partial_schedule(schedule, partial)
    end

    if schedule == :nothing
        error("Mis-defined initial schedule! Check that instructions are properly nested")
    end

    return schedule
end


"""
PARSING FUNCTIONS
ISL -> Julia
"""

"""
parse an isl AST back into julia code
"""
function parse_ast(ast::Ptr{ISL.API.isl_ast_node}, kernel::LoopKernel)::Expr
    node_type = ISL.API.isl_ast_node_get_type(ast)

    if node_type == ISL.API.isl_ast_node_for
        # for loop
        expr = parse_ast_for(ast, kernel)

    elseif node_type == ISL.API.isl_ast_node_if
        # if statement
        expr = parse_ast_if(ast, kernel)

    elseif node_type == ISL.API.isl_ast_node_block
        # compound node- multiple elements in block
        expr = parse_ast_block(ast, kernel)

    elseif node_type == ISL.API.isl_ast_node_mark
        # mark in the ast- comment

    elseif node_type == ISL.API.isl_ast_node_user
        # expression statement
        expr = parse_ast_user(ast, kernel)

    elseif node_type == ISL.API.isl_ast_node_error
        # error- quit

    end
end


"""
parse an ast node representing a for loop into an expression
"""
function parse_ast_for(ast::Ptr{ISL.API.isl_ast_node}, kernel::LoopKernel)::Expr
    executes_once = ISL.API.isl_ast_node_for_is_degenerate(ast) # if loop executes only once

    # loop qualities
    iter = ISL.API.isl_ast_node_for_get_iterator(ast)
    id = ISL.API.isl_ast_expr_get_id(iter) # id of iterator
    ISL.API.isl_id_free(id)
    name = Base.unsafe_convert(Ptr{Cchar}, ISL.API.isl_id_get_name(id))
    name = Symbol(Base.unsafe_string(name)) # iterator symbol
    lb = parse_ast_expr(ISL.API.isl_ast_node_for_get_init(ast)) # initial condition
    body = parse_ast(ISL.API.isl_ast_node_for_get_body(ast), kernel) # body of loop

    if executes_once == ISL.API.isl_bool_true
        # loop only runs once, just set name = init and run body
        expr = quote
            $name = $init
            $body
        end
        return expr
    else
        step = parse_ast_expr(ISL.API.isl_ast_node_for_get_inc(ast)) # incremental step
        cond = parse_ast_expr(ISL.API.isl_ast_node_for_get_cond(ast)) # final condition
        ub = cond.args[3]
        if cond.args[1] == :(<)
            ub = :($ub + 1)
        elseif cond.args[1] != :(<=)
            error("unexpected conditional operator ", cond.head)
        end

        if step == 1
            # for name = lb:up
                # body
            expr = quote
                for $name = $lb:$ub
                    $body
                end
            end
            return expr

        else
            # for name = lb:step:up
                # body
            expr = quote
                for $name = $lb:$step:$ub
                    $body
                end
            end
            return expr
        end
    end
end


"""
parse an ast node representing an if block into an expression
"""
function parse_ast_if(ast::Ptr{ISL.API.isl_ast_node}, kernel::LoopKernel)::Expr
    cond = ISL.API.isl_ast_node_if_get_cond(ast)
    cond = parse_ast_expr(cond)
    if_body = ISL.API.isl_ast_node_if_get_then(ast)
    if_body = parse_ast(if_body, kernel)
    has_else = ISL.API.isl_ast_node_if_has_else(ast)
    if has_else == ISL.API.isl_bool_true
        # if cond
            # if_body
        # else
            # else_body
        # end
        else_body = ISL.API.isl_ast_node_if_get_else(ast)
        else_body = parse_ast(else_body, kernel)
        expr = quote
            if $cond
                $if_body
            else
                $else_body
            end
        end
        return expr
    else
        # if cond
            # if_body
        # end
        expr = quote
            if $cond
                $if_body
            end
        end
        return expr
    end
end


"""
parse an ast node representing a block into an expression
"""
function parse_ast_block(ast::Ptr{ISL.API.isl_ast_node}, kernel::LoopKernel)::Expr
    children = ISL.API.isl_ast_node_block_get_children(ast)
    n = ISL.API.isl_ast_node_list_n_ast_node(children)
    exs = []
    for i=0:n-1
        child = ISL.API.isl_ast_node_list_get_at(children, i)
        push!(exs, parse_ast(child, kernel))
    end
    return quote $(exs...) end
end

"""
helpers for replacing symbols in expression
"""
function replace_expr_syms(ex::Expr, sym_dict::Dict{Symbol, Union{Symbol, Number, Expr}})::Expr
    new_args = []
    for arg in ex.args
        push!(new_args, replace_expr_syms(arg, sym_dict))
    end
    new_ex = Expr(ex.head, new_args...)
    return new_ex
end

function replace_expr_syms(sym::Symbol, sym_dict::Dict{Symbol, Union{Symbol, Number, Expr}})
    if haskey(sym_dict, sym)
        return sym_dict[sym]
    end
    return sym
end

function replace_expr_syms(num::Number, sym_dict::Dict{Symbol, Union{Symbol, Number, Expr}})
    return num
end

function replace_expr_syms(l::LineNumberNode, sym_dict::Dict{Symbol, Union{Symbol, Number, Expr}})
    return l
end


"""
parse an ast node representing an expression
"""
function parse_ast_user(ast::Ptr{ISL.API.isl_ast_node}, kernel::LoopKernel)::Expr
    expr = ISL.API.isl_ast_node_user_get_expr(ast)
    ISL.API.isl_ast_expr_free(expr)
    first_expr = ISL.API.isl_ast_expr_op_get_arg(expr, 0) # first arg is "function" name (symbol cooresponding to instruction)
    n = ISL.API.isl_ast_expr_get_op_n_arg(expr)
    id = ISL.API.isl_ast_expr_id_get_id(first_expr)
    ISL.API.isl_ast_expr_free(first_expr)
    name = Base.unsafe_convert(Ptr{Cchar}, ISL.API.isl_id_get_name(id)) # name of identifier
    name = Base.unsafe_string(name)
    ISL.API.isl_id_free(id)

    # search for matching instruction
    for instruction in kernel.instructions
        if sym_to_str(instruction.iname) == name
            # found instruction
            # replace original iterators with new iterators
            old_ids = []
            for domain in kernel.domains
                for iname in instruction.dependencies
                    if domain.iname == iname
                        push!(old_ids, iname)
                    end
                end
            end
            new_ids = Dict{Symbol, Union{Symbol, Number, Expr}}()
            for i=1:Int(n)-1
                arg = ISL.API.isl_ast_expr_op_get_arg(expr, i)
                arg = parse_ast_expr(arg)
                new_ids[old_ids[i]] = arg
            end
            ex = replace_expr_syms(instruction.body, new_ids)
            return :(@inbounds $ex)
        end
    end

    return :()
end


"""
parse an ast expression condition
"""
function parse_ast_expr(expr::Ptr{ISL.API.isl_ast_expr})::Union{Symbol, Expr, Number}
    type = ISL.API.isl_ast_expr_get_type(expr)
    if type == ISL.API.isl_ast_expr_id
        id = ISL.API.isl_ast_expr_id_get_id(expr)
        name = Base.unsafe_convert(Ptr{Cchar}, ISL.API.isl_id_get_name(id)) # name of identifier
        name = Base.unsafe_string(name)
        ISL.API.isl_id_free(id)
        return :($(Symbol(name)))
    elseif type == ISL.API.isl_ast_expr_int
        val = ISL.API.isl_ast_expr_int_get_val(expr)
        num = ISL.API.isl_val_get_d(val)
        if Int(num) == num
            num = Int(num)
        end
        ISL.API.isl_val_free(val)
        return :($num)
    else
        if ISL.API.isl_ast_expr_get_op_n_arg(expr) == 1
            arg1 = parse_ast_expr(ISL.API.isl_ast_expr_op_get_arg(expr, 0))
            arg2 = :()
        else
            arg1 = parse_ast_expr(ISL.API.isl_ast_expr_op_get_arg(expr, 0))
            arg2 = parse_ast_expr(ISL.API.isl_ast_expr_op_get_arg(expr, 1))
        end

        op_type = ISL.API.isl_ast_expr_op_get_type(expr)
        if op_type == ISL.API.isl_ast_expr_op_and || op_type == ISL.API.isl_ast_expr_op_and_then
            return :($arg1 && $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_or || op_type == ISL.API.isl_ast_expr_op_or_else
            return :($arg1 || $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_max
            return :(max($arg1, $arg2))
        elseif op_type == ISL.API.isl_ast_expr_op_min
            return :(min($arg1, $arg2))
        elseif op_type == ISL.API.isl_ast_expr_op_minus
            return :(-$arg1)
        elseif op_type == ISL.API.isl_ast_expr_op_add
            return :($arg1 + $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_sub
            return :($arg1 - $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_mul
            return :($arg1 * $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_fdiv_q
            return :(Int(floor($arg1 / $arg2)))
        elseif op_type == ISL.API.isl_ast_expr_op_pdiv_q || op_type == ISL.API.isl_ast_expr_op_div
            return :(Int($arg1 / $arg2))
        elseif op_type == ISL.API.isl_ast_expr_op_pdiv_r || op_type == ISL.API.isl_ast_expr_op_zdiv_r
            return :($arg1 % $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_eq
            return :($arg1 == $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_le
            return :($arg1 <= $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_ge
            return :($arg1 >= $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_lt
            return :($arg1 < $arg2)
        elseif op_type == ISL.API.isl_ast_expr_op_gt
            return :($arg1 > $arg2)
        else
            # invalid
            return :()
        end
    end
end
