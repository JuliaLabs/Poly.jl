using MacroTools
using JuLoop

export compile

"""
DEPENDENCIES FUNCTIONS
"""

"""
modifies the dependencies of the kernel's instructions to include those on loops
"""
function add_loop_deps(kernel::LoopKernel)
    for inst in kernel.instructions
        for domain in kernel.domains
            if inexpr(inst.body, domain.iname)
                # instruction references loop iname
                push!(inst.dependencies, domain.iname)
                push!(domain.instructions, inst)
            elseif domain.iname in inst.dependencies
                # add instructions for loops declared as dependencies
                push!(domain.instructions, inst)
            end
        end
    end
    for i = 1:length(kernel.domains)
        domain1 = kernel.domains[i]
        for j = i+1:length(kernel.domains)
            domain2 = kernel.domains[j]
            if inexpr(domain1.recurrence, domain2.iname)
                push!(domain2.dependencies, domain1.iname)
                push!(domain1.instructions, domain2)
            elseif isexpr(domain1.lowerbound) && inexpr(domain1.lowerbound, domain2.iname)
                push!(domain2.dependencies, domain1.iname)
                push!(domain1.instructions, domain2)
            elseif isexpr(domain1.upperbound) && inexpr(domain1.upperbound, domain2.iname)
                push!(domain2.dependencies, domain1.iname)
                push!(domain1.instructions, domain2)
            end
        end
    end
    # TODO: add dependencies between loops due to instructions?
end


"""
uses the implicit order of kernel instructions to add dependencies between instructions
"""
function add_impl_deps(kernel::LoopKernel)
    for i = 1:length(kernel.instructions)
        # look only in instructions following i
        for j = i+1:length(kernel.instructions)
            # check if there is a read/write dependency
            # write in i, read in j -> dependency
            # read in i, write in j -> dependency (stale read in i otherwise)
            # write in i, write in j -> dependency (wrong final value otherwise)

            # first check for function call (and if so, assume dependency)
            # TODO more intelligent dependencies with function calls
            if kernel.instructions[j].body.head == :call
                push!(kernel.instructions[j].dependencies, kernel.instructions[i].iname)
            else
                # need to check if LHS symbol that is being written in i is read in j
                lhs_i = kernel.instructions[i].body.args[1]
                while typeof(lhs_i) != Symbol
                    lhs_i = lhs_i.args[1]
                end
                rhs_j = kernel.instructions[j].body
                if rhs_j.head == :(=)
                    # only need rhs, since not +=, etc
                    rhs_j = rhs_j.args[2]
                end
                if inexpr(rhs_j, lhs_i)
                    push!(kernel.instructions[j].dependencies, kernel.instructions[i].iname)
                end
                # need to check if RHS read in i is written to in j
                lhs_j = kernel.instructions[j].body.args[1]
                while typeof(lhs_j) != Symbol
                    lhs_j = lhs_j.args[1]
                end
                rhs_i = kernel.instructions[i].body
                if rhs_i.head == :(=)
                    # only need rhs, since not +=, etc
                    rhs_i = rhs_i.args[2]
                end
                if inexpr(rhs_i, lhs_j)
                    push!(kernel.instructions[j].dependencies, kernel.instructions[i].iname)
                end
                # need to check if i and j write to same thing
                if inexpr(lhs_i, lhs_j)
                    push!(kernel.instructions[j].dependencies, kernel.instructions[i].iname)
                end
            end
        end
    end
end


"""
AST FUNCTIONS
"""

"""
constructs an AST using existing dependencies in kernel
"""
function construct_ast(kernel::LoopKernel)::AST
    add_impl_deps(kernel)
    add_loop_deps(kernel)

    ast = AST(nothing, AST[], AST[]) # head node
    remaining_items = Union{Instruction, Domain}[]
    nodes = Dict{Symbol, AST}()

    # add base instructions
    for instruction in kernel.instructions
        if isempty(instruction.dependencies)
            node = AST(instruction, AST[], AST[])
            push!(ast.children, node)
            push!(node.parents, ast)
            nodes[instruction.iname] = node
        else
            push!(remaining_items, instruction)
        end

    end

    # add base domains
    for domain in kernel.domains
        if isempty(domain.dependencies)
            node = AST(domain, AST[], AST[])
            push!(ast.children, node)
            push!(node.parents, ast)
            nodes[domain.iname] = node
        else
            push!(remaining_items, domain)
        end
    end

    # add remaining items as their dependencies are satisfied
    while !isempty(remaining_items)
        i = findfirst(remaining_items) do item
          good = true
          for dep in item.dependencies
              if !(dep in keys(nodes))
                  good = false
                  break
              end
          end
          good
        end

        if i == nothing
            error("items left but dependencies not satisfied: ", [s.iname for s in remaining_items])
        end

        item = remaining_items[i]
        deleteat!(remaining_items, i)

        node = AST(item, AST[], AST[])
        for dep in item.dependencies
            push!(node.parents, nodes[dep])
            push!(nodes[dep].children, node)
            nodes[item.iname] = node
        end
    end

    return ast

end


"""
topologically sort the AST into an order of instructions
modifies ast
"""
function topological_sort_order(ast::AST)::Vector{Union{Domain, Instruction}}
    order = Union{Domain, Instruction}[]
    sources = AST[head for head in ast.children]

    while !isempty(sources)
        # pop first to preserve original order if no dependencies
        n = popfirst!(sources)
        push!(order, n.self)
        while !isempty(n.children)
            child = popfirst!(n.children)
            filter!(e->e != n, child.parents)
            # child has no more dependencies
            if isempty(child.parents)
                push!(sources, child)
            end
        end
    end

    return order
end


"""
topologically sort the AST into a partial order of instructions
modifies ast
items in lists can run in any order
ex: [1, 2, [3, 4, 5], 6] describes that 3, 4, and 5 can be in any order
"""
function topological_sort_partial_order(ast::AST)::Vector{Union{Domain, Instruction, Vector}}
    order = Union{Domain, Instruction, Vector}[]
    temp_order = Union{Domain, Instruction, Vector}[]
    sources = AST[head for head in ast.children]
    second_sources = AST[]

    while !isempty(sources) || !isempty(second_sources)
        if isempty(sources)
            # add all in second_sources to sources
            append!(sources, second_sources)
            # clear second sources
            second_sources = AST[]
            # add temp_order to order
            if length(temp_order) == 1
                push!(order, temp_order[1])
            else
                # push as sub list
                push!(order, temp_order)
            end
            # clear temp_order
            temp_order = Union{Domain, Instruction, Vector}[]
        end
        # pop first to preserve original order if no dependencies
        n = popfirst!(sources)
        push!(temp_order, n.self)
        while !isempty(n.children)
            child = popfirst!(n.children)
            filter!(e->e != n, child.parents)
            # child has no more dependencies
            if isempty(child.parents)
                push!(second_sources, child)
            end
        end
    end

    if length(temp_order) == 1
        push!(order, temp_order[1])
    elseif length(temp_order) > 1
        # push as sub list
        push!(order, temp_order)
    end

    return order
end


"""
check if domains share instructions, and return those instructions
"""
function domains_shared_inst(domain1::Domain, domain2::Domain)::Vector{Union{Instruction, Domain}}
    shared = Union{Instruction, Domain}[]
    for inst1 in domain1.instructions
        for inst2 in domain2.instructions
            if inst1.iname == inst2.iname
                push!(shared, inst1)
            end
        end
    end
    return shared
end


"""
nest loops
"""
function nest_loops(kernel::LoopKernel, order::Vector{Union{Domain, Instruction}})::Vector{Union{Domain, Instruction}}

    new_order = Union{Domain, Instruction}[]

    nested_domains = Dict{Symbol, Bool}(s.iname=> false for s in kernel.domains)
    parents = Dict{Symbol, Domain}()

    for item in order
        if typeof(item) == Domain
            nested = false
            for domain in kernel.domains
                if domain != item
                    # get shared instructions
                    shared = domains_shared_inst(domain, item)
                    if length(shared) > 0
                        i = findfirst(e->e==item, order)
                        j = findfirst(e->e==domain, order)
                        # only modify if item is before domain (avoid repeats)
                        if i < j
                            # remove overlapping instructions from earlier domain
                            indices = findall(e->e in shared, item.instructions)
                            deleteat!(item.instructions, indices)
                            # add later domain into earlier domain's instructions
                            after = indices[1]
                            insert!(item.instructions, after, domain)
                            # append to new order if not already subset or in order
                            if domain.iname in keys(parents)
                                # domain is already nested, replace in nesting
                                parent = parents[domain.iname]
                                index = findall(e->e==domain, parent.instructions)[1]
                                # remove domain from old parent
                                deleteat!(parent.instructions, index)
                                # add item as new child
                                if !(item.iname in keys(parents))
                                    # nest where domain was
                                    insert!(parent.instructions, index, item)
                                    parents[item.iname] = parent
                                else
                                    # need to propogate up until no more parents to nest
                                    # find super parent of domain
                                    super_parent = parent
                                    while super_parent.iname in keys(parents)
                                        super_parent = parents[super_parent.iname]
                                    end
                                    item_parent = parents[item.iname]
                                    while item_parent.iname in keys(parents) && parents[item_parent.iname] != super_parent
                                        item_parent = parents[item_parent.iname]
                                    end

                                    if item_parent != parent
                                        # add highest up parent that isn't super parent to domain's immediate parent
                                        if !(item_parent in parent.instructions)
                                            insert!(parent.instructions, index, item_parent)
                                        end
                                        # remove item parent from old parent if one
                                        if item_parent.iname in keys(parents)
                                            p = parents[item_parent.iname]
                                            if p != parent
                                                deleteat!(p.instructions, findall(e->e==item_parent, p.instructions))
                                            end
                                        end
                                        parents[item_parent.iname] = parent

                                    end

                                    if item_parent in new_order
                                        deleteat!(new_order, findall(e->e==item_parent, new_order)[1])
                                    end
                                end
                            elseif !(item.iname in keys(parents)) && !(item in new_order)
                                # item is not nested in another or in order
                                push!(new_order, item)
                            end

                            # store parents
                            nested = true
                            nested_domains[domain.iname] = true
                            nested_domains[item.iname] = true
                            parents[domain.iname] = item
                        end
                    end
                end
            end
            if !nested && !nested_domains[item.iname]
                # shares no instructions with other loops
                push!(new_order, item)
            end
        else
            # instruction not in loop
            if length(item.dependencies) == 0
                push!(new_order, item)
            end
        end
    end

    return new_order

end


"""
SYMBOLS FUNCTIONS
"""

"""
get all symbols in a symbol
"""
get_all_symbols(s::Symbol)::Set{Symbol} = Set([s])


"""
get all symbols in a number
"""
get_all_symbols(n::Number)::Set{Symbol} = Set{Symbol}()

"""
get all symbols in a LineNumberNode
"""
get_all_symbols(n::LineNumberNode)::Set{Symbol} = Set{Symbol}()


"""
get all symbols in an expression
"""
function get_all_symbols(expr::Expr)::Set{Symbol}
    symbols = Set{Symbol}()
    if expr.head == :call
        # if function call, include all function args
        for arg in expr.args[2:end]
            union!(symbols, get_all_symbols(arg))
        end
    elseif expr.head == :macrocall
        for arg in expr.args[2:end]
            union!(symbols, get_all_symbols(arg))
        end
    else
        for arg in expr.args
            union!(symbols, get_all_symbols(arg))
        end
    end
    return symbols
end


"""
get all symbols in a kernel
"""
function get_all_symbols(kernel::LoopKernel)::Set{Symbol}
    symbols = Set{Symbol}()

    for instruction in kernel.instructions
        union!(symbols, get_all_symbols(instruction.body))
    end

    for domain in kernel.domains
        union!(symbols, get_all_symbols(domain.recurrence))
        union!(symbols, get_all_symbols(domain.lowerbound))
        union!(symbols, get_all_symbols(domain.upperbound))
        push!(symbols, domain.iname)
    end

    return symbols
end


"""
get all function arguments to kernel
"""
function get_kernel_args(kernel::LoopKernel)::Set{Symbol}
    all_symbols = get_all_symbols(kernel)
    defined_symbols = Set{Symbol}()

    # get all symbols defined in LHS of instructions
    for instruction in kernel.instructions
        lhs = instruction.body.args[1]
        if typeof(lhs) == Symbol && instruction.body.head == :(=)
            push!(defined_symbols, lhs)
        end
    end

    # get all loop domain symbols
    for domain in kernel.domains
        push!(defined_symbols, domain.iname)
    end

    # args are all symbols that are not defined symbols
    args = Set([s for s in all_symbols if !(s in defined_symbols)])

    append!(kernel.args, args)
    return args
end


"""
CONSTRUCTION FUNCTIONS
"""

"""
construct a domain
"""
function construct(domain::Domain)
    body = map(expr->construct(expr), domain.instructions)
    iname = domain.iname
    quote
        let $iname = $(domain.lowerbound)
            while $iname <= $(domain.upperbound)
                $(body...)
                $(domain.recurrence)
            end
        end
    end
end


"""
construct an instruction
"""
construct(instruction::Instruction) = instruction.body


"""
construct a list of instructions/domains
"""
function construct(list::Vector{Union{Instruction, Domain}})
    body = map(expr->construct(expr), list)
    quote
        $(body...)
    end
end


"""
COMPILATION
"""

"""
compile native julia code given a kernel
"""
function compile(kernel::LoopKernel)
    ast = construct_ast(kernel)
    order = topological_sort_order(ast)
    order = nest_loops(kernel, order)

    # construct from order
    body = construct(order)

    # kernel args
    args = get_kernel_args(kernel)

    # run_polyhedral_model(kernel)

    expr = quote
        function $(gensym(:JuLoop))(;$(args...))
            $(body)
        end
    end
    eval(expr)
end


# function printall(d)
#     if typeof(d) == Instruction
#         print(" ", d, " ")
#     else
#         print(" ", d, "-> ( ")
#         for i in d.instructions
#             printall(i)
#         end
#         print(" ) <-", d, " ")
#     end
# end

"""
compile native julia code given a kernel to an expression
"""
function compile_expr(kernel::LoopKernel)::Expr
    ast = construct_ast(kernel)
    order = topological_sort_order(ast)
    # ast = construct_ast(kernel)
    # order2 = topological_sort_partial_order(ast)
    order = nest_loops(kernel, order)

    # for i in order
    #     printall(i)
    # end

    # construct from order
    expr = construct(order)

    # kernel "args" for isl
    get_kernel_args(kernel)

    run_polyhedral_model(kernel)

    return expr
end
