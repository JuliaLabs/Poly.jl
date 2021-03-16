using MacroTools
using JuLoop

export compile_native_julia

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
            # if inexpr(kernel.instructions[i].body, kernel.instructions[j].iname)
            #     # j uses i directly, add dep
            #     push!(kernel.instructions[j].dependencies, kernel.instructions[i].iname)
            # end
            # check if there is a read/write dependency
            # read in i, write in j -> no dependency
            # read/read -> no dependency
            # write in i, read in j -> dependency
            # write in i, write in j -> no dependency
            # first check for function call (and if so, assume dependency)
            # TODO more intelligent dependencies with function calls
            if kernel.instructions[j].body.head == :call
                push!(kernel.instructions[j].dependencies, kernel.instructions[i].iname)
            else
                # need to check if LHS symbol that is being written in i is read in j
                lhs = kernel.instructions[j].body.args[1]
                while typeof(lhs) != Symbol
                    lhs = lhs.args[1]
                end
                rhs = kernel.instructions[j].body
                if rhs.head == :(=)
                    # only need rhs, since not +=, etc
                    rhs = rhs.args[2]
                end
                if inexpr(rhs, lhs)
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
constructs an AST using existing dependencies
"""
function construct_ast(kernel::LoopKernel)::AST
    add_impl_deps(kernel)
    add_loop_deps(kernel)

    ast = AST(nothing, AST[], AST[]) # head node
    remaining_items = Union{Instruction, Domain}[]
    nodes = Dict{Symbol, AST}()

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
            error("items left but dependencies not satisfied")
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
        n = popfirst!(sources)
        push!(order, n.self)
        while !isempty(n.children)
            child = popfirst!(n.children)
            filter!(e->e != n, child.parents)
            if isempty(child.parents)
                push!(sources, child)
            end
        end
    end

    return order
end


"""
check if domains share instructions
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

    for item in order
        if typeof(item) == Domain
            for domain in kernel.domains
                if domain != item
                    shared = domains_shared_inst(domain, item)
                    if length(shared) > 0
                        i = findfirst(e->e==item, order)
                        j = findfirst(e->e==domain, order)
                        if i < j
                            # remove overlapping instructions from earlier domain
                            indices = findall(e->e in shared, item.instructions)
                            deleteat!(item.instructions, indices)
                            # add later domain into earlier domain's instructions
                            after = indices[1]
                            insert!(item.instructions, after, domain)
                            # append to new order if not already subset
                            append = true
                            for base in new_order
                                if typeof(base) == Domain
                                    if item in base.instructions
                                        append = false
                                    end
                                end
                            end
                            if append
                                push!(new_order, item)
                            end
                        end
                    end
                end
            end
        else
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
get all symbols in symbol
"""
get_all_symbols(s::Symbol)::Set{Symbol} = Set([s])


"""
get all symbols in number
"""
get_all_symbols(n::Number)::Set{Symbol} = Set{Symbol}()


"""
get all symbols in expression
"""
function get_all_symbols(expr::Expr)::Set{Symbol}
    symbols = Set{Symbol}()
    if expr.head == :call
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
get all symbols in kernel
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

    for instruction in kernel.instructions
        lhs = instruction.body.args[1]
        if typeof(lhs) != Symbol
            continue
        end
        push!(defined_symbols, lhs)
    end

    for domain in kernel.domains
        push!(defined_symbols, domain.iname)
    end

    args = Set([s for s in all_symbols if !(s in defined_symbols)])

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
        $iname = $(domain.lowerbound)
        while $iname <= $(domain.upperbound)
            $(body...)
            $iname = $(domain.recurrence)
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
function compile_native_julia(kernel::LoopKernel)
    ast = construct_ast(kernel)
    @show ast
    order = topological_sort_order(ast)
    @show [e.iname for e in order]
    order = nest_loops(kernel, order)
    @show [e.iname for e in order]

    for d in kernel.domains
        @show d.iname, [s.iname for s in d.instructions]
    end

    # construct from order
    body = construct(order)

    # kernel args
    args = get_kernel_args(kernel)

    expr = quote
        function $(gensym(:JuLoop))(;$(args...))
            $(body)
        end
    end
    @show expr
    eval(expr)
end
