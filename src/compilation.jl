using MacroTools
import datatypes.jl
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
                append!(inst.dependencies, domain.iname)
            end
        end
    end
    for i = 1:len(kernel.domains)
        domain1 = kernel.domains[i]
        for j = i+1:len(kernel.domains)
            domain2 = kernel.domains[j]
            if inexpr(domain1.recurrence, domain2.iname)
                append!(domain2.dependencies, domain1.iname)
            elseif isexpr(domain1.lowerbound) && inexpr(domain1.lowerbound, domain2.iname)
                append!(domain2.dependencies, domain1.iname)
            elseif isexpr(domain1.upperbound) && inexpr(domain1.upperbound, domain2.iname)
                append!(domain2.dependencies, domain1.iname)
            end
        end
    end
    # TODO: add dependencies between loops due to instructions?
end


"""
uses the implicit order of kernel instructions to add dependencies between instructions
"""
function add_impl_deps(kernel::LoopKernel)
    for i = 1:len(kernel.instructions)
        # look only in instructions following i
        for j = i+1:len(kernel.instructions)
            # if inexpr(kernel.instructions[i].body, kernel.instructions[j].id)
            #     # j uses i directly, add dep
            #     append!(kernel.instructions[j].dependencies, kernel.instructions[i].id)
            # end
            # check if there is a read/write dependency
            # read in i, write in j -> no dependency
            # read/read -> no dependency
            # write in i, read in j -> dependency
            # write in i, write in j -> no dependency
            # first check for function call (and if so, assume dependency)
            if kernel.instructions[j].head == :call
                append!(kernel.instructions[j].dependencies, kernel.instructions[i].id)
            else
                # need to check if LHS symbol that is being written in i is read in j
                lhs = kernel.instructions[j].body.args[1]
                while typeof(lhs) != :Symbol
                    lhs = lhs.args[1]
                end
                rhs = kernel.instructions[j].body
                if rhs.head == :(=)
                    # only need rhs, since not +=, etc
                    rhs = rhs.args[2]
                end
                if inexpr(rhs, lhs)
                    append!(kernel.instructions[j].dependencies, kernel.instructions[i].id)
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
    add_impl_deps()
    add_loop_deps()

    ast = AST{:none, Vector{AST}[], Vector{AST}[]} # head node
    remaining_items = Vector{Union{Instruction, Domain}}[]
    nodes = Dict{Symbol, AST}()

    for instruction in kernel.instructions
        if isempty(instruction.dependencies)
            node = AST{instruction, Vector{AST}[], Vector{AST}[]}
            append!(ast.children, node)
            append!(node.parents, ast)
            nodes[instruction] = node
        else
            append!(remaining_items, instruction)
        end

    end

    for domain in kernel.domains
        if isempty(domain.dependencies)
            node = AST{domain, Vector{AST}[], Vector{AST}[]}
            append!(ast.children, node)
            append!(node.parents, ast)
            nodes[domain] = node
        else
            append(remaining_items, domain)
        end
    end


    while !isempty(remaining_items)
        next = :none
        for item in remaining_items
            good = true
            for dep in item.dependencies
                if !(dep in keys(nodes))
                    good = false
                end
            end
            if good
                next = pop!(remaining_items, item)
                break
            end
        end

        if next != :none
            node = AST{next, Vector{AST}[], Vector{AST}[]}
            for dep in next.dependencies
                append!(node.parents, nodes[dep])
                append!(nodes[dep].children, node)
            end
        else
            break
        end
    end

    return ast

end


"""
topologically sort the AST into an order of instructions
will modify ast
"""
function topological_sort_order(ast::AST)::Vector{Union{Domain, Instruction}}
    order = Vector{Union{Domain, Instruction}[]
    sources = Vector{Union{Domain, Instruction}[head for head in ast.children]

    while !isempty(sources)
        n = pop!(sources, sources[1])
        append!(order, n.self)
        while !isempty(n.children)
            child = pop!(n.children, n.children[1])
            remove!(child.parents, n)
            if isempty(child.parents)
                if typeof(n.self) == Domain
                    append!(n.self.instructions, child.self)
                else
                    append!(sources, child)
                end
            end
        end
    end

    return order
end


"""
SYMBOLS FUNCTIONS
"""

"""
get all symbols in symbol
"""
get_all_symbols(s::Symbol)::Vector{Symbol} = [s]


"""
get all symbols in number
"""
get_all_symbols(n::Number)::Vector{Symbol} = []


"""
get all symbols in expression
"""
function get_all_symbols(expr::Expr)::Vector{Symbol}
    symbols = Vector{Symbol}[]
    for arg in expr.args
        symbols = vcat(symbols, get_all_symbols(arg))
    end
    return symbols
end


"""
get all symbols in kernel
"""
function get_all_symbols(kernel::LoopKernel)::Vector{Symbol}
    symbols = Vector{Symbol}[]

    for instruction in kernel.instructions
        symbols = vcat(symbols, get_all_symbols(instruction.body))
    end

    for domain in kernel.domains
        symbols = vcat(symbols, get_all_symbols(domain.recurrence))
        symbols = vcat(symbols, get_all_symbols(domain.lowerbound))
        symbols = vcat(symbols, get_all_symbols(domain.upperbound))
        append!(symbols, domain.iname)
    end

    return symbols
end


"""
get all function arguments to kernel
"""
function get_kernel_args(kernel::LoopKernel)::Vector{Symbol}
    all_symbols = get_all_symbols(kernel)
    defined_symbols = Vector{Symbol}[]

    for instruction in kernel.instructions
        lhs = instruction.body.args[1]
        append!(defined_symbols, lhs)
    end

    for domain in kernel.domains
        append!(defined_symbols, domain.iname)
    end

    args = [s for symbol in all_symbols if !(symbol in defined_symbols)]

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
    order = topological_sort_order(ast)

    # construct from order
    body = construct(order)

    # kernel args
    args = get_kernel_args(kernel)

    expr = quote
        function $(gensym())($(args...))
            $(body)
        end
    end
    @show expr
    eval(expr)
end
