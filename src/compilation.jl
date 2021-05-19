using MacroTools
using Poly

export compile


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
get all symbols in Any
"""
get_all_symbols(::Any)::Set{Symbol} = Set{Symbol}()

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
        union!(symbols, get_all_symbols(domain.step))
        union!(symbols, get_all_symbols(domain.lowerbound))
        union!(symbols, get_all_symbols(domain.upperbound))
        push!(symbols, domain.iname)
    end

    return symbols
end


"""
KERNEL ARG FUNCTIONS
"""

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
CONSTANT FINDING FUNCTIONS
"""

"""
determine all constants in kernel
"""
function set_kernel_consts(kernel::LoopKernel)
    const_symbols = Set{Symbol}()
    domain_inames = [domain.iname for domain in kernel.domains]

    # get all symbols in loop bounds (since must be constant)
    for domain in kernel.domains
        union!(const_symbols, get_kernel_consts(domain.lowerbound, domain_inames))
        union!(const_symbols, get_kernel_consts(domain.upperbound, domain_inames))
        union!(const_symbols, get_kernel_consts(domain.step, domain_inames))
    end

    append!(kernel.consts, const_symbols)
end

function get_kernel_consts(ex::Expr, domain_inames::Vector{Symbol})::Set{Symbol}
    const_symbols = Set{Symbol}()
    for arg in ex.args[2:end]
        union!(const_symbols, get_kernel_consts(arg, domain_inames))
    end
    return const_symbols
end

function get_kernel_consts(sym::Symbol, domain_inames::Vector{Symbol})::Set{Symbol}
    if !(sym in domain_inames)
        return Set([sym])
    end
    return Set()
end

function get_kernel_consts(any, domain_inames::Vector{Symbol})::Set{Symbol}
    return Set()
end


"""
COMPILATION
"""

"""
compile native julia code given a kernel
"""
function compile(kernel::LoopKernel; debug=false, verbose=0, tile=-1)
    # kernel args
    args = get_kernel_args(kernel)
    set_kernel_consts(kernel)

    body = run_polyhedral_model(kernel, debug=debug, verbose=verbose, tile=tile)

    expr = quote
        function $(gensym(:Poly))(;$(args...))
            $(body)
        end
    end
    eval(expr)
end

"""
compile native julia code given a kernel to an expression
"""
function compile_expr(kernel::LoopKernel; debug=false, verbose=0, tile=-1)::Expr
    # kernel "args" for isl
    get_kernel_args(kernel)
    set_kernel_consts(kernel)

    expr = run_polyhedral_model(kernel, debug=debug, verbose=verbose, tile=tile)

    return expr
end
