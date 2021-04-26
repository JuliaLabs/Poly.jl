using JuLoop
export @poly_loop, @depends_on
using MacroTools


"""
helpers for evaluating an expression in mod
will error if symbols are not in global scope
"""
function eval_ex_in_mod(ex::Expr, mod::Module)::Number
    op = ex.args[1]
    new_args = []
    for arg in ex.args[2:end]
        push!(new_args, eval_ex_in_mod(arg, mod))
    end
    new_ex = Expr(ex.head, op, new_args...)
    return eval(new_ex)
end

function eval_ex_in_mod(sym::Symbol, mod::Module)
    return eval(quote $mod.$sym end)
end

function eval_ex_in_mod(num::Number, mod::Module)::Number
    return num
end

"""
helper for macro
"""
function poly_loop(ex::Expr, mod::Module)::LoopKernel
    # generate loop domain
    bounds = ex.args[1]
    symbol = bounds.args[1]
    space = bounds.args[2]
    if typeof(space) != Expr || !(length(space.args) == 3 || length(space.args) == 4)
        error("expected a loop in the form bound:bound or bound:step:bound")
    end
    if length(space.args) == 3
        # bound : bound
        lowerbound = 0
        upperbound = 0
        recurrence = 0
        try
            lowerbound = eval_ex_in_mod(space.args[2], mod)
            upperbound = eval_ex_in_mod(space.args[3], mod)
            recurrence = :($symbol += 1)
        catch e
            # symbols not in global scope
            lowerbound = space.args[2]
            upperbound = space.args[3]
            recurrence = :($symbol += 1)
        end
    else
        # bound : step : bound
        lowerbound = 0
        upperbound = 0
        recurrence = 0
        try
            lowerbound = eval_ex_in_mod(space.args[2], mod)
            upperbound = eval_ex_in_mod(space.args[4], mod)
            recurrence = :($symbol += eval_ex_in_mod(space.args[3], mod))
        catch e
            # symbols not in global scope
            lowerbound = space.args[2]
            upperbound = space.args[3]
            recurrence = :($symbol += 1)
        end
    end
    domains = [Domain(symbol, lowerbound, upperbound, recurrence, Set(), [])]

    # generate loop body instructions
    body = ex.args[2]
    instructions = []
    for line in body.args
        if typeof(line) == Expr && line.head == :for
            # nested loops
            subkern = poly_loop(line, mod)
            append!(instructions, subkern.instructions)
            append!(domains, subkern.domains)
        elseif typeof(line) == Expr && line.head == :macrocall && line.args[1] == Symbol("@depends_on")
            # dependency macro
            if length(line.args) != 4 || typeof(line.args[3]) != Expr || length(line.args[3].args) != 2
                error("use dependency with @depends_on elem=<Symbol> <Expr> only")
            end
            dependency = line.args[3].args[2]
            if typeof(dependency) != Symbol
                error("use dependency @depends_on with symbol only")
            end
            expr = line.args[4]
            symbol = expr
            while typeof(symbol) != Symbol
                symbol = symbol.args[1]
            end
            iname = gensym(string(symbol))
            inst = Instruction(iname, expr, Set([dependency]))
            push!(instructions, inst)
        elseif typeof(line) != LineNumberNode
            symbol = line
            while typeof(symbol) != Symbol
                symbol = symbol.args[1]
            end
            iname = gensym(string(symbol))
            inst = Instruction(iname, line, Set())
            push!(instructions, inst)
        end
    end

    kern = LoopKernel(instructions, domains, [], [], [])
    return kern
end


"""
Used to allow for code generation and restructuring of a for loop. May reorder some instructions (including loop orderings) and vectorize code. Only the outermost loop should use the macro. The other loops will be handled by the outermost macro

Example:
@poly_loop for i = 1:size(out, 1)
    for j = 1:size(out, 2)
        for k = 1:size(A, 2)
            out[i, j] += A[i, k] * B[k, j]
        end
        out[i, j] *= 2
    end
end
(constructs a matrix multiplication and double)

Notes:

@poly_loop currently requires a normal for loop (i.e i = lowerbound:upperbound or i = lowerbound:step:upperbound)

Loop bounds must NOT be function calls, UNLESS @poly_loop is in the global scope of a module. In the global scope, function calls will be evalutated in advance. If not in the global scope, expressions can be used like:
    n = size(out, 1)
    @poly_loop for i = 1:n ...

Dependencies are inferred by accesses. If a dependency is required that cannot be inferred from variable reads and writes, please use the @depends_on macro to specify a depencence

    Example:
    @poly_loop for i = 1:10
        @depends_on elem=i println("hello world") # otherwise, there is no easy way to tell if this can't be elevated from the loop, and transformation is agressive
    end

"""
macro poly_loop(ex0...)
    mod = __module__ # get module calling macro
    if length(ex0) != 1
        error("expected only a loop, got: ", length(ex0), "elements")
    end
    ex0 = ex0[1] # expression is first arg

    if ex0.head == :for
        kernel = poly_loop(ex0, mod)
    else
        error("expected a for loop, got: ", ex0.head)
    end

    # make kernel
    expr = MacroTools.prettify(compile_expr(kernel))
    esc(quote
        $(expr)
    end)
end


"""
placeholder macro for @poly_loop manual dependencies
Usage: @depends_on elem=<Symbol> <Expr> (not for loops)
see @poly_loop for more details
"""
macro depends_on(ex0...)

end
