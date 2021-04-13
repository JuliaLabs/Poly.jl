using JuLoop
export @poly_loop, @depends_on
using MacroTools

"""
helper for macro
"""
function poly_loop(ex::Expr)::LoopKernel
    ex = MacroTools.rmlines(ex)
    # generate loop domain
    bounds = ex.args[1]
    symbol = bounds.args[1]
    space = bounds.args[2]
    if typeof(space) != Expr || !(length(space.args) == 3 || length(space.args) == 4)
        error("expected a loop in the form bound:bound or bound:step:bound")
    end
    if length(space.args) == 3
        # bound : bound
        lowerbound = space.args[2]
        upperbound = space.args[3]
        recurrence = :($symbol += 1)
    else
        # bound : step : bound
        lowerbound = space.args[2]
        upperbound = space.args[4]
        recurrence = :($symbol += space.args[3])
    end
    domains = [Domain(symbol, lowerbound, upperbound, recurrence, Set(), [])]

    # generate loop body instructions
    body = ex.args[2]
    instructions = []
    for line in body.args
        if typeof(line) == Expr && line.head == :for
            # nested loops
            subkern = poly_loop(line)
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

    kern = LoopKernel(instructions, domains, [], [])
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

Dependencies are inferred by accesses. If a dependency is required that cannot be inferred from variable reads and writes, please use the @depends_on macro to specify a depencence

    Example:
    @poly_loop for i = 1:10
        @depends_on elem=i println("hello world")
    end

"""
macro poly_loop(ex0...)
    if length(ex0) != 1
        error("expected only a loop, got", length(ex0), "elements")
    end
    ex0 = ex0[1] # expression is first arg

    if ex0.head == :for
        loopkern = poly_loop(ex0)
    else
        error("expected a for loop, got: ", ex0.head)
    end

    # make kernel
    expr = MacroTools.prettify(compile_expr(loopkern))
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
