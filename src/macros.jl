using JuLoop
export @poly_loop

"""
helper for macro
"""
function poly_loop(ex::Expr)::LoopKernel
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
        elseif typeof(line) != LineNumberNode
            iname = gensym(:poly_instruction)
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
    expr = compile_expr(loopkern)
    esc(quote
        $(expr)
    end)
end
