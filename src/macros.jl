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
function poly_loop(ex::Expr, mod::Module, parent_doms=[])::LoopKernel
    # generate loop domain
    bounds = ex.args[1]
    iname = bounds.args[1]
    space = bounds.args[2]
    if typeof(space) != Expr || !(length(space.args) == 3 || length(space.args) == 4)
        error("expected a loop in the form bound:bound or bound:step:bound")
    end
    if length(space.args) == 3
        # bound : bound
        lowerbound = space.args[2]
        upperbound = space.args[3]
        step = 1
    else
        # bound : step : bound
        lowerbound = space.args[2]
        upperbound = space.args[4]
        step = space.args[3]
    end
    dom = Domain(iname, lowerbound, upperbound, step, Set(), [])
    domains = [dom]
    new_parent_doms = copy(parent_doms)
    append!(new_parent_doms, [dom])

    # generate loop body instructions
    body = ex.args[2]
    instructions = []
    for line in body.args
        if typeof(line) == Expr && line.head == :for
            # nested loops
            subkern = poly_loop(line, mod, new_parent_doms)
            append!(instructions, subkern.instructions)
            append!(domains, subkern.domains)
            push!(dom.instructions, subkern.domains[1])
            for d in parent_doms
                push!(subkern.domains[1].dependencies, d.iname)
            end
            push!(subkern.domains[1].dependencies, dom.iname)
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
            for d in parent_doms
                push!(inst.dependencies, d.iname)
            end
            push!(inst.dependencies, dom.iname)
            push!(dom.instructions, inst)
        elseif typeof(line) != LineNumberNode
            symbol = line
            while typeof(symbol) != Symbol
                symbol = symbol.args[1]
            end
            iname = gensym(string(symbol))
            inst = Instruction(iname, line, Set())
            push!(instructions, inst)
            for d in parent_doms
                push!(inst.dependencies, d.iname)
            end
            push!(inst.dependencies, dom.iname)
            push!(dom.instructions, inst)
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

Loop bounds must NOT be function calls (such as size), since ISL cannot evaluate them. One easy workaround is something like:
    n = size(out, 1)
    @poly_loop for i = 1:n ...

Dependencies are inferred by accesses. If a dependency is required that cannot be inferred from variable reads and writes, please use the @depends_on macro to specify a depencence

    Example:
    @poly_loop for i = 1:10
        @depends_on elem=i println("hello world") # otherwise, there is no easy way to tell if this can't be elevated from the loop, and transformation is agressive
    end

Since dependencies are inferred from accesses, function calls are currently NOT supported. If a function modifies any inputs, then there is no way to know which inputs are modified or how those inputs are modified. If a function returns an output such as a matrix, the write to each individual index of the matrix (for example A[i, j]) cannot be inferrred by ISL. In addition, @inbounds is not neccessary as the aggresive transformation will add @inbounds automatically.

All loop iterators must have unique names, even if in different scopes. This is so an original schedule can be extracted from the code.

All array indices must be in terms of loop bounds or constants, not of variables. For example,
        c = i + j
        A[c]
    Must become
        A[i + j]

Multilpicative indexing (i.e. A[i*t, j]) is not supported. If needed, stride over the iterator.
    So, this:
        @poly_loop for t=1:NUM_TILES
            A[TILE_SIZE*t] = 1
        end
    Would need to become this:
        @poly_loop for t=1:TILE_SIZE:NUM_TILES
            A[t] = 1
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
