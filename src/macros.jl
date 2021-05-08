using JuLoop
export @poly_loop, @depends_on, sub_sym_in_ex
using MacroTools


"""
helpers for substituting values for symbols in an expression
    will sub an interpolated symbol for the value as well as a non-interpolated symbol
"""
function sub_sym_in_ex(ex::Expr, sub_sym::Symbol, val::Number)
    if ex.head == :$ && ex.args[1] == sub_sym
        return val
    end
    new_args = []
    for arg in ex.args
        push!(new_args, sub_sym_in_ex(arg, sub_sym, val))
    end
    return Expr(ex.head, new_args...)
end

function sub_sym_in_ex(sym::Symbol, sub_sym::Symbol, val::Number)
    if sym == sub_sym
        return val
    end
    return sym
end

function sub_sym_in_ex(any, sub_sym::Symbol, val::Number)
    return any
end

"""
get instructions out of block with cond
"""
function get_insts_in_block(args, cond::Expr, parent_doms::Vector{Domain}, dom::Domain)::Vector{Instruction}
    instructions = []
    for arg in args
        if typeof(arg) != LineNumberNode
            symbol = arg
            while typeof(symbol) != Symbol
                symbol = symbol.args[1]
            end
            iname = gensym(string(symbol))
            inst = Instruction(iname, arg, Set(), cond)
            for d in parent_doms
                push!(inst.dependencies, d.iname)
            end
            push!(inst.dependencies, dom.iname)
            push!(dom.instructions, inst)
            push!(instructions, inst)
        end
    end
    return instructions
end

"""
helper for macro
"""
function poly_loop(ex::Expr; parent_doms=Domain[])::LoopKernel
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
            subkern = poly_loop(line, parent_doms=new_parent_doms)
            append!(instructions, subkern.instructions)
            append!(domains, subkern.domains)
            push!(dom.instructions, subkern.domains[1])
            for d in parent_doms
                push!(subkern.domains[1].dependencies, d.iname)
            end
            push!(subkern.domains[1].dependencies, dom.iname)
        elseif typeof(line) == Expr && line.head == :if
            # if block
            ifcond = line.args[1]
            # instructions in if block
            append!(instructions, get_insts_in_block(line.args[2].args, ifcond, parent_doms, dom))
            # elseif or else or both
            if length(line.args) == 3
                rest = line.args[3]
                if rest.head == :elseif
                    # elseif block
                    elifcond = rest.args[1].args[2] # since line number node
                    elifcond = :($elifcond && !($ifcond))
                    # instructions in elseif block
                    append!(instructions, get_insts_in_block(rest.args[2].args, elifcond, parent_doms, dom))
                    # else block also
                    if length(rest.args) == 3
                        # else block
                        elsecond = :(!($ifcond) && !($elifcond))
                        # instructions in else block
                        append!(instructions, get_insts_in_block(rest.args[3].args, elsecond, parent_doms, dom))
                    end
                else
                    # else block
                    elsecond = :(!($ifcond))
                    # instructions in else block
                    append!(instructions, get_insts_in_block(rest.args, elsecond, parent_doms, dom))
                end
            end
        elseif typeof(line) != LineNumberNode
            symbol = line
            while typeof(symbol) != Symbol
                symbol = symbol.args[1]
            end
            iname = gensym(string(symbol))
            inst = Instruction(iname, line, Set(), :())
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
helper to get the value in a value
"""
function __get(::Type{Val{N}}) where N
    return N
end


"""
Used to allow for code generation and restructuring of a for loop. May reorder some instructions (including loop orderings) and vectorize code. Only the outermost loop should use the macro. The other loops will be handled by the outermost macro

Example:
@poly_loop for i = 1:n
    for j = 1:m
        for k = 1:r
            out[i, j] += A[i, k] * B[k, j]
        end
        out[i, j] *= 2
    end
end
(constructs a matrix multiplication and double)

Notes:

@poly_loop requires a normal for loop (i.e i = lowerbound:upperbound or i = lowerbound:step:upperbound)

Loop bounds:
    Loop bounds must NOT be function calls (such as size), since ISL cannot evaluate them. One easy workaround is something like:
        n = size(out, 1)
        @poly_loop for i = 1:n ...

Loop step:

Debugging:
    Debugging and verbosity options (see run_polyhedral_model for details) can be passed as:
        @poly_loop debug=true verbose=2 for ...

Functions:
    Since dependencies are inferred from accesses, function calls are currently NOT supported.
        If a function modifies any inputs, then there is no way to know which inputs are modified or how those inputs are modified. If a function returns an output such as a matrix, the write to each individual index of the matrix (for example A[i, j]) cannot be inferrred by ISL. In addition, @inbounds is not neccessary as the aggresive transformation will add @inbounds automatically.
        The only exceptions are the following (which can be used in loop bounds and in instructions):
            min():
                @poly_loop for i=1:min(n, m) ...
            max():
                @poly_loop for i=1:max(n, m) ...
            floor() -> DO NOT cast to Int() (will be added later):
                @poly_loop for i=1:floor(n/2) ...

Loop names:
    All loop iterators must have unique names, even if in different scopes. This is so an original schedule can be extracted from the code.

Interpolation:
    ISL does not work with non-numerical striding. So, this macro supports interpolating the values into the function so that ISL has access to the numerical values (and non-numerical strides MUST be interpolated). This delays the analysis to runtime, so it is advised that interpolated values are typically constant so that compilation does not occur over and over
    Ex:
        @poly_loop for i=1:\$c:n
            arr[i] = 1
        end
    In addition, any constants can be interpolated, and sometimes this may improve the performance of ISL. Numerical bounds can help with tiling and inference on spaces. However, it is important that interpolated values are actually constants, or very rarely change, so that this delayed compilation does not occur regularly during runtime.

Indexing:
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
            @poly_loop for t=1:\$TILE_SIZE:NUM_TILES
                A[t] = 1
            end

If blocks:
    If blocks must contain conditions on the loop iterators and constants only (not on any elements of an array, for example)

"""
macro poly_loop(ex0...)
    if length(ex0) > 3
        error("expected at most 2 options and a loop, got: ", length(ex0), "elements")
    end
    debug = false
    verbose = 0
    if length(ex0) == 2
        arg = ex0[1]
        if arg.args[1] == :debug
            debug = arg.args[2]
        elseif arg.args[1] == :verbose
            verbose = arg.args[2]
        end
    elseif length(ex0) == 3
        arg1 = ex0[1]
        if arg1.args[1] == :debug
            debug = arg1.args[2]
        elseif arg1.args[1] == :verbose
            verbose = arg1.args[2]
        end
        arg2 = ex0[2]
        if arg2.args[1] == :debug
            debug = arg2.args[2]
        elseif arg2.args[1] == :verbose
            verbose = arg2.args[2]
        end
    end
    if typeof(debug) != Bool
        error("expected boolean for debug, got: ", debug)
    end
    if typeof(verbose) != Int
        error("expected integer for verbose, got: ", verbose)
    end

    ex0 = ex0[end] # expression is last arg

    if ex0.head != :for
        error("expected a for loop, got: ", ex0.head)
    end

    liftedex = copy(ex0)
    letargs = Base._lift_one_interp!(liftedex)
    if letargs != []
        # have interpolation, need to evaluate at runtime to get values of symbols
        sub_symbols = []
        # gather all the symbols that need to be replaced
        for arg in letargs
            push!(sub_symbols, arg.args[2].args[1])
        end
        # compile the Val() expressions that will be passed to the generated function
        exprs = []
        for sym in sub_symbols
            push!(exprs, Expr(:call, :Val, Expr(:escape, sym)))
        end
        quote
            # generated function will have the values of the symbols for replacing at runtime
            @generated function delay_poly_loop(vals...)
                expr = $(QuoteNode(ex0))
                syms = $(QuoteNode(sub_symbols))
                d = $(QuoteNode(debug))
                v = $(QuoteNode(verbose))
                for (i, val) in enumerate(vals)
                    N = __get(val)
                    expr = sub_sym_in_ex(expr, syms[i], N)
                end
                kernel = poly_loop(expr)
                expr = MacroTools.prettify(compile_expr(kernel, debug=d, verbose=v))
                # use the compiled expression
                expr
            end
            delay_poly_loop($(exprs...))
        end
    else
        # make kernel
        kernel = poly_loop(ex0)
        expr = MacroTools.prettify(compile_expr(kernel, debug=debug, verbose=verbose))
        esc(quote
            $(expr)
        end)
    end
end
