using JuLoop

"""
return the striding level of iter in ex
if i is the first dimension (i.e A[i, j], with i), returns 0
    etc for larger dimensions

ex: A[i, j] doesn't stride in i and strides in dimension 1 with j
"""
function striding_level(ex::Expr, iter::Symbol)::Int
    level = 0
    if ex.head == :ref
        for (i, sym) in enumerate(ex.args[2:end])
            if sym == iter
                level += i-1
            elseif typeof(sym) == Expr
                if inexpr(sym, iter)
                    level += i-1
                end
            end
        end
    else
        for arg in ex.args[2:end]
            level += striding_level(arg, iter)
        end
    end
end

function striding_level(ex, iter::Symbol)::Int
    return 0
end


"""
Returns the overall striding level of domains
"""
function get_loop_stride_counts(kernel::LoopKernel)::Dict{Symbol, Int}
    stride_counts = Dict{Symbol, Int}()
    for domain in kernel.domains
        stride_counts[domain.iname] = 0
        for instruction in domain.instructions
            stride_counts[domain.iname] += striding_level(instruction.body, domain.iname)
        end
    end
    return stride_counts
end
