using JuLoop
using MacroTools

export get_best_nesting_order

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
        for arg in ex.args
            level += striding_level(arg, iter)
        end
    end
    return level
end

function striding_level(ex, iter::Symbol)::Int
    return 0
end


"""
Returns the overall striding level of domains for grouping of domains
"""
function get_loop_stride_counts(kernel::LoopKernel)::Dict{Symbol, Int}
    stride_counts = Dict{Symbol, Int}()
    for domain in kernel.domains
        stride_counts[domain.iname] = 0
        for instruction in kernel.instructions
            if typeof(instruction) == Instruction
                stride_counts[domain.iname] += striding_level(instruction.body, domain.iname)
            end
        end
    end
    return stride_counts
end

"""
Determines if domain can be inserted given list of already inserted domains
"""
function can_insert_domain(domain::Domain, order::Vector{Symbol}, kernel::LoopKernel)::Bool
    symbols = get_all_symbols(domain.lowerbound)
    union!(symbols, get_all_symbols(domain.upperbound))
    union!(symbols, get_all_symbols(domain.step))
    for sym in symbols
        if !(sym in order || sym in kernel.consts)
            return false
        end
    end
    return true
end

"""
helper for nesting order, returns nesting order for a grouping of domains starting at parent domain
"""
function get_best_nesting_order_on_domains(parent_domain::Domain, stride_order::Vector{Symbol}, kernel::LoopKernel)::Vector{Symbol}
    nesting_order = Symbol[]
    # first domain
    parent = parent_domain
    child = :none
    for instruction in parent.instructions
        if typeof(instruction) == Domain
            # found sub domain child
            child = instruction
            break
        end
    end
    # keep track of domains that can be reordered
    can_reorder = Set{Domain}()
    # add parent to right place
    if length(parent.instructions) == 1 && child != :none
        # can reorder later
        push!(can_reorder, parent)
    else
        # add in order now
        push!(nesting_order, parent.iname)
    end
    # keep reordering while sub domains exist
    while child != :none
        # can reorder if only one child
        if length(parent.instructions) == 1 # just one child domain
            push!(can_reorder, child)
        else # cannot reorder
            while length(can_reorder) > 0
                # find first in stride order and can reorder that can be added
                i = findfirst(stride_order) do item
                  good = false
                  for domain in can_reorder
                      if domain.iname == item && can_insert_domain(domain, nesting_order, kernel)
                          good = true
                          delete!(can_reorder, domain)
                          break
                      end
                  end
                  good
                end

                push!(nesting_order, stride_order[i])
            end
            # add children in right order
            for sub in parent.instructions
                if typeof(sub) == Domain
                    push!(nesting_order, sub.iname)
                end
            end
        end
        # get next parent and child
        parent = child
        child = :none
        for instruction in parent.instructions
            if typeof(instruction) == Domain
                child = instruction
                break
            end
        end
    end
    # done, add remaining parent and rest in reorder set
    while length(can_reorder) > 0
        # find first in stride order and can reorder that can be added
        i = findfirst(stride_order) do item
          good = false
          for domain in can_reorder
              if domain.iname == item && can_insert_domain(domain, nesting_order, kernel)
                  good = true
                  delete!(can_reorder, domain)
                  break
              end
          end
          good
        end

        push!(nesting_order, stride_order[i])
    end
    # add last domain
    if !(parent.iname in nesting_order)
        push!(nesting_order, parent.iname)
    end
    return nesting_order
end


"""
determines the best valid loop nesting order for groupings of domains (all in continuous nest (one sub domain per domain))
"valid" if loops have no other instructions between them
"best" if innermost loop has the least amount of striding
"""
function get_best_nesting_orders(kernel::LoopKernel)::Set{Vector{Symbol}}
    stride_counts = get_loop_stride_counts(kernel)
    stride_order = sort(collect(keys(stride_counts)), by=(x -> stride_counts[x]), rev=true) # sort least to most striding
    visited = Symbol[]
    nesting_orders = Set{Vector{Symbol}}()
    for domain in kernel.domains
        if !(domain.iname in visited)
            count = 0
            for child in domain.instructions
                if typeof(child) == Domain
                    count += 1
                end
            end
            if count == 1
                order = get_best_nesting_order_on_domains(domain, stride_order, kernel)
                append!(visited, order)
                push!(nesting_orders, order)
            end
        end
    end
    return nesting_orders
end

"""
reorders the loops in the original program to the best valid order
"""
function reorder_loops(kernel::LoopKernel)::Set{Vector{Symbol}}

end
