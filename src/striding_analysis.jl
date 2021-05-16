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
reorder the loops of a band node's partial schedule
"""
function reorder_band(band::Ptr{ISL.API.isl_schedule_node}, loop_ordering::Set{Vector{Symbol}}, context::Ptr{ISL.API.isl_ctx})
    # get schedule of band
    sched = ISL.API.isl_schedule_node_band_get_partial_schedule(band)
    # make a copy of the schedule (isl_copy does not return a copy but the actual reference)
    sched_str = ISL.API.isl_multi_union_pw_aff_to_str(sched)
    sched_str = Base.unsafe_convert(Ptr{Cchar}, sched_str)
    sched_str = Base.unsafe_string(sched_str)
    sched_copy = ISL.API.isl_multi_union_pw_aff_read_from_str(context, sched_str)
    # number of upa in sched
    n = ISL.API.isl_multi_union_pw_aff_dim(sched_copy, ISL.API.isl_dim_type(3))
    count = 0
    upas = []
    for i=0:n-1
        # read from copy
        upa = ISL.API.isl_multi_union_pw_aff_get_union_pw_aff(sched_copy, i)
        push!(upas, (false, upa))
    end

    # warm up
    # this is messy, but there is no way to set a partial schedule of a band
    # this sets the first element so that there are no more references known to the mupa
    # then, we can modify the mupa of the band directly
    upa_str = ISL.API.isl_union_pw_aff_to_str(upas[1][2])
    upa_str = Base.unsafe_convert(Ptr{Cchar}, upa_str)
    upa_str = Base.unsafe_string(upa_str)
    new_upa = ISL.API.isl_union_pw_aff_read_from_str(context, upa_str)
    ISL.API.isl_multi_union_pw_aff_set_union_pw_aff(sched, 0, new_upa)

    # add upas in order of loop orderings
    for ordering in loop_ordering
        for iname in ordering
            for (i, val) in enumerate(upas)
                used = val[1]
                upa = val[2]
                if used
                    continue
                end

                upa_str = ISL.API.isl_union_pw_aff_to_str(upa)
                upa_str = Base.unsafe_convert(Ptr{Cchar}, upa_str)
                upa_str = Base.unsafe_string(upa_str)
                # which iterator it is is between the -> and :
                # name[i0, i1, ...] -> [(i0)] : conditions...
                upa_iname_str = Base.split(upa_str, "->")[end]
                upa_iname_str = Base.split(upa_iname_str, ":")[1]
                # determine if iname in upa
                if occursin(sym_to_str(iname), upa_iname_str)
                    # make new upa (to avoid copying issues)
                    new_upa = ISL.API.isl_union_pw_aff_read_from_str(context, upa_str)
                    # set count location to upa in new schedule
                    ISL.API.isl_multi_union_pw_aff_set_union_pw_aff(sched, count, new_upa)
                    # increment count
                    count += 1
                    upas[i] = (true, upa)
                end
            end
        end
    end

    # add remaining upas
    for val in upas
        used = val[1]
        upa = val[2]
        if used
            continue
        end
        upa_str = ISL.API.isl_union_pw_aff_to_str(upa)
        upa_str = Base.unsafe_convert(Ptr{Cchar}, upa_str)
        upa_str = Base.unsafe_string(upa_str)
        # make new upa (to avoid copying issues)
        new_upa = ISL.API.isl_union_pw_aff_read_from_str(context, upa_str)
        # set count location to upa in new schedule
        ISL.API.isl_multi_union_pw_aff_set_union_pw_aff(sched, count, new_upa)
        # increment count
        count += 1
    end

end

"""
reorder the loops in a schedule

input:
    - kernel: kernel to reorder
    - schedule: ISL schedule of kernel
    - context: ISL context
    - loop_ordering: orderings of loops to use
    - loop_order_valid: if orderings are valid

returns reordered loop schedule if loop orderings are valid
"""
function reorder_schedule_loops(kernel::LoopKernel, schedule::Ptr{ISL.API.isl_schedule}, context::Ptr{ISL.API.isl_ctx}, loop_ordering::Set{Vector{Symbol}}, loop_order_valid::Bool)

    if !loop_order_valid
        return schedule
    end

    n = get_tile_dim(kernel)
    root = ISL.API.isl_schedule_get_root(schedule)
    node = root
    next_nodes = Set()
    # search for the outermost band node to reorder
    # (only band nodes can be reordered)
    while ISL.API.isl_schedule_node_n_children(node) > 0
        # go through all children of current node
        for i=0:ISL.API.isl_schedule_node_n_children(node)-1
            band = ISL.API.isl_schedule_node_get_child(node, i)
            if ISL.API.isl_schedule_node_get_type(band) == ISL.API.isl_schedule_node_type(0) # band node
                # in a band node, tile band with dimension n
                reorder_band(band, loop_ordering, context)
                # need to add child of band, since split band
                push!(next_nodes, ISL.API.isl_schedule_node_get_child(band, 0))
                break
            else
                # add to next nodes
                push!(next_nodes, band)
            end
        end

        if length(next_nodes) == 0
            break # done, no children left
        end
        node = pop!(next_nodes) # get next node
    end

end
