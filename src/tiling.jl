using JuLoop
using Hwloc
using ISL
import ISL.API

"""
get the best tiling of arrays of n dimensions, on band_node, using l1 cache sizes only
assumes float64 arrays (worst-case sizing)
trys to find a tile size tile_1 x tile_2 x .... x tile_t where tile_1 x ... x tile_n <= l1cachesize (the array being tiled will fit in the l1 cache) and tile_i <= l1linesize for all i (all tiled dimensions will fit in a cache line)
finds square tiles

returns: isl_multi_val representation of tiling sizes
"""
function tiling_sizes(n, band::Ptr{ISL.API.isl_schedule_node}, context::Ptr{ISL.API.isl_ctx})::Ptr{ISL.API.isl_multi_val}
    size = Hwloc.l1cache_sizes()[1]/sizeof(Float64)
    linesize = Hwloc.l1cache_linesizes()[1]/sizeof(Float64)
    tiles = Int[]
    tile_size = Int(floor(min(linesize, size^(1/n))))
    band_space = ISL.API.isl_schedule_node_band_get_space(band)
    dim = ISL.API.isl_space_dim(band_space, ISL.API.isl_dim_type(3))
    tiles = ISL.API.isl_val_list_alloc(context, dim)
    for i=1:dim
        ISL.API.isl_val_list_add(tiles, ISL.API.isl_val_int_from_si(context, tile_size))
    end
    tiles = ISL.API.isl_multi_val_from_val_list(band_space, tiles)
    return tiles
end

"""
shift the band to start at 0 so that full tiling occurs
"""
function shift_band_zero(band::Ptr{ISL.API.isl_schedule_node}) ::Tuple{Ptr{ISL.API.isl_multi_union_pw_aff}, Ptr{ISL.API.isl_multi_union_pw_aff}}
    # band domain
    domain = ISL.API.isl_schedule_node_get_domain(band)
    println("BAND DOM:")
    ISL.API.isl_union_set_dump(domain)
    # mupa partial schedule
    partial_schedule = ISL.API.isl_schedule_node_band_get_partial_schedule(band)
    println("CURRENT PARTIAL:")
    mupa = ISL.API.isl_multi_union_pw_aff_intersect_domain(ISL.API.isl_multi_union_pw_aff_copy(partial_schedule), ISL.API.isl_union_set_copy(domain))
    ISL.API.isl_multi_union_pw_aff_dump(mupa)
    # get shift mupa
    multi_val = ISL.API.isl_multi_union_pw_aff_min_multi_val(mupa)
    println("MV:")
    ISL.API.isl_multi_val_dump(multi_val)
    shift = ISL.API.isl_multi_union_pw_aff_multi_val_on_domain(domain, multi_val)
    shift_neg = ISL.API.isl_multi_union_pw_aff_neg(ISL.API.isl_multi_union_pw_aff_copy(shift))
    # shift partial schedule
    partial_schedule = ISL.API.isl_multi_union_pw_aff_add(partial_schedule, shift_neg)
    println("NEW PARTIAL:")
    ISL.API.isl_multi_union_pw_aff_dump(partial_schedule)
    # # shift back
    # partial_schedule_back = ISL.API.isl_multi_union_pw_aff_add(partial_schedule, ISL.API.isl_multi_union_pw_aff_copy(shift))
    # println("NEW NEW PARTIAL:")
    # ISL.API.isl_multi_union_pw_aff_dump(partial_schedule_back)
    # put new schedule in band
    # band = ISL.API.isl_schedule_node_band_set_partial_schedule(partial_schedule)
    # keep shift (mupa) to shift back up later
    return partial_schedule, shift
end

"""
tile a partial schedule
"""
function tile_partial_schedule(sched::Ptr{ISL.API.isl_multi_union_pw_aff}, sizes::Ptr{ISL.API.isl_multi_val})::Ptr{ISL.API.isl_multi_union_pw_aff}
    ctx = ISL.API.isl_multi_val_get_ctx(sizes)
    scale = ISL.API.isl_options_get_tile_scale_tile_loops(ctx)
    n = ISL.API.isl_multi_union_pw_aff_size(sched)
    for i=0:n-1
        upa = ISL.API.isl_multi_union_pw_aff_get_union_pw_aff(sched, i)
        v = ISL.API.isl_multi_val_get_val(sizes, i)
        upa = ISL.API.isl_union_pw_aff_scale_down_val(upa, ISL.API.isl_val_copy(v))
		upa = ISL.API.isl_union_pw_aff_floor(upa)
		if scale == ISL.API.isl_bool_true
			upa = ISL.API.isl_union_pw_aff_scale_val(upa, ISL.API.isl_val_copy(v))
        end
        ISL.API.isl_val_free(v)
        sched = ISL.API.isl_multi_union_pw_aff_set_union_pw_aff(sched, i, upa)
    end
    return sched
end

"""
tile a band node with tile dimension n (max dimension of tile)
"""
function tile_band(n, band::Ptr{ISL.API.isl_schedule_node}, context::Ptr{ISL.API.isl_ctx})::Ptr{ISL.API.isl_schedule}
    # shift to zero
    partial_schedule, shift = shift_band_zero(band)
    # tile
    multi_val = tiling_sizes(n, band, context)
    partial_schedule = tile_partial_schedule(partial_schedule, multi_val)
    # band = ISL.API.isl_schedule_node_band_tile(band, multival)
    # shift back
    # partial_schedule = ISL.API.isl_schedule_node_band_get_partial_schedule(band)
    partial_schedule = ISL.API.isl_multi_union_pw_aff_add(partial_schedule, ISL.API.isl_multi_union_pw_aff_copy(shift))
    # insert tiled schedule below band
    band = ISL.API.isl_schedule_node_insert_partial_schedule(band, partial_schedule)
    # get new schedule
    schedule = ISL.API.isl_schedule_node_get_schedule(band)
    return schedule
end

"""
get the max tile dimension
"""
function get_tile_dim(kernel::LoopKernel)::Int
    max_dim = 0
    for instruction in kernel.instructions
        dim = get_tile_dim(instruction.body)
        max_dim = max(max_dim, dim)
    end
    return max_dim
end

function get_tile_dim(ex::Expr)::Int
    max_dim = 0
    if ex.head == :ref
        # care about dimension of arrays, so references
        return length(ex.args[2:end]) # number of dimensions in ref
    else
        for arg in ex.args
            max_dim = max(get_tile_dim(arg), max_dim)
        end
    end
    return max_dim
end

function get_tile_dim(any)::Int
    return 0
end

"""
tile a schedule
"""
function tile_schedule(kernel::LoopKernel, schedule::Ptr{ISL.API.isl_schedule}, context::Ptr{ISL.API.isl_ctx}, tile_dim::Int)::Ptr{ISL.API.isl_schedule}
    tile_schedule = schedule
    n = get_tile_dim(kernel)
    root = ISL.API.isl_schedule_get_root(schedule)
    node = root
    next_nodes = Set()
    count = 0
    tiled_level = 0
    while ISL.API.isl_schedule_node_n_children(node) > 0
        count += 1
        # go through all children of current node
        tiled = false
        for i=0:ISL.API.isl_schedule_node_n_children(node)-1
            band = ISL.API.isl_schedule_node_get_child(node, i)
            if ISL.API.isl_schedule_node_get_type(band) == ISL.API.isl_schedule_node_type(0) # band node
                # in a band node, tile band with dimension n
                tile_schedule = tile_band(n, band, context)
                # need to add child of band, since split band
                push!(next_nodes, ISL.API.isl_schedule_node_get_child(band, 0))
                tiled = true
                break
            else
                # add to next nodes
                push!(next_nodes, band)
            end
        end

        if tiled
            # tiled_level += 1
            break
        end

        if tiled_level == tile_dim
            break # reached max tiling depth, done
        end

        if length(next_nodes) == 0
            break # done, no children left
        end
        node = pop!(next_nodes) # get next node
        if count > 10
            break
        end
    end

    return tile_schedule
end
