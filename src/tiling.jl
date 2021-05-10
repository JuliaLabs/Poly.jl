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
tile a band node with tile dimension n (max dimension of tile)
"""
function tile_band(n, band::Ptr{ISL.API.isl_schedule_node}, context::Ptr{ISL.API.isl_ctx})::Ptr{ISL.API.isl_schedule}
    multival = tiling_sizes(n, band, context)
    band = ISL.API.isl_schedule_node_band_tile(band, multival)
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
            else
                # add to next nodes
                push!(next_nodes, band)
            end
        end

        if tiled
            tiled_level += 1
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
