using JuLoop
using Test
using StaticArrays
using BenchmarkTools

function poly_mul(A, B, out)
    n = size(out, 1)
    m = size(out, 2)
    r = size(A, 2)
    @poly_loop for i = 1:n
        for j = 1:m
            for k = 1:r
                out[i, j] += A[i, k] * B[k, j]
            end
        end
    end
end

function mul(A, B, out)
    n = size(out, 1)
    m = size(out, 2)
    r = size(A, 2)
    @inbounds for i = 1:n
        for j = 1:m
            for k = 1:r
                out[i, j] += A[i, k] * B[k, j]
            end
        end
    end
end

function mul_right_order(A, B, out)
    n = size(out, 1)
    m = size(out, 2)
    r = size(A, 2)
    @inbounds for j = 1:m
        for k = 1:r
            for i = 1:n
                out[i, j] += A[i, k] * B[k, j]
            end
        end
    end
end

function tiled_mul(A, B, C)
    N = size(C)[1]
    TILE_DIM = 32
    NUM_TILES = Int(N/TILE_DIM)

    @inbounds for gj = 1:NUM_TILES
        for gi = 1:NUM_TILES
            # loop over tiles needed for this calculation
            tile1 = @MArray zeros(TILE_DIM, TILE_DIM)
            tile2 = @MArray zeros(TILE_DIM, TILE_DIM)
            for t = 0:NUM_TILES-1
                # load tiles needed for calculation
                for j = 1:TILE_DIM
                    for i = 1:TILE_DIM
                        # global tile
                        I = (gi-1) * TILE_DIM + i
                        J = (gj-1) * TILE_DIM + j
                        # get tile1 and tile2 values
                        tile1[i, j] = A[I, t*TILE_DIM + j]
                        tile2[i, j] = B[t*TILE_DIM + i, J]
                    end
                end
                # synchronize
                # loop over tiles to calculate for I, J spot
                for jj in 1:TILE_DIM
                    JJ = (gj-1) * TILE_DIM + jj
                    # loop over row/col in tiles
                    for k = 1:TILE_DIM
                        for ii = 1:TILE_DIM
                            # global tile
                            II = (gi-1) * TILE_DIM + ii
                            # add tile1 * tile2
                            C[II, JJ] += tile1[ii, k] * tile2[k, jj]
                        end
                    end
                end
            end
        end
    end
end

@testset "JuLoop.jl" begin

    # test matrix mul wrong loop order (check automatic reordering)
    DIM = 256
    A = rand(DIM, DIM)
    B = rand(DIM, DIM)
    out = zeros(DIM, DIM)

    t_poly = @benchmark poly_mul($A, $B, out) setup=(out = zeros(256, 256))

    t_orig = @benchmark mul($A, $B, out) setup=(out = zeros(256, 256))

    t_right_order = @benchmark mul_right_order($A, $B, out) setup=(out = zeros(256, 256))

    r = ratio(minimum(t_orig), minimum(t_poly))
    j = judge(minimum(t_poly), minimum(t_right_order))

    @test r.time > 2 # want at least 2x speedup over bad order
    @test j.time != :regression # want comparable or better than right order


    # test performance nears tiled matrix multiplication (check automatic tiling)
    DIM = 512
    A = rand(DIM, DIM)
    B = rand(DIM, DIM)
    out = zeros(DIM, DIM)

    t_poly = @benchmark poly_mul($A, $B, out) setup=(out = zeros(512, 512))

    t_tile = @benchmark tiled_mul($A, $B, out) setup=(out = zeros(512, 512))

    j = judge(minimum(t_poly), minimum(t_tile))

    @test j.time != :regression # want improvement or invariant (not regression)

end
