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

function poly_array_double(A, out)
    n = size(out, 1)

    @poly_loop for i = 1:n
        out[i] = A[i]*2
    end
end

function array_double(A, out)
    n = size(out, 1)

    for i = 1:n
        out[i] = A[i]*2
    end
end

function poly_matmul_double(A, B, out)
    n = size(out, 1)
    m = size(out, 2)
    r = size(A, 2)

    @poly_loop for i = 1:n
        for j = 1:m
            for k = 1:r
                out[i, j] += A[i, k] * B[k, j]
            end
            out[i, j] *= 2
        end
    end
end

function matmul_double(A, B, out)
    n = size(out, 1)
    m = size(out, 2)
    r = size(A, 2)

    for i = 1:n
        for j = 1:m
            for k = 1:r
                out[i, j] += A[i, k] * B[k, j]
            end
            out[i, j] *= 2
        end
    end
end

function poly_dd_stride(arr)
    n = size(arr, 1)
    @poly_loop for i = 1:2:n
        for j = 1:2:n
            arr[i, j] += 1
        end
    end
end

function dd_stride(arr)
    n = size(arr, 1)
    for i = 1:2:n
        for j = 1:2:n
            arr[i, j] += 1
        end
    end
end

function poly_interloop(arr)
    n = size(arr, 1)
    @poly_loop for i = 1:n
        for j = i:n
            arr[i, j] = 1
        end
    end
end

function interloop(arr)
    n = size(arr, 1)
    for i = 1:n
        for j = i:n
            arr[i, j] = 1
        end
    end
end

function poly_interloop_nonuniform(arr)
    n = size(arr, 1)
    @poly_loop for i = 1:n
        for j = 1:n
            for k = i:n
                arr[i, j] += arr[j, k] + k
            end
        end
    end
end

function interloop_nonuniform(arr)
    n = size(arr, 1)
    for i = 1:n
        for j = 1:n
            for k = i:n
                arr[i, j] += arr[j, k] + k
            end
        end
    end
end

function poly_nonuniform(arr)
    n = size(arr, 1)
    @poly_loop for i = 1:n
        for j = 3:n
            arr[i, j] = arr[j, i] + arr[i, j-1]
        end
    end
end

function nonuniform(arr)
    n = size(arr, 1)
    @poly_loop for i = 1:n
        for j = 3:n
            arr[i, j] = arr[j, i] + arr[i, j-1]
        end
    end
end

function poly_lu(PA, L, U)
    n = size(PA, 1)
    @poly_loop for j=1:n
        L[j, j] = 1.0
        for i=1:j
            s1 = 0.0
            for k=1:i
                s1 += U[k, j] * L[i, k]
            end
            U[i, j] = PA[i, j] - s1
        end

        for ii=j+1:n
            s2 = 0.0
            for kk=1:j
                s2 += U[kk, j] * L[ii, kk]
            end
            L[ii, j] = (PA[ii, j] - s2) / U[j, j]
        end
    end
end

function lu(PA, L, U)
    n = size(PA, 1)
    for j=1:n
        L[j, j] = 1.0
        for i=1:j
            s1 = 0.0
            for k=1:i
                s1 += U[k, j] * L[i, k]
            end
            U[i, j] = PA[i, j] - s1
        end

        for ii=j+1:n
            s2 = 0.0
            for kk=1:j
                s2 += U[kk, j] * L[ii, kk]
            end
            L[ii, j] = (PA[ii, j] - s2) / U[j, j]
        end
    end
end

@testset "performance tests" begin

    @testset "automatic loop reordering" begin
        # test matrix mul wrong loop order (check automatic reordering)
        DIM = 512
        A = rand(DIM, DIM)
        B = rand(DIM, DIM)
        out = zeros(DIM, DIM)

        poly_mul(A, B, out) # allow for compiling once

        t_poly = @benchmark poly_mul($A, $B, out) setup=(out = zeros($DIM, $DIM))
        t_orig = @benchmark mul($A, $B, out) setup=(out = zeros($DIM, $DIM))
        t_right_order = @benchmark mul_right_order($A, $B, out) setup=(out = zeros($DIM, $DIM))

        jorig = judge(minimum(t_poly), minimum(t_orig))
        jright = judge(minimum(t_poly), minimum(t_right_order))

        @test jorig.time == :improvement # want an improvement over bad order
        @test jright.time != :regression # want comparable or better than right order
    end

    # @testset "automatic tiling" begin
    #     # test performance nears tiled matrix multiplication (check automatic tiling)
    #     DIM = 1024
    #     A = rand(DIM, DIM)
    #     B = rand(DIM, DIM)
    #     out = zeros(DIM, DIM)
    #
    #     poly_mul(A, B, out) # allow for compiling once
    #
    #     t_poly = @benchmark poly_mul($A, $B, out) setup=(out = zeros($DIM, $DIM))
    #     t_tile = @benchmark tiled_mul($A, $B, out) setup=(out = zeros($DIM, $DIM))
    #
    #     j = judge(minimum(t_poly), minimum(t_tile))
    #
    #     @test j.time != :regression # want improvement or invariant (not regression)
    # end

    """
    All of the following tests check that code never gets slower than original
    """

    @testset "array double" begin
        DIM = 512
        A = rand(DIM)
        out = zeros(DIM)

        poly_array_double(A, out) # allow for compiling once

        t_poly = @benchmark poly_array_double($A, out) setup=(out = zeros($DIM))
        t_orig = @benchmark array_double($A, out) setup=(out = zeros($DIM))

        j = judge(minimum(t_poly), minimum(t_orig))

        @test j.time != :regression
    end

    @testset "matmul double" begin
        DIM = 512
        A = rand(DIM, DIM)
        B = rand(DIM, DIM)
        out = zeros(DIM, DIM)

        poly_matmul_double(A, B, out) # allow for compiling once

        t_poly = @benchmark poly_matmul_double($A, $B, out) setup=(out = zeros($DIM, $DIM))
        t_orig = @benchmark matmul_double($A, $B, out) setup=(out = zeros($DIM, $DIM))

        j = judge(minimum(t_poly), minimum(t_orig))

        @test j.time != :regression
    end

    @testset "2D stride test" begin
        DIM = 512
        arr = ones(DIM, DIM)

        poly_dd_stride(arr) # allow for compiling once

        t_poly = @benchmark poly_dd_stride(arr) setup=(arr = ones($DIM, $DIM))
        t_orig = @benchmark dd_stride(arr) setup=(arr = ones($DIM, $DIM))

        j = judge(minimum(t_poly), minimum(t_orig))

        @test j.time != :regression

    end

    @testset "inter-loop dependence test" begin
        DIM = 512
        arr = zeros(DIM, DIM)

        poly_interloop(arr) # allow for compiling once

        t_poly = @benchmark poly_interloop(arr) setup=(arr = zeros($DIM, $DIM))
        t_orig = @benchmark interloop(arr) setup=(arr = zeros($DIM, $DIM))

        j = judge(minimum(t_poly), minimum(t_orig))

        @test j.time != :regression
    end

    # @testset "inter-loop non-uniform dependence test" begin
    #     DIM = 512
    #     arr = zeros(DIM, DIM)
    #
    #     poly_interloop_nonuniform(arr) # allow for compiling once
    #
    #     t_poly = @benchmark poly_interloop_nonuniform(arr) setup=(arr = zeros($DIM, $DIM))
    #     t_orig = @benchmark interloop_nonuniform(arr) setup=(arr = zeros($DIM, $DIM))
    #
    #     j = judge(minimum(t_poly), minimum(t_orig))
    #
    #     @test j.time != :regression
    # end

    @testset "non-uniform dependence" begin
        DIM = 128
        arr = ones(DIM, DIM)

        poly_nonuniform(arr) # allow for compiling once

        t_poly = @benchmark poly_nonuniform(arr) setup=(arr = ones($DIM, $DIM))
        t_orig = @benchmark nonuniform(arr) setup=(arr = ones($DIM, $DIM))

        j = judge(minimum(t_poly), minimum(t_orig))

        @test j.time != :regression
    end

    @testset "lu decomposition" begin
        DIM = 512
        A = rand(DIM, DIM)*10
        P = zeros(DIM, DIM)
        for i = 1:DIM
            P[i, i] = 1.0
        end
        for j=1:DIM
            i = findmax(A[:, j])[2]
            if j != i
                P[j, :], P[i, :] = P[i, :], P[j, :]
            end
        end

        PA = P*A
        L = zeros(DIM, DIM)
        U = zeros(DIM, DIM)

        poly_lu(PA, L, U) # allow for compiling once

        t_poly = @benchmark poly_lu($PA, L, U) setup=(L = zeros($DIM, $DIM); U = zeros($DIM, $DIM))
        t_orig = @benchmark lu($PA, L, U) setup=(L = zeros($DIM, $DIM); U = zeros($DIM, $DIM))

        j = judge(minimum(t_poly), minimum(t_orig))

        @test j.time != :regression
    end

end
