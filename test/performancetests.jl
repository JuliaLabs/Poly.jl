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

function poly_mul_no_tile(A, B, out)
    n = size(out, 1)
    m = size(out, 2)
    r = size(A, 2)
    @poly_loop tile=0 for i = 1:n
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
    TILE_DIM = 64
    NUM_TILES = Int(N/TILE_DIM)
    # tile1 = @MArray zeros(TILE_DIM, TILE_DIM)
    # tile2 = @MArray zeros(TILE_DIM, TILE_DIM)

    @inbounds for gj = 0:TILE_DIM:N-1
        for gk = 0:TILE_DIM:N-1
            for gi = 0:TILE_DIM:N-1
                # for j = 1:TILE_DIM
                #     for i = 1:TILE_DIM
                #         tile1[i, j] = A[gi + i, gk + j]
                #         tile2[i, j] = B[gk + i, gj + j]
                #     end
                # end
                for jj in 1:TILE_DIM
                    for k = 1:TILE_DIM
                        for ii = 1:TILE_DIM
                            # C[gi + ii, gj + jj] += tile1[ii, k] * tile2[k, jj]
                            C[gi + ii, gj + jj] += A[gi + ii, gk + k] * B[gk + k, gj + jj]
                        end
                    end
                end
            end
        end
    end
end

function tiled_mul_notsquare(A, B, C)
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)
    TILE_DIM = 64

    @inbounds for gj = 0:TILE_DIM:M-1
        for gk = 0:TILE_DIM:R-1
            for gi = 0:TILE_DIM:N-1
                for jj in 1:TILE_DIM
                    for k = 1:TILE_DIM
                        for ii = 1:TILE_DIM
                            C[gi + ii, gj + jj] += A[gi + ii, gk + k] * B[gk + k, gj + jj]
                        end
                    end
                end
            end
        end
    end
end

function tiled_mul_any(A, B, C)
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)
    TILE_DIM = 64

    @inbounds for gj = 1:TILE_DIM:M
        for gk = 1:TILE_DIM:R
            for gi = 1:TILE_DIM:N
                for jj in gj:min(gj + TILE_DIM, M)
                    for k = gk:min(gk + TILE_DIM, R)
                        for ii = gi:min(gi + TILE_DIM, N)
                            C[ii, jj] += A[ii, k] * B[k, jj]
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
        for j = 1:n
            out[i, j] = A[i, j]*2
        end
    end
end

function array_double(A, out)
    n = size(out, 1)

    for i = 1:n
        for j = 1:n
            out[i, j] = A[i, j]*2
        end
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

        poly_mul_no_tile(A, B, out) # allow for compiling once

        t_poly = @benchmark poly_mul_no_tile($A, $B, out) setup=(out = zeros($DIM, $DIM))
        t_orig = @benchmark mul($A, $B, out) setup=(out = zeros($DIM, $DIM))
        t_right_order = @benchmark mul_right_order($A, $B, out) setup=(out = zeros($DIM, $DIM))

        jorig = judge(minimum(t_poly), minimum(t_orig))
        jright = judge(minimum(t_poly), minimum(t_right_order))

        @test jorig.time == :improvement # want an improvement over bad order
        @test jright.time != :regression # want comparable or better than right order
    end

    @testset "automatic tiling square" begin
        # test performance nears tiled matrix multiplication (check automatic tiling)
        DIM = 1024
        A = rand(DIM, DIM)
        B = rand(DIM, DIM)
        out = zeros(DIM, DIM)

        poly_mul(A, B, out) # allow for compiling once

        t_poly = @benchmark poly_mul($A, $B, out) setup=(out = zeros($DIM, $DIM))
        t_tile = @benchmark tiled_mul($A, $B, out) setup=(out = zeros($DIM, $DIM))

        j = judge(minimum(t_poly), minimum(t_tile))

        @test j.time != :regression # want improvement or invariant (not regression)
    end

    @testset "automatic tiling non-square" begin
        # test performance nears tiled matrix multiplication (check automatic tiling), for non-square matricies that are multiples of tile dim
        N = 1024
        R = 2048
        M = 512
        A = rand(N, R)
        B = rand(R, M)
        out = zeros(N, M)

        poly_mul(A, B, out) # allow for compiling once

        t_poly = @benchmark poly_mul($A, $B, out) setup=(out = zeros($N, $M))
        t_tile = @benchmark tiled_mul_notsquare($A, $B, out) setup=(out = zeros($N, $M))

        j = judge(minimum(t_poly), minimum(t_tile))

        @test j.time != :regression # want improvement or invariant (not regression)
    end

    @testset "automatic tiling any dimensions" begin
        # test performance nears tiled matrix multiplication (check automatic tiling), for any size matrices
        N = 1020
        R = 900
        M = 981
        A = rand(N, R)
        B = rand(R, M)
        out = zeros(N, M)

        poly_mul(A, B, out) # allow for compiling once

        t_poly = @benchmark poly_mul($A, $B, out) setup=(out = zeros($N, $M))
        t_tile = @benchmark tiled_mul_any($A, $B, out) setup=(out = zeros($N, $M))

        j = judge(minimum(t_poly), minimum(t_tile))

        @test j.time != :regression # want improvement or invariant (not regression)
    end

    """
    All of the following tests check that code always gets faster than the original
    All of the code tested has misordered loops or can improve from tiling
    """

    @testset "array double" begin
        DIM = 512
        A = rand(DIM, DIM)
        out = zeros(DIM, DIM)

        poly_array_double(A, out) # allow for compiling once

        t_poly = @benchmark poly_array_double($A, out) setup=(out = zeros($DIM, $DIM))
        t_orig = @benchmark array_double($A, out) setup=(out = zeros($DIM, $DIM))

        j = judge(minimum(t_poly), minimum(t_orig))

        @test j.time == :improvement
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

        @test j.time == :improvement
    end

    @testset "2D stride" begin
        DIM = 512
        arr = ones(DIM, DIM)

        poly_dd_stride(arr) # allow for compiling once

        t_poly = @benchmark poly_dd_stride(arr) setup=(arr = ones($DIM, $DIM))
        t_orig = @benchmark dd_stride(arr) setup=(arr = ones($DIM, $DIM))

        j = judge(minimum(t_poly), minimum(t_orig))

        @test j.time == :improvement

    end

    @testset "lu decomposition" begin
        DIM = 1024

        # this step is hard in plain julia (no functions), so we just do it ahead of time
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

        @test j.time == :improvement
    end

end
