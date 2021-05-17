using JuLoop
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using SIMD
# so BLAS can't cheat
LinearAlgebra.BLAS.set_num_threads(1)


function mul_basic(A::Array{T,2}, B::Array{T,2}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)

    for i = 1:N
        for j = 1:M
            for k = 1:R
                C[i, j] += A[i, k] * B[k, j]
            end
        end
    end
end

function mul_optimized(A::Array{T,2}, B::Array{T,2}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)
    TILE_DIM = 64

    @simd for gj = 0:TILE_DIM:M-1
        @simd for t = 0:TILE_DIM:R-1
            @simd for gi = 0:TILE_DIM:N-1
            tile1 = @MArray zeros(TILE_DIM, TILE_DIM)
            tile2 = @MArray zeros(TILE_DIM, TILE_DIM)
                @simd for j = 1:TILE_DIM
                    @simd for i = 1:TILE_DIM
                        @inbounds tile1[i, j] = A[gi + i, t + j]
                        @inbounds tile2[i, j] = B[t + i, gj + j]
                    end
                end
                @simd for j = 1:TILE_DIM
                    @simd for k = 1:TILE_DIM
                        @simd for i = 1:TILE_DIM
                            @inbounds C[gi + i, gj + j] += tile1[i, k] * tile2[k, j]
                        end
                    end
                end
            end
        end
    end
end

function mul_poly(A::Array{T,2}, B::Array{T,2}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)

    @poly_loop for i = 1:N
        for j = 1:M
            for k = 1:R
                C[i, j] += A[i, k] * B[k, j]
            end
        end
    end
end

n = 1024
r = 512
m = 1024
A = rand(n, r)
B = rand(r, m)
out = zeros(n, m)

basic = @benchmark mul_basic($A, $B, out) setup=(out=zeros($n, $m))
optimized = @benchmark mul_optimized($A, $B, out) setup=(out=zeros($n, $m))
expert = @benchmark LinearAlgebra.mul!(out, $A, $B) setup=(out=zeros($n, $m))
poly = @benchmark mul_poly($A, $B, out) setup=(out=zeros($n, $m))

@show basic
@show optimized
@show expert
@show poly

rbasic = ratio(minimum(basic), minimum(poly))
roptimized = ratio(minimum(optimized), minimum(poly))
rexpert = ratio(minimum(expert), minimum(poly))

@show rbasic
@show roptimized
@show rexpert

jbasic = judge(minimum(poly), minimum(basic))
joptimized = judge(minimum(poly), minimum(optimized))
jexpert = judge(minimum(poly), minimum(expert))

@show jbasic
@show joptimized
@show jexpert
