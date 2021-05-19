using Poly
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using SIMD
using Base.Threads

num_threads = Base.Threads.nthreads()
LinearAlgebra.BLAS.set_num_threads(num_threads)
thread = false
if num_threads != 1
    thread = true
end


function mul_naive(A::Array{T,2}, B::LinearAlgebra.Adjoint{T,Array{T,2}}, C::Array{T,2}) where {T}
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

function mul_simple(A::Array{T,2}, B::LinearAlgebra.Adjoint{T,Array{T,2}}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)

    @inbounds for k = 1:R
        for j = 1:M
            for i = 1:N
                C[i, j] += A[i, k] * B[k, j]
            end
        end
    end
end

function mul_optimized(A::Array{T,2}, B::LinearAlgebra.Adjoint{T,Array{T,2}}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)
    TILE_DIM = 64

    @simd for t = 0:TILE_DIM:R-1
        @simd for gj = 0:TILE_DIM:M-1
            @simd for gi = 0:TILE_DIM:N-1
            tile1 = @MArray zeros(TILE_DIM, TILE_DIM)
            tile2 = @MArray zeros(TILE_DIM, TILE_DIM)
                @simd for j = 1:TILE_DIM
                    @simd for i = 1:TILE_DIM
                        @inbounds tile1[i, j] = A[gi + i, t + j]
                        @inbounds tile2[i, j] = B[t + i, gj + j]
                    end
                end
                @simd for k = 1:TILE_DIM
                    @simd for j = 1:TILE_DIM
                        @simd for i = 1:TILE_DIM
                            @inbounds C[gi + i, gj + j] += tile1[i, k] * tile2[k, j]
                        end
                    end
                end
            end
        end
    end
end

function mul_optimized_v2(A::Array{T,2}, B::LinearAlgebra.Adjoint{T,Array{T,2}}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)
    TILE_DIM = 64

    @simd for t = 0:TILE_DIM:R-1
        @simd for gj = 0:TILE_DIM:M-1
            @simd for gi = 0:TILE_DIM:N-1
                @simd for k = 1:TILE_DIM
                    @simd for j = 1:TILE_DIM
                        @simd for i = 1:TILE_DIM
                            @inbounds C[gi + i, gj + j] += A[gi + i, t + k] * B[t + k, gj + j]
                        end
                    end
                end
            end
        end
    end
end

function mul_optimized_thread(A::Array{T,2}, B::LinearAlgebra.Adjoint{T,Array{T,2}}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)
    TILE_DIM = 64

    @threads for t = 0:TILE_DIM:R-1
        @simd for gj = 0:TILE_DIM:M-1
            @simd for gi = 0:TILE_DIM:N-1
                @simd for k = 1:TILE_DIM
                    @simd for j = 1:TILE_DIM
                        @simd for i = 1:TILE_DIM
                            @inbounds C[gi + i, gj + j] += A[gi + i, t + k] * B[t + k, gj + j]
                        end
                    end
                end
            end
        end
    end
end

function mul_poly(A::Array{T,2}, B::LinearAlgebra.Adjoint{T,Array{T,2}}, C::Array{T,2}) where {T}
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

function mul_poly_rt(A::Array{T,2}, B::LinearAlgebra.Adjoint{T,Array{T,2}}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)

    @poly_loop for i = 1:$N
        for j = 1:$M
            for k = 1:$R
                C[i, j] += A[i, k] * B[k, j]
            end
        end
    end
end

function mul_poly_thread(A::Array{T,2}, B::LinearAlgebra.Adjoint{T,Array{T,2}}, C::Array{T,2}) where {T}
    N = size(C, 1)
    M = size(C, 2)
    R = size(A, 2)

    @poly_loop thread=true for i = 1:N
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
B = rand(m, r)
B = B'
out = zeros(n, m)

basic = @benchmark mul_naive($A, $B, out) setup=(out=zeros($n, $m))
simple = @benchmark mul_simple($A, $B, out) setup=(out=zeros($n, $m))
optimized = @benchmark mul_optimized($A, $B, out) setup=(out=zeros($n, $m))
if thread
    optimized2 = @benchmark mul_optimized_thread($A, $B, out) setup=(out=zeros($n, $m))
    poly = @benchmark mul_poly_thread($A, $B, out) setup=(out=zeros($n, $m))
    polyrt = nothing
else
    optimized2 = @benchmark mul_optimized_v2($A, $B, out) setup=(out=zeros($n, $m))
    poly = @benchmark mul_poly($A, $B, out) setup=(out=zeros($n, $m))
    polyrt = @benchmark mul_poly_rt($A, $B, out) setup=(out=zeros($n, $m))
end
expert = @benchmark LinearAlgebra.mul!(out, $A, $B) setup=(out=zeros($n, $m))

@show allocs(basic)
@show allocs(simple)
@show allocs(optimized)
@show allocs(optimized2)
@show allocs(expert)
@show allocs(poly)
if !thread
    @show allocs(polyrt)
end
@show minimum(basic)
@show minimum(simple)
@show minimum(optimized)
@show minimum(optimized2)
@show minimum(expert)
@show minimum(poly)
@show minimum(polyrt)

rbasic = ratio(minimum(basic), minimum(poly))
rsimple = ratio(minimum(simple), minimum(poly))
roptimized = ratio(minimum(optimized), minimum(poly))
roptimized2 = ratio(minimum(optimized2), minimum(poly))
rexpert = ratio(minimum(expert), minimum(poly))

if !thread
    rbasicrt = ratio(minimum(basic), minimum(polyrt))
    rsimplert = ratio(minimum(simple), minimum(polyrt))
    roptimizedrt = ratio(minimum(optimized), minimum(polyrt))
    roptimized2rt = ratio(minimum(optimized2), minimum(polyrt))
    rexpertrt = ratio(minimum(expert), minimum(polyrt))
end

@show rbasic
@show rsimple
@show roptimized
@show roptimized2
@show rexpert

if !thread
    @show rbasicrt
    @show rsimplert
    @show roptimizedrt
    @show roptimized2rt
    @show rexpertrt
end

jbasic = judge(minimum(poly), minimum(basic))
jsimple = judge(minimum(poly), minimum(simple))
joptimized = judge(minimum(poly), minimum(optimized))
joptimized2 = judge(minimum(poly), minimum(optimized2))
jexpert = judge(minimum(poly), minimum(expert))

if !thread
    jbasicrt = judge(minimum(polyrt), minimum(basic))
    jsimplert = judge(minimum(polyrt), minimum(simple))
    joptimizedrt = judge(minimum(polyrt), minimum(optimized))
    joptimized2rt = judge(minimum(polyrt), minimum(optimized2))
    jexpertrt = judge(minimum(polyrt), minimum(expert))
end

@show jbasic
@show jsimple
@show joptimized
@show joptimized2
@show jexpert

if !thread
    @show jbasicrt
    @show jsimplert
    @show joptimizedrt
    @show joptimized2rt
    @show jexpertrt
end
