using JuLoop
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using SIMD

# run with --threads auto for jacobitithreading
if length(ARGS) == 0
    # so BLAS can't cheat
    LinearAlgebra.BLAS.set_num_threads(1)
end


function jacobi_naive(A::Array{T,2}, B::Array{T,2}, timesteps, n) where {T}
    for t = 1:timesteps-1
        for i = 3:n-2
            b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        for ii = 3:n-2
            a[ii] = b[ii]
        end
    end
end

function jacobi_simple(A::Array{T,2}, B::Array{T,2}, timesteps, n) where {T}
    @inbounds for t = 1:timesteps-1
        for i = 3:n-2
            b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        for ii = 3:n-2
            a[ii] = b[ii]
        end
    end
end

function jacobi_optimized(A::Array{T,2}, B::Array{T,2}, timesteps, n) where {T}
    @simd for t = 1:timesteps-1
        @simd for i = 3:n-2
            @inbounds b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        @simd for ii = 3:n-2
            @inbounds a[ii] = b[ii]
        end
    end
end

function jacobi_poly(A::Array{T,2}, B::Array{T,2}, timesteps, n) where {T}
    @poly_loop tile=0 for t = 1:timesteps-1
        for i = 3:n-2
            b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        for ii = 3:n-2
            a[ii] = b[ii]
        end
    end
end

n = 1024
A = rand(n)
B = zeros(n)
t = 512

basic = @benchmark jacobi_naive(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))
simple = @benchmark jacobi_simple(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))
optimized = @benchmark jacobi_optimized(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))
poly = @benchmark jacobi_poly(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))

@show allocs(basic)
@show allocs(simple)
@show allocs(optimized)
@show allocs(poly)
@show minimum(basic)
@show minimum(simple)
@show minimum(optimized)
@show minimum(poly)

rbasic = ratio(minimum(basic), minimum(poly))
rsimple = ratio(minimum(simple), minimum(poly))
roptimized = ratio(minimum(optimized), minimum(poly))

@show rbasic
@show rsimple
@show roptimized

jbasic = judge(minimum(poly), minimum(basic))
jsimple = judge(minimum(poly), minimum(simple))
joptimized = judge(minimum(poly), minimum(optimized))

@show jbasic
@show jsimple
@show joptimized
