using Poly
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using SIMD
using Base.Threads

LinearAlgebra.BLAS.set_num_threads(Base.Threads.nthreads())

function jacobi_naive(a::Array{T,1}, b::Array{T,1}, timesteps, n) where {T}
    for t = 1:timesteps-1
        for i = 3:n-2
            b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        for ii = 3:n-2
            a[ii] = b[ii]
        end
    end
end

function jacobi_simple(a::Array{T,1}, b::Array{T,1}, timesteps, n) where {T}
    @inbounds for t = 1:timesteps-1
        for i = 3:n-2
            b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        for ii = 3:n-2
            a[ii] = b[ii]
        end
    end
end

function jacobi_optimized(a::Array{T,1}, b::Array{T,1}, timesteps, n) where {T}
    @simd for t = 1:timesteps-1
        @simd for i = 3:n-2
            @inbounds b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        @simd for ii = 3:n-2
            @inbounds a[ii] = b[ii]
        end
    end
end

function jacobi_poly(a::Array{T,1}, b::Array{T,1}, timesteps, n) where {T}
    @poly_loop tile=0 for t = 1:timesteps-1
        for i = 3:n-2
            b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        for ii = 3:n-2
            a[ii] = b[ii]
        end
    end
end

function jacobi_poly_rt(a::Array{T,1}, b::Array{T,1}, timesteps, n) where {T}
    @poly_loop tile=0 for t = 1:$timesteps-1
        for i = 3:$n-2
            b[i] = (a[i - 1] + a[i] + a[i + 1])/3
        end
        for ii = 3:$n-2
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
polyrt = @benchmark jacobi_poly_rt(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))

@show allocs(basic)
@show allocs(simple)
@show allocs(optimized)
@show allocs(poly)
@show allocs(polyrt)
@show minimum(basic)
@show minimum(simple)
@show minimum(optimized)
@show minimum(poly)
@show minimum(polyrt)

rbasic = ratio(minimum(basic), minimum(poly))
rsimple = ratio(minimum(simple), minimum(poly))
roptimized = ratio(minimum(optimized), minimum(poly))

rbasicrt = ratio(minimum(basic), minimum(polyrt))
rsimplert = ratio(minimum(simple), minimum(polyrt))
roptimizedrt = ratio(minimum(optimized), minimum(polyrt))

@show rbasic
@show rsimple
@show roptimized

@show rbasicrt
@show rsimplert
@show roptimizedrt

jbasic = judge(minimum(poly), minimum(basic))
jsimple = judge(minimum(poly), minimum(simple))
joptimized = judge(minimum(poly), minimum(optimized))

jbasicrt = judge(minimum(polyrt), minimum(basic))
jsimplert = judge(minimum(polyrt), minimum(simple))
joptimizedrt = judge(minimum(polyrt), minimum(optimized))

@show jbasic
@show jsimple
@show joptimized

@show jbasicrt
@show jsimplert
@show joptimizedrt
