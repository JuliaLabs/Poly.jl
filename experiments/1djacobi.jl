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

function jacobi_optimized_thread(a::Array{T,1}, b::Array{T,1}, timesteps, n) where {T}
    @threads for t = 1:timesteps-1
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

function jacobi_poly_thread(a::Array{T,1}, b::Array{T,1}, timesteps, n) where {T}
    @poly_loop thread=true tile=0 for t = 1:timesteps-1
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
if thread
    optimized = @benchmark jacobi_optimized_thread(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))
    poly = @benchmark jacobi_poly_thread(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))
    polyrt = nothing
else
    optimized = @benchmark jacobi_optimized(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))
    poly = @benchmark jacobi_poly(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))
    polyrt = @benchmark jacobi_poly_rt(in, B, $t, $n) setup=(in=copy(A); B=zeros(n))
end

@show allocs(basic)
@show allocs(simple)
@show allocs(optimized)
@show allocs(poly)
if !thread
    @show allocs(polyrt)
end
@show minimum(basic)
@show minimum(simple)
@show minimum(optimized)
@show minimum(poly)
if !thread
    @show minimum(polyrt)
end

rbasic = ratio(minimum(basic), minimum(poly))
rsimple = ratio(minimum(simple), minimum(poly))
roptimized = ratio(minimum(optimized), minimum(poly))

if !thread
    rbasicrt = ratio(minimum(basic), minimum(polyrt))
    rsimplert = ratio(minimum(simple), minimum(polyrt))
    roptimizedrt = ratio(minimum(optimized), minimum(polyrt))
end

@show rbasic
@show rsimple
@show roptimized

if !thread
    @show rbasicrt
    @show rsimplert
    @show roptimizedrt
end

jbasic = judge(minimum(poly), minimum(basic))
jsimple = judge(minimum(poly), minimum(simple))
joptimized = judge(minimum(poly), minimum(optimized))

if !thread
    jbasicrt = judge(minimum(polyrt), minimum(basic))
    jsimplert = judge(minimum(polyrt), minimum(simple))
    joptimizedrt = judge(minimum(polyrt), minimum(optimized))
end

@show jbasic
@show jsimple
@show joptimized

if !thread
    @show jbasicrt
    @show jsimplert
    @show joptimizedrt
end
