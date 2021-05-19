using Poly
using LinearAlgebra
using BenchmarkTools
using Base.Threads

num_threads = Base.Threads.nthreads()
LinearAlgebra.BLAS.set_num_threads(num_threads)
thread = false
if num_threads != 1
    thread = true
end


function lu_naive(PA::Array{T,2}, L::Array{T,2}, U::Array{T,2}) where {T}
    n = size(PA, 1)
    for k=1:n
        for i=k+1:n
            PA[i, k] /= PA[k, k]
        end
        for ii=k+1:n
            for j=k+1:n
                PA[ii, j] -= PA[ii,k]*PA[k,j]
            end
        end
    end

    for j=1:n
        for i=1:j
            U[i, j] = PA[i, j]
        end
    end

    for j=1:n
        L[j, j] = 1.0
        for i=j+1:n
            L[i, j] = PA[i, j]
        end
    end
end

function lu_simple(PA::Array{T,2}, L::Array{T,2}, U::Array{T,2}) where {T}
    n = size(PA, 1)
    @inbounds for k=1:n
        for i=k+1:n
            PA[i, k] /= PA[k, k]
        end
        for j=k+1:n
            for ii=k+1:n
                PA[ii, j] -= PA[ii,k]*PA[k,j]
            end
        end
    end

    @inbounds for j=1:n
        for i=1:j
            U[i, j] = PA[i, j]
        end
    end

    @inbounds for j=1:n
        L[j, j] = 1.0
        for i=j+1:n
            L[i, j] = PA[i, j]
        end
    end
end

function lu_poly(PA::Array{T,2}, L::Array{T,2}, U::Array{T,2}) where {T}
    n = size(PA, 1)
    @poly_loop for k=1:n
        for i=k+1:n
            PA[i, k] /= PA[k, k]
        end
        for ii=k+1:n
            for j=k+1:n
                PA[ii, j] -= PA[ii,k]*PA[k,j]
            end
        end
    end

    @poly_loop for j=1:n
        for i=1:j
            U[i, j] = PA[i, j]
        end
    end

    @poly_loop for j=1:n
        L[j, j] = 1.0
        for i=j+1:n
            L[i, j] = PA[i, j]
        end
    end
end

function lu_poly_rt(PA::Array{T,2}, L::Array{T,2}, U::Array{T,2}) where {T}
    n = size(PA, 1)
    @poly_loop for k=1:$n
        for i=k+1:$n
            PA[i, k] /= PA[k, k]
        end
        for ii=k+1:$n
            for j=k+1:$n
                PA[ii, j] -= PA[ii,k]*PA[k,j]
            end
        end
    end

    @poly_loop for j=1:$n
        for i=1:j
            U[i, j] = PA[i, j]
        end
    end

    @poly_loop for j=1:$n
        L[j, j] = 1.0
        for i=j+1:$n
            L[i, j] = PA[i, j]
        end
    end
end

function lu_poly_thread(PA::Array{T,2}, L::Array{T,2}, U::Array{T,2}) where {T}
    n = size(PA, 1)
    @poly_loop thread=true for k=1:n
        for i=k+1:n
            PA[i, k] /= PA[k, k]
        end
        for ii=k+1:n
            for j=k+1:n
                PA[ii, j] -= PA[ii,k]*PA[k,j]
            end
        end
    end

    @poly_loop thread=true for j=1:n
        for i=1:j
            U[i, j] = PA[i, j]
        end
    end

    @poly_loop thread=true for j=1:n
        L[j, j] = 1.0
        for i=j+1:n
            L[i, j] = PA[i, j]
        end
    end
end

dim = 128
A = rand(dim, dim)*10
P = zeros(dim, dim)
for i = 1:dim
    P[i, i] = 1.0
end
for j=1:dim
    i = findmax(A[:, j])[2]
    if j != i
        P[j, :], P[i, :] = P[i, :], P[j, :]
    end
end

PA = P*A
L = zeros(dim, dim)
U = zeros(dim, dim)

basic = @benchmark lu_naive(in, L, U) setup=(in = copy(PA); L = zeros($dim, $dim); U = zeros($dim, $dim))
simple = @benchmark lu_simple(in, L, U) setup=(in = copy(PA); L = zeros($dim, $dim); U = zeros($dim, $dim))
expert = @benchmark LinearAlgebra.lu!(in,  Val(false)) setup=(in = copy(PA))
if thread
    poly = @benchmark lu_poly_thread(in, L, U) setup=(in = copy(PA); L = zeros($dim, $dim); U = zeros($dim, $dim))
else
    poly = @benchmark lu_poly(in, L, U) setup=(in = copy(PA); L = zeros($dim, $dim); U = zeros($dim, $dim))
    polyrt = @benchmark lu_poly_rt(in, L, U) setup=(in = copy(PA); L = zeros($dim, $dim); U = zeros($dim, $dim))
end

@show allocs(basic)
@show allocs(simple)
@show allocs(expert)
@show allocs(poly)
if !thread
    @show allocs(polyrt)
end
@show minimum(basic)
@show minimum(simple)
@show minimum(expert)
@show minimum(poly)
if !thread
    @show minimum(polyrt)
end

rbasic = ratio(minimum(basic), minimum(poly))
rsimple = ratio(minimum(simple), minimum(poly))
rexpert = ratio(minimum(expert), minimum(poly))

if !thread
    rbasicrt = ratio(minimum(basic), minimum(polyrt))
    rsimplert = ratio(minimum(simple), minimum(polyrt))
    rexpertrt = ratio(minimum(expert), minimum(polyrt))
end

@show rbasic
@show rsimple
@show rexpert

if !thread
    @show rbasicrt
    @show rsimplert
    @show rexpertrt
end

jbasic = judge(minimum(poly), minimum(basic))
jsimple = judge(minimum(poly), minimum(simple))
jexpert = judge(minimum(poly), minimum(expert))

if !thread
    jbasicrt = judge(minimum(polyrt), minimum(basic))
    jsimplert = judge(minimum(polyrt), minimum(simple))
    jexpertrt = judge(minimum(polyrt), minimum(expert))
end

@show jbasic
@show jsimple
@show jexpert

if !thread
    @show jbasicrt
    @show jsimplert
    @show jexpertrt
end
