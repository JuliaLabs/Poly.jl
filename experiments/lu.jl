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

function lu_simple(PA::Array{T,2}, L::Array{T,2}, U::Array{T,2}) where {T}
    n = size(PA, 1)
    @inbounds for j=1:n
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

function lu_poly(PA::Array{T,2}, L::Array{T,2}, U::Array{T,2}) where {T}
    n = size(PA, 1)
    @poly_loop thread=thread for j=1:n
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

function lu_poly_rt(PA::Array{T,2}, L::Array{T,2}, U::Array{T,2}) where {T}
    n = size(PA, 1)
    s1 = 0.0
    s2 = 0.0
    @poly_loop thread=thread for j=1:$n
        L[j, j] = 1.0
        for i=1:j
            s1 = 0.0
            for k=1:i
                s1 += U[k, j] * L[i, k]
            end
            U[i, j] = PA[i, j] - s1
        end

        for ii=j+1:$n
            s2 = 0.0
            for kk=1:j
                s2 += U[kk, j] * L[ii, kk]
            end
            L[ii, j] = (PA[ii, j] - s2) / U[j, j]
        end
    end
end

dim = 1024
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

basic = @benchmark lu_naive($PA, L, U) setup=(L = zeros($dim, $dim); U = zeros($dim, $dim))
simple = @benchmark lu_simple($PA, L, U) setup=(L = zeros($dim, $dim); U = zeros($dim, $dim))
expert = @benchmark LinearAlgebra.lu($PA,  Val(false))
poly = @benchmark lu_poly($PA, L, U) setup=(L = zeros($dim, $dim); U = zeros($dim, $dim))
polyrt = @benchmark lu_poly_rt($PA, L, U) setup=(L = zeros($dim, $dim); U = zeros($dim, $dim))

@show allocs(basic)
@show allocs(simple)
@show allocs(expert)
@show allocs(poly)
@show allocs(polyrt)
@show minimum(basic)
@show minimum(simple)
@show minimum(expert)
@show minimum(poly)
@show minimum(polyrt)

rbasic = ratio(minimum(basic), minimum(poly))
rsimple = ratio(minimum(simple), minimum(poly))
rexpert = ratio(minimum(expert), minimum(poly))

rbasicrt = ratio(minimum(basic), minimum(polyrt))
rsimplert = ratio(minimum(simple), minimum(polyrt))
rexpertrt = ratio(minimum(expert), minimum(polyrt))

@show rbasic
@show rsimple
@show rexpert

@show rbasicrt
@show rsimplert
@show rexpertrt

jbasic = judge(minimum(poly), minimum(basic))
jsimple = judge(minimum(poly), minimum(simple))
jexpert = judge(minimum(poly), minimum(expert))

jbasicrt = judge(minimum(polyrt), minimum(basic))
jsimplert = judge(minimum(polyrt), minimum(simple))
jexpertrt = judge(minimum(polyrt), minimum(expert))

@show jbasic
@show jsimple
@show jexpert

@show jbasicrt
@show jsimplert
@show jexpertrt
