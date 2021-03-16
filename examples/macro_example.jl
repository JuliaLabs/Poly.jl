using JuLoop

A = rand(10)
out = zeros(10)

@poly_loop for i = 1:size(out, 1)
    out[i] = A[i]*2
end

@show isapprox(out, A*2)


A = rand(10, 10)
B = rand(10, 10)
out = zeros(10, 10)

@poly_loop for i = 1:size(out, 1)
    for j = 1:size(out, 2)
        for k = 1:size(A, 2)
            out[i, j] += A[i, k] * B[k, j]
        end
        out[i, j] *= 2
    end
end

@show isapprox(A*B*2, out)
