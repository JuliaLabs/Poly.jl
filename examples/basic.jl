using JuLoop

# multiply arrays and double
instructions = [Instruction(:mult, :(out[i, j] += A[i, k] * B[k, j]), Set()),
                Instruction(:double, :(out[i, j] *= 2), Set())]

domains = [Domain(:i, 1, :(size(out, 1)), :(i += 1), Set(), []),
           Domain(:j, 1, :(size(out, 2)), :(j += 1), Set(), []),
           Domain(:k, 1, :(size(A, 2)), :(k += 1), Set(), [])]


kern = compile_native_julia(LoopKernel(instructions, domains, [:out, :A, :B], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}]))

A = rand(10, 10)
B = rand(10, 10)
out = zeros(10, 10)

kern(out=out, A=A, B=B)

@show isapprox(out, A*B*2)


# add arrays with global offset
instructions = [Instruction(:offset, :(offset = 4), Set()),
                Instruction(:add, :(out[i, j] = A[i, j] + B[i, j] + offset), Set())]

domains = [Domain(:i, 1, :(size(out, 1)), :(i += 1), Set(), []),
           Domain(:j, 1, :(size(out, 2)), :(j += 1), Set(), [])]


kern = compile_native_julia(LoopKernel(instructions, domains, [:out, :A, :B], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}]))

A = rand(10, 10)
B = rand(10, 10)
out = zeros(10, 10)

kern(out=out, A=A, B=B)

@show isapprox(out, A+B.+4)
