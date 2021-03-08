
"""
Represents a kernel
Args:
    instructions: list of instructions in kernel. Includes assignments, variables, function calls, etc
    domains: list of ISL SimpleSet domains
    args: arguments to kernel
    argtypes: types of arguments to kernel

Example LoopKernel:
>>> LoopKernel(instructions, domains, [:out, :A, :B], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}])
"""
struct LoopKernel
    instructions::Vector{Instruction}
    domains::Vector{Domain}
    args::Vector{Symbol}
    argtypes::Vector{Type}
end


"""
Represents an ISL SimpleSet
Args:
    iname: symbol representing the simple set
    lowerbound: lowerbound of the simple set
    upperbound: upperbound of the simple set
    recurrence: recurrence relation of elements in the set (i.e, change to loop bounds)

Example Domain:
>>> domains = [Domain(:i, 1, :(size(out, 1)), :(i += 1)),
               Domain(:j, 1, :(size(out, 2)), :(j += 1)),
               Domain(:k, 1, :(size(A, 2)), :(k += 1))]
"""
struct Domain
    iname::Symbol
    lowerbound::Union{Number, Expr}
    upperbound::Union{Number, Expr}
    recurrence::Expr
end


"""
Represents an instruction
Args:
    id: symbol representing instruction
    body: expression in the instruction
    dependencies: ids of instructions that this instruction depends on

Example Instruction:
>>> instructions = [Instruction(:mult, :(out[i, j] += A[i, k] * B[k, j])), Instruction(:double, :(out[i, j] *= 2), dependencies=[:mult])]
"""
struct Instruction
    id::Symbol
    body::Expr
    dependencies::Vector{Symbol}
end

"""
@poly_loop for ()....
    @poly_loop for ()...
        @poly_loop for ()...
            out[i, j] += A[i, k] * B[k, j]
        out[i, j] *= 2

"""
