module JuLoop

export LoopKernel, Domain, Instruction


"""
Represents an instruction
Args:
    iname: symbol representing instruction
    body: expression in the instruction
    dependencies: ids of instructions that this instruction depends on

Example Instruction:
>>> instructions = [Instruction(:mult, :(out[i, j] += A[i, k] * B[k, j]), Set()),
                    Instruction(:double, :(out[i, j] *= 2), Set([:mult]))]
"""
struct Instruction
    iname::Symbol
    body::Expr
    dependencies::Set{Symbol}
end

Base.show(io::IO, d::Instruction) = print(io, d.iname)


"""
Represents an iteration domain
Args:
    iname: symbol representing the simple set
    lowerbound: lowerbound of the simple set (inclusive)
    upperbound: upperbound of the simple set (inclusive)
    step: step of loop iteration
    dependencies: dependencies of this domain
    instructions: instructions in the body of this domain (assigned after initialization)

Example Domain:
>>> domains = [Domain(:i, 1, :(size(out, 1)), :(i += 1), Set(), []),
               Domain(:j, 1, :(size(out, 2)), :(j += 1), Set(), []),
               Domain(:k, 1, :(size(A, 2)), :(k += 1), Set(), [])]
"""
struct Domain
    iname::Symbol
    lowerbound::Union{Number, Symbol, Expr}
    upperbound::Union{Number, Symbol, Expr}
    step::Union{Number, Symbol, Expr}
    dependencies::Set{Symbol}
    instructions::Vector{Union{Instruction, Domain}}
end

Base.show(io::IO, d::Domain) = print(io, d.iname)


"""
Represents a kernel
Args:
    instructions: list of instructions in kernel. Includes assignments, variables, function calls, etc
    domains: list of ISL SimpleSet domains
    args: arguments to kernel
    argtypes: types of arguments to kernel

Example LoopKernel:
>>> kern = LoopKernel(instructions, domains, [:out, :A, :B], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}], [])
"""
struct LoopKernel
    instructions::Vector{Instruction}
    domains::Vector{Domain}
    args::Vector{Symbol}
    argtypes::Vector{Type}
    consts::Vector{Symbol}
end


"""
Represents the AST of instruction and loop ordering for a kernel
Args:
    self: either Instruction or Domain object (or nothing)
    children: list of dependent Instructions or Domains
    parents: list of Instructions or Domains dependent on
"""
struct AST
    self::Union{Domain, Instruction, Nothing}
    children::Vector{AST}
    parents::Vector{AST}
end

function get_id(ast::AST)
    if ast.self == nothing
        return nothing
    end
    return ast.self.iname
end

using AbstractTrees
AbstractTrees.children(d::AST) = d.children
AbstractTrees.parentlinks(::AST) = AbstractTrees.StoredParents()
AbstractTrees.printnode(io::IO, d::AST) = print(io, get_id(d))
Base.show(io::IO, d::AST) = print_tree(io, d)

include("compilation.jl")
include("macros.jl")
include("polyhedral.jl")
include("striding_analysis.jl")

end
