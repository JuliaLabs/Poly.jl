using JuLoop
using Test
using StaticArrays

@testset "JuLoop.jl" begin
    # multiply arrays and double
    instructions = [Instruction(:mult, :(out[i, j] += A[i, k] * B[k, j]), Set()),
                    Instruction(:double, :(out[i, j] *= 2), Set())]

    domains = [Domain(:i, 1, :(size(out, 1)), :(i += 1), Set(), []),
               Domain(:j, 1, :(size(out, 2)), :(j += 1), Set(), []),
               Domain(:k, 1, :(size(A, 2)), :(k += 1), Set(), [])]


    kern = compile(LoopKernel(instructions, domains, [:out, :A, :B], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}]))

    A = rand(10, 10)
    B = rand(10, 10)
    out = zeros(10, 10)

    kern(out=out, A=A, B=B)

    @test isapprox(out, A*B*2)


    # add arrays with global offset
    instructions = [Instruction(:offset, :(offset = 4), Set()),
                    Instruction(:add, :(out[i, j] = A[i, j] + B[i, j] + offset), Set())]

    domains = [Domain(:i, 1, :(size(out, 1)), :(i += 1), Set(), []),
               Domain(:j, 1, :(size(out, 2)), :(j += 1), Set(), [])]


    kern = compile(LoopKernel(instructions, domains, [:out, :A, :B], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}]))

    A = rand(10, 10)
    B = rand(10, 10)
    out = zeros(10, 10)

    kern(out=out, A=A, B=B)

    @test isapprox(out, A+B.+4)


    # test macros simple
    A = rand(10)
    out = zeros(10)

    @poly_loop for i = 1:size(out, 1)
        out[i] = A[i]*2
    end

    @test isapprox(out, A*2)

    # test macros nested loops
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

    @test isapprox(out, A*B*2)

    # test macros complicated (tiled matrix multiplication)
    A = rand(128, 128)
    B = rand(128, 128)
    C = zeros(128, 128)

    N = size(C)[1]
    TILE_DIM = 32
    NUM_TILES = Int(N/TILE_DIM)

    @poly_loop for gj = 1:NUM_TILES
        for gi = 1:NUM_TILES
            # loop over tiles needed for this calculation
            tile1 = @MArray zeros(TILE_DIM, TILE_DIM)
            tile2 = @MArray zeros(TILE_DIM, TILE_DIM)
            for t = 0:NUM_TILES-1
                # load tiles needed for calculation
                for i = 1:TILE_DIM
                    for j = 1:TILE_DIM
                        # global tile
                        I = (gi-1) * TILE_DIM + i
                        J = (gj-1) * TILE_DIM + j
                        # get tile1 and tile2 values
                        tile1[i, j] = A[I, t*TILE_DIM + j]
                        tile2[i, j] = B[t*TILE_DIM + i, J]
                    end
                end
                # synchronize
                # loop over tiles to calculate for I, J spot
                for jj in 1:TILE_DIM
                    JJ = (gj-1) * TILE_DIM + jj
                    # loop over row/col in tiles
                    for k = 1:TILE_DIM
                        for ii = 1:TILE_DIM
                            # global tile
                            II = (gi-1) * TILE_DIM + ii
                            # add tile1 * tile2
                            @depends_on elem=t C[II, JJ] += tile1[ii, k] * tile2[k, jj]
                        end
                    end
                end
            end
        end
    end

    @test isapprox(C, A*B)


    # simple @depends_on test
    count = 0
    @poly_loop for i = 1:10
        @depends_on elem=i count += 1
    end

    @test count == 10

end
