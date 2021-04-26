using JuLoop
using Test
using StaticArrays

@testset "JuLoop.jl" begin
    # multiply arrays and double
    instructions = [Instruction(:mult, :(out[i, j] += A[i, k] * B[k, j]), Set()),
                    Instruction(:double, :(out[i, j] *= 2), Set())]

    domains = [Domain(:i, 1, :(n), :(i += 1), Set(), []),
               Domain(:j, 1, :(m), :(j += 1), Set(), []),
               Domain(:k, 1, :(r), :(k += 1), Set(), [])]

    kern = compile(LoopKernel(instructions, domains, [:out, :A, :B, :n, :r, :m], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}], []))

    A = rand(10, 10)
    B = rand(10, 10)
    out = zeros(10, 10)
    n = size(out, 1)
    r = size(A, 2)
    m = size(out, 2)

    kern(out=out, A=A, B=B, n=n, r=r, m=m)

    @test isapprox(out, A*B*2)


    # add arrays with global offset
    offset = 4
    instructions = [Instruction(:add, :(out[i, j] = A[i, j] + B[i, j] + offset), Set())]

    domains = [Domain(:i, 1, :(n), :(i += 1), Set(), []),
               Domain(:j, 1, :(m), :(j += 1), Set(), [])]


    kern = compile(LoopKernel(instructions, domains, [:out, :A, :B, :n, :m, :offset], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}], []))

    A = rand(10, 10)
    B = rand(10, 10)
    out = zeros(10, 10)
    n = size(out, 1)
    m = size(out, 2)

    kern(out=out, A=A, B=B, n=n, m=m, offset=offset)

    @test isapprox(out, A+B.+4)


    # test macros simple
    A = rand(10)
    out = zeros(10)
    n = size(out, 1)

    @poly_loop for i = 1:n
        out[i] = A[i]*2
    end

    @test isapprox(out, A*2)

    # test macros nested loops
    A = rand(10, 10)
    B = rand(10, 10)
    out = zeros(10, 10)
    n = size(out, 1)
    m = size(out, 2)
    r = size(A, 2)

    @poly_loop for i = 1:n
        for j = 1:m
            for k = 1:r
                out[i, j] += A[i, k] * B[k, j]
            end
            out[i, j] *= 2
        end
    end

    @test isapprox(out, A*B*2)

    # test loop iterators depending on other iterators
    arr = zeros(3, 3)
    expected = [1 1 1;
                0 1 1;
                0 0 1]
    @poly_loop for i = 1:3
        for j = i:3
            arr[i, j] = 1
        end
    end

    @test arr == expected

    # simple @depends_on test
    count = 0
    @poly_loop for i = 1:10
        @depends_on elem=i count += 1
    end

    @test count == 10

    # test macros complicated (tiled matrix multiplication)
    # A = rand(128, 128)
    # B = rand(128, 128)
    # C = zeros(128, 128)
    #
    # N = size(C)[1]
    # TILE_DIM = 32
    # NUM_TILES = Int(N/TILE_DIM)
    #
    # @poly_loop for gj = 1:NUM_TILES
    #     for gi = 1:NUM_TILES
    #         # loop over tiles needed for this calculation
    #         tile1 = @MArray zeros(TILE_DIM, TILE_DIM)
    #         tile2 = @MArray zeros(TILE_DIM, TILE_DIM)
    #         for t = 0:NUM_TILES-1
    #             # load tiles needed for calculation
    #             for i = 1:TILE_DIM
    #                 for j = 1:TILE_DIM
    #                     # global tile
    #                     I = (gi-1) * TILE_DIM + i
    #                     J = (gj-1) * TILE_DIM + j
    #                     # get tile1 and tile2 values
    #                     tile1[i, j] = A[I, t*TILE_DIM + j]
    #                     tile2[i, j] = B[t*TILE_DIM + i, J]
    #                 end
    #             end
    #             # synchronize
    #             # loop over tiles to calculate for I, J spot
    #             for jj in 1:TILE_DIM
    #                 JJ = (gj-1) * TILE_DIM + jj
    #                 # loop over row/col in tiles
    #                 for k = 1:TILE_DIM
    #                     for ii = 1:TILE_DIM
    #                         # global tile
    #                         II = (gi-1) * TILE_DIM + ii
    #                         # add tile1 * tile2
    #                         @depends_on elem=t C[II, JJ] += tile1[ii, k] * tile2[k, jj]
    #                     end
    #                 end
    #             end
    #         end
    #     end
    # end
    #
    # @test isapprox(C, A*B)

end
