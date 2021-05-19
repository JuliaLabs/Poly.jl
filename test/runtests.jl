using Poly
using Test
using StaticArrays

@testset "correctness tests" begin
    @testset "basic non-macro test" begin
        # multiply arrays and double
        instructions = [Instruction(:mult, :(out[i, j] += A[i, k] * B[k, j]), Set([:i, :j, :k]), :()),
                        Instruction(:double, :(out[i, j] *= 2), Set([:i, :j]), :())]

        dom_k = Domain(:k, 1, :(r), 1, Set([:j]), [instructions[1]])
        dom_j = Domain(:j, 1, :(m), 1, Set([:i]), [dom_k, instructions[2]])
        dom_i = Domain(:i, 1, :(n), 1, Set(), [dom_j])
        domains = [dom_i, dom_j, dom_k]

        kern = compile(LoopKernel(instructions, domains, [:out, :A, :B, :n, :r, :m], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}], []))

        A = rand(10, 10)
        B = rand(10, 10)
        out = zeros(10, 10)
        n = size(out, 1)
        r = size(A, 2)
        m = size(out, 2)

        kern(out=out, A=A, B=B, n=n, r=r, m=m)

        @test isapprox(out, A*B*2)
    end

    @testset "basic non-macro test 2" begin
        # add arrays with global offset
        offset = 4
        instructions = [Instruction(:add, :(out[i, j] = A[i, j] + B[i, j] + offset), Set([:i, :j]), :())]

        dom_j = Domain(:j, 1, :(m), 1, Set([:i]), [instructions[1]])
        dom_i = Domain(:i, 1, :(n), 1, Set(), [dom_j])
        domains = [dom_i, dom_j]


        kern = compile(LoopKernel(instructions, domains, [:out, :A, :B, :n, :m, :offset], [Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}], []))

        A = rand(10, 10)
        B = rand(10, 10)
        out = zeros(10, 10)
        n = size(out, 1)
        m = size(out, 2)

        kern(out=out, A=A, B=B, n=n, m=m, offset=offset)

        @test isapprox(out, A+B.+4)
    end

    @testset "simple macro test" begin
        # test macros simple
        A = rand(10)
        out = zeros(10)
        n = size(out, 1)

        @poly_loop for i = 1:n
            out[i] = A[i]*2
        end

        @test isapprox(out, A*2)
    end

    @testset "nested macro test" begin
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
    end

    @testset "macro stride test" begin
        # test loops with non-1 stride
        arr = ones(6)
        @poly_loop for i = 1:3:6
            arr[i] += 1
        end
        @test arr == [2, 1, 1, 2, 1, 1]
    end

    @testset "macro stride 2D test" begin
        # test loops with non-1 stride
        arr = ones(6, 6)
        @poly_loop for i = 1:2:6
            for j = 1:2:6
                arr[i, j] += 1
            end
        end

        arr2 = ones(6, 6)
        for i = 1:2:6
            for j = 1:2:6
                arr2[i, j] += 1
            end
        end

        @test arr == arr2
    end

    @testset "inter-loop dependence test" begin
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
    end

    @testset "inter-loop dependence test 2" begin
        # test loop iterators depending on other iterators
        arr = zeros(3, 3)
        expected = [1 0 0;
                    1 1 0;
                    1 1 1]
        @poly_loop for i = 1:3
            for j = 1:i
                arr[i, j] = 1
            end
        end

        @test arr == expected
    end

    @testset "inter-loop dependence test 3" begin
        # test loop iterators depending on other iterators
        arr = zeros(3, 3)
        expected = [6 6 6;
                    17 5 5;
                    9 8 3]
        @poly_loop for i = 1:3
            for j = 1:3
                for k = i:3
                    arr[i, j] += arr[j, k] + k
                end
            end
        end

        @test arr == expected
    end

    @testset "dependence test" begin
        # simple @depends_on test
        count = 0
        @poly_loop for i = 1:10
            count += 1
        end

        @test count == 10
    end

    @testset "min/max test" begin
        # checks that min and max functions work
        n = 10
        m = 12
        arr = zeros(n, m)

        @poly_loop for i=1:max(n, 8)
            for j=1:min(15, m)
                arr[i, j] = 1
            end
        end

        arr2 = zeros(n, m)

        for i=1:max(n, 8)
            for j=1:min(15, m)
                arr2[i, j] = 1
            end
        end

        @test arr == arr2
    end

    @testset "floor test" begin
        # tests that floor function works
        n = 10
        arr = zeros(n, n)

        @poly_loop for i=1:floor(n/2)
            for j=1:floor(n/2)
                arr[i, j] = 1
            end
        end

        arr2 = zeros(n, n)

        for i=1:Int(floor(n/2))
            for j=1:Int(floor(n/2))
                arr2[i, j] = 1
            end
        end

        @test arr == arr2
    end

    @testset "if block test" begin
        # tests that an if block works
        n = 10
        arr = zeros(n)

        @poly_loop for i=1:n
            if i % 2 == 0
                arr[i] = 1
            end
        end

        arr2 = zeros(n)

        for i=1:n
            if i % 2 == 0
                arr2[i] = 1
            end
        end

        @test arr == arr2
    end

    @testset "if/elseif block test" begin
        # tests that an if/elseif block works
        n = 10
        arr = zeros(n, n)

        @poly_loop for i=1:n
            for j=1:n
                if i % 2 == 0
                    arr[i, j] = 1
                elseif j % 2 == 0
                    arr[i, j] = 2
                end
            end
        end

        arr2 = zeros(n, n)

        for i=1:n
            for j=1:n
                if i % 2 == 0
                    arr2[i, j] = 1
                elseif j % 2 == 0
                    arr2[i, j] = 2
                end
            end
        end

        @test arr == arr2
    end

    @testset "if/elseif/else block test" begin
        # tests that an if/elseif/else block works
        n = 10
        arr = zeros(n, n)

        @poly_loop for i=1:n
            for j=1:n
                if i % 2 == 0
                    arr[i, j] = 1
                elseif (j + n) % 2 == 0
                    arr[i, j] = 2
                else
                    arr[i, j] = 3
                end
            end
        end

        arr2 = zeros(n, n)

        for i=1:n
            for j=1:n
                if i % 2 == 0
                    arr2[i, j] = 1
                elseif (j + n) % 2 == 0
                    arr2[i, j] = 2
                else
                    arr2[i, j] = 3
                end
            end
        end

        @test arr == arr2
    end

    @testset "if/else block test" begin
        # tests that an if/else block works
        n = 10
        arr = zeros(n, n)

        @poly_loop for i=1:n
            for j=1:n
                if n % 2 == 0
                    arr[i, j] = 1
                else
                    arr[i, j] = 2
                end
            end
        end

        arr2 = zeros(n, n)

        for i=1:n
            for j=1:n
                if n % 2 == 0
                    arr2[i, j] = 1
                else
                    arr2[i, j] = 2
                end
            end
        end

        @test arr == arr2
    end

    @testset "interpolation test 1" begin
        # checks that interpolation works
        c = 2
        n = 10
        arr = zeros(n)

        @poly_loop for i=1:$c:n
            arr[i] = 1
        end

        arr2 = zeros(n)
        for i=1:c:n
            arr2[i] = 1
        end

        @test arr == arr2
    end

    @testset "interpolation test 2" begin
        # checks that multiple 2d interpolation works
        c = 2
        d = 3
        n = 10
        m = 9
        arr = zeros(n, m)

        @poly_loop for i=1:$c:$n
            for j=1:$d:$m
                arr[i, j] = 1
            end
        end

        arr2 = zeros(n, m)
        for i=1:c:n
            for j=1:d:m
                arr2[i, j] = 1
            end
        end

        @test arr == arr2
    end

    @testset "non-uniform dependence" begin
        n = 64
        arr = ones(n, n)
        @poly_loop tile=0 for i = 1:n
            for j = 3:n
                arr[i, j] = arr[j, i] + arr[i, j-1]
            end
        end

        arr2 = ones(n, n)
        for i = 1:n
            for j = 3:n
                arr2[i, j] = arr2[j, i] + arr2[i, j-1]
            end
        end

        @test isapprox(arr, arr2)
    end

    @testset "basic matmul" begin
        # test macros basic matmul
        n = 128
        r = 128
        m = 128
        A = rand(n, r)
        B = rand(r, m)
        out = zeros(n, m)
        n = size(out, 1)
        m = size(out, 2)
        r = size(A, 2)

        @poly_loop for i = 1:n
            for j = 1:m
                for k = 1:r
                    out[i, j] += A[i, k] * B[k, j]
                end
            end
        end

        @test isapprox(out, A*B)
    end

    @testset "1d jacobi" begin
        N = 50
        T = 10
        a = rand(N)
        aa = copy(a)
        b = zeros(N)

        @poly_loop tile=0 for t = 1:T-1
            for i = 3:N-2
                b[i] = (a[i - 1] + a[i] + a[i + 1])/3
            end
            for ii = 3:N-2
                a[ii] = b[ii]
            end
        end

        bb = zeros(N)
        for t = 1:T-1
            for i = 3:N-2
                bb[i] = (aa[i - 1] + aa[i] + aa[i + 1])/3
            end
            for ii = 3:N-2
                aa[ii] = bb[ii]
            end
        end

        @test isapprox(a, aa)
    end

    @testset "lu decomposition" begin
        n = 128
        A = rand(n, n)*10
        L = zeros(n, n)
        U = zeros(n, n)
        P = zeros(n, n)
        for i = 1:n
            P[i, i] = 1.0
        end
        for j=1:n
            i = findmax(A[:, j])[2]
            if j != i
                P[j, :], P[i, :] = P[i, :], P[j, :]
            end
        end
        PA = P*A

        @poly_loop for j=1:n
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

        @test isapprox(PA, L*U)
    end

    @testset "tiled matrix multiplication" begin
        # test macros complicated (tiled matrix multiplication)
        N = 128

        A = rand(N, N)
        B = rand(N, N)
        C = zeros(N, N)

        TILE_DIM = 32
        tile1 = @MArray zeros(TILE_DIM, TILE_DIM)
        tile2 = @MArray zeros(TILE_DIM, TILE_DIM)

        @poly_loop tile=0 for gj = 0:$TILE_DIM:N-1
            for gi = 0:$TILE_DIM:N-1
                for t = 0:$TILE_DIM:N-1
                    for j = 1:$TILE_DIM
                        for i = 1:$TILE_DIM
                            tile1[i, j] = A[gi + i, t + j]
                            tile2[i, j] = B[t + i, gj + j]
                        end
                    end
                    for jj in 1:$TILE_DIM
                        for k = 1:$TILE_DIM
                            for ii = 1:$TILE_DIM
                                C[gi + ii, gj + jj] += tile1[ii, k] * tile2[k, jj]
                            end
                        end
                    end
                end
            end
        end

        @test isapprox(C, A*B)
    end

end

# uncomment to run performance testing (takes a while, so commented out for github)
# include("performancetests.jl")
