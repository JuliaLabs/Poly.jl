var documenterSearchIndex = {"docs":
[{"location":"#Poly.jl-Documentation","page":"Poly.jl Documentation","title":"Poly.jl Documentation","text":"","category":"section"},{"location":"","page":"Poly.jl Documentation","title":"Poly.jl Documentation","text":"","category":"page"},{"location":"#Functions","page":"Poly.jl Documentation","title":"Functions","text":"","category":"section"},{"location":"","page":"Poly.jl Documentation","title":"Poly.jl Documentation","text":"@poly_loop(ex0...)","category":"page"},{"location":"#Poly.@poly_loop-Tuple","page":"Poly.jl Documentation","title":"Poly.@poly_loop","text":"Runs a Polyhedral model on the input for loop. May reorder instructions (including loop orderings) and vectorize code. Only the outermost loop should use the macro. By default, will also tile code. See below for more details.\n\nExample\n\n@poly_loop for i = 1:n\n    for j = 1:m\n        for k = 1:r\n            out[i, j] += A[i, k] * B[k, j]\n        end\n        out[i, j] *= 2\n    end\nend\n\n(runs a matrix multiplication and double)\n\nnote: Loop Types\n@poly_loop requires a normal for loop (i.e i = lowerbound:upperbound or i = lowerbound:step:upperbound)\n\nTiling\n\nLoops will be automatically tiled based on l1 cache sizes unless tile=0 is set.\n\nExample No Tiling:\n\n@poly_loop tile=0 for i=1:n\n    arr[i] = 1.0\nend\n\nTiling is usually faster, but in cases where the loop bounds are small it will just introduce more loop overhead. Unless tile=0 is set, the outermost band (grouping of permutable loops) will be tiled.\n\nwarning: Invalid Tiling\nSome non-uniform access patterns can result in invalid tiling. This is pretty rare. If this occurs, code will likely error during compilation.\n\nThreading\n\nTo enable multithreading of the outermost loop, pass thread=true (disabled by default):\n\nExample:\n\n@poly_loop thread=true for i=1:n\n    arr[i] = 1.0\nend\n\nnote: Threading Condition\nJulia must be run with julia --threads=auto (or a number of threads instead of auto) to enable threading at all.\n\nDebugging\n\nDebugging and verbosity options (see runpolyhedralmodel for details) can be passed like:\n\njulia> @poly_loop debug=true verbose=2 for i=1:4:n-1\n    arr[i] *= arr[i + 1]\nend\n\nThis will print out the following information, depending on the verbosity level set:\n\n0: no printing (unless debug=true, where ISL errors are printed)\n1: final Julia expression is printed\n2: initial C code, loop orderings, new schedule, final C code, and final Julia expression printed\n3: all of 2 plus domain, access relations, dependence relations, original schedule, and new schedule constraints printed\n\nFunctions\n\nwarning: Functions\nFunction calls are not currently supported in most cases\n\nIf a function modifies any inputs, then there is no way to know which inputs are modified or how those inputs are modified. If a function returns an output such as a matrix, the write to each individual index of the matrix (for example A[i, j]) cannot be inferrred by ISL. In addition, @inbounds is not neccessary (or allowed) as the aggresive transformation will add @inbounds automatically.\n\nLoop Bounds\n\nThis also applies to loop bounds. One easy workaround is something like:\n\njuila> n = size(out, 1)\njulia> @poly_loop for i = 1:n\n    out[i] += 1.0\nend\n\nExceptions\n\nThe only exceptions are the following functions (which can be used in loop bounds and in instructions): min():\n\n@poly_loop for i=1:min(n, m)\n    arr[i, i] = 1.0\nend\n\nmax():\n\n@poly_loop for i=1:max(n, m)\n    arr[i, i] -= 1.0\nend\n\nfloor() -> DO NOT cast to Int() (will be added later):\n\n@poly_loop for i=1:floor(n/2)\n    arr[i] += 1.0\nend\n\nnote: Loop Names\nAll loop iterators must have unique symbols, even if in different scopes. This is so an original schedule can be extracted from the code.\n\nExample\n\n@poly_loop for i=1:n\n    for j=1:floor(n/2)\n        arr[i, j] += 1.0\n    end\n    for jj=floor(n/2):n\n        arr[i, j] *= 2.0\n    end\nend\n\nInterpolation\n\nwarning: Striding\nISL does not work with non-numerical striding, i.e. for i=1:TILE:n.\n\nAs a workaround, this macro supports interpolating the values into the function so that ISL has access to the numerical values (and non-numerical strides MUST be interpolated). This delays the analysis to runtime, so it is advised that interpolated values are typically constant so that compilation does not occur over and over.\n\nExamples:\n\n@poly_loop for i=1:$c:n\n    arr[i] = 1\nend\n\n@poly_loop for i=1:$c:$n\n    arr[i] = 1\nend\n\nIn addition, any constants can be interpolated, and sometimes this may improve the performance of ISL. Numerical bounds can help with tiling and inference on spaces. However, it is important that interpolated values are actually constants, or very rarely change, so that this delayed compilation does not occur regularly during runtime.\n\nnote: Indexing\nAll array indices must be in terms of loop bounds or constants, not of variables.\n\nExample:\n\njulia> c = i + j\njulia> A[c]\n\nMust become\n\njulia> A[i + j]\n\nwarning: Multiplicative Indexing\nMultilpicative indexing (i.e. A[i*t, j]) is not supported. If needed, stride over the iterator.\n\nExample:\n\nNUM_TILES = Int(N/TILE_SIZE)\n@poly_loop for t=1:NUM_TILES\n    A[TILE_SIZE*t] = 1\nend\n\nWould need to become this:\n\n@poly_loop for t=1:$TILE_SIZE:N\n    A[t] = 1\nend\n\nIf/elseif/else blocks:\n\nnote: If Block Conditions\nIf blocks must contain conditions on the loop iterators and constants only (not on any elements of an array, for example)\n\nExample:\n\n@poly_loop for i=1:n\n    if i%2 == 0\n        A[i] = 1\n    else\n        A[i] = 2\n    end\nend\n\n\n\n\n\n","category":"macro"},{"location":"#Index","page":"Poly.jl Documentation","title":"Index","text":"","category":"section"},{"location":"","page":"Poly.jl Documentation","title":"Poly.jl Documentation","text":"","category":"page"}]
}
