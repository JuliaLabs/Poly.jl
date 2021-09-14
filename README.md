# Poly.jl

[![Build Status](https://github.com/mjulian31/Poly.jl/workflows/CI/badge.svg)](https://github.com/mjulian31/Poly.jl/actions)
[![Coverage](https://codecov.io/gh/mjulian31/Poly.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mjulian31/Poly.jl)
[![Documentation]][http://julia.mit.edu/Poly.jl/dev/]

## Code Generation
Poly can reorder loops, change iteration spaces, reorder instructions, and restructure code. The main goal is to allow for vectorization and efficient loop execution

### Macros
@poly_loop can be used to tag a loop for code restructuring by Poly, using the polyhedral model of compilation. This can result in automatic loop reordering, domain coalescing, tiling, and more, specifically with memory locality and vectorization in mind. See the documentation for details on how to use @poly_loop.
