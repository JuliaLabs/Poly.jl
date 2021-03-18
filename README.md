# JuLoop

[![Build Status](https://github.com/mjulian31/JuLoop.jl/workflows/CI/badge.svg)](https://github.com/mjulian31/JuLoop.jl/actions)
[![Coverage](https://codecov.io/gh/mjulian31/JuLoop.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mjulian31/JuLoop.jl)

## Code Generation
JuLoop can reorder loops, change iteration spaces, reorder instructions, and restructure code. The main goal is to allow for vectorization and efficient loop execution

### Macros
@poly_loop can be used to tag a loop for code restructuring by JuLoop. Only the outermost loop should use the macro.
@depends_on marks an instruction dependency in order to declare explicit dependencies for code restructuring.
