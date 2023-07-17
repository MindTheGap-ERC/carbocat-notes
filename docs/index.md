---
title: CarboKitten
subtitle: project notes
---

# Notes

These are project notes and exercises.

``` {.julia file=src/MindTheGap.jl}
module MindTheGap

include("./Stencil.jl")
include("./Burgess2013.jl")
include("./Utility.jl")
include("./RandomFields.jl")
include("./BS92.jl")

end
```

``` {.julia file=test/runtests.jl}
using Test
using MindTheGap.Stencil
using MindTheGap.Utility

@testset "Mind The Gap tests" begin
    <<spec>>
end
```
