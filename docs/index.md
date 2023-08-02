---
title: CarboKitten
subtitle: project notes
---

# Index

- [Carbonate Production](./bosscher-1992.md) a reproduction of the Bosscher et al. '92 result
- [CarboKitten (summary)](./carbocat.md) a Julia clone of CarboCAT
  - [Cellular Automaton](./carbocat-ca.md) the CA component of CarboKitten
  - [Transport](./carbocat-transport.md) sediment transport
  - Composite model:
    - [CA + Production + Transport](./carbocat-cpt.md)
  - Under the hood:
    - [Stencil operations](./stencils.md)

!include README.md

# Scaffold

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
