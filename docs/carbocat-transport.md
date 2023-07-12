---
title: CarboKitten
subtitle: Sediment Transport
---

## Transport
The sediment that is produced is distributed into lower lying neighbour cells that are not occupied by a producer. A user defined fraction of sediment from a producer is transported, first divided equally to lower neighbours, cascading to its neighbours by splitting in half and so on. The cascade stops when the sediment reaches a minimal threshold.

Thus, this step has two free parameters: the transported fraction of produced carbonate and the lower threshold.

Apparent from the illustration in B13 Figure 4, a 8-cell neighbourhood is used. Nothing is mentioned about the order in which the transport is computed. We may tag transported sediment with a bit flip and assign a new lythofacies to transported sediment.

``` {.julia file=src/Burgess2013/Transport.jl}
module Transport

using MindTheGap.Stencil

function transport(::Type{B}, production::Matrix{Int}, n_species::Int) where {B <: Boundary{2}}

end

end
```

