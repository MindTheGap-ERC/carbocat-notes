---
title: CarboKitten
subtitle: CarboCAT in Julia
---

# CarboCAT 2013
CarboCAT is primarily based on a very simple cellular automaton (CA). We may explore this CA as a first step in implementing the model in Julia.

``` {.julia file=src/Burgess2013.jl}
module Burgess2013

include("Burgess2013/Config.jl")
include("Burgess2013/CA.jl")
include("Burgess2013/Production.jl")
include("Burgess2013/Transport.jl")

end
```

``` {.julia file=src/Burgess2013/Config.jl}
module Config

struct Species
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::Float64
    extinction_coefficient::Float64
    saturation_intensity::Float64
end

Iₖ(s::Species) = s.saturation_intensity
k(s::Species) = s.extinction_coefficient
gₘ(s::Species) = s.maximum_growth_rate

end
```

## Cellular Automaton
The paper talks about cycling the order of preference for occupying an empty cell at each iteration. This means that the rules change slighly every iteration.

``` {.julia #cycle-permutation}
cycle_permutation(n_species::Int) =
    (circshift(1:n_species, x) for x in Iterators.countfrom(0))
```

The `stencil` function has an `args...` variadic arguments that are forwarded to the given rule. This means we can create a `rules` function that we pass the preference order as a second argument.

``` {.julia #burgess2013-rules}
function rules(neighbourhood::Matrix{Int}, order::Vector{Int})
    cell_species = neighbourhood[3, 3]
    neighbour_count(species) = sum(neighbourhood .== species)
    if cell_species == 0
        for species in order
            n = neighbour_count(species)
            if 6 <= n && n <= 10
                return species
            end
        end
        0
    else
        n = neighbour_count(cell_species) - 1
        (4 <= n && n <= 10 ? cell_species : 0)
    end
end
```

This function is not yet adaptible to the given rule set. Such a modification is not so hard to make. 

The paper talks about a 50x50 grid initialized with uniform random values.

``` {.julia file=src/Burgess2013/CA.jl}
module CA

using MindTheGap.Stencil

<<cycle-permutation>>
<<burgess2013-rules>>

function run(::Type{B}, init::Matrix{Int}, n_species::Int) where {B <: Boundary{2}}
    Channel{Matrix{Int}}() do ch
        target = Matrix{Int}(undef, size(init))
        put!(ch, init)
        stencil_op = stencil(Int, B, (5, 5), rules)
        for perm in cycle_permutation(n_species)
            stencil_op(init, target, perm)
            init, target = target, init
            put!(ch, init)
        end
    end
end

end
```

First, let us reproduce Figure 3 in Burgess 2013.

![First 8 generations](fig/burgess2013-fig3.svg)

By eye comparison seems to indicate that this CA is working the same. I'm curious to the behaviour after more iterations. Let's try 10, 100, 1000 and so on.

![Assymptotic behaviour](fig/burgess2013-long-times.svg)

The little qualitative change between 100 and 1000 iterations would indicate that this CA remains "interesting" for a long time.

On my laptop I can run about 150 iterations per second with current code. When using periodic boundaries, I get to 1500 iterations per second, which is peculiar. A lot can still be optimized.

<details><summary>Plotting code</summary>

``` {.julia file=src/figures/ca.jl}
using MindTheGap.Burgess2013.CA
using MindTheGap.Stencil: Reflected
using MindTheGap.Utility
using GnuplotLite

function plot_array(w::Int, h::Int, msgs; xtics="set xtics", ytics="set ytics")
    ch = Channel{String}() do ch
        g = Gnuplot(ch)

        top_margin = 0.002
        bottom_margin = 0.1
        left_margin = 0.05
        right_margin = 0.05
        inner_margin = 0.02
        vert_inner_margin = 0.025

        # h * ph + (h-1) * inner_margin + top_margin + bottom_margin = 1
        plot_height = (1.0 - (h-1)*inner_margin - top_margin - bottom_margin) / h
        plot_width = (1.0 - (w-1)*vert_inner_margin - left_margin - right_margin) / w

        g |> send("set multiplot; unset xtics")
        for (i, msg) in enumerate(msgs)
            grid_x = (i - 1) % w
            grid_y = (i - 1) ÷ w
            if grid_x == 0
                tmargin = 1.0 - top_margin - (plot_height + vert_inner_margin) * grid_y
                bmargin = tmargin - plot_height
                g |> send("set tmargin at screen $(tmargin); set bmargin at screen $(bmargin)")
                g |> send(ytics)
            end
            if grid_x == 1
                g |> send("unset ytics")
            end
            if grid_y == h-1 && grid_x == 0
                g |> send(xtics)
            end
            lmargin = left_margin + (plot_width + inner_margin) * grid_x
            rmargin = lmargin + plot_width
            g |> send("set lmargin at screen $(lmargin); set rmargin at screen $(rmargin)")
            g |> msg
        end
        g |> send("unset multiplot")
    end

    send(join(ch, "\n"))
end

function plot(output::String)
    init = rand(0:3, 50, 50)
    result = Iterators.take(CA.run(Reflected{2}, init, 3), 8)

    gnuplot() do g
        g |>
            send("set term svg size 900, 440") |>
            send("set output '$(output)'") |>
            send("load 'data/blue-to-red.pal'") |>
            send("set size square") |>
            send("set xrange [0:50]; set yrange [0:50]") |>
            send("unset colorbox") |>
            plot_array(4, 2, (send("data" => r) *
                              send("plot \$data matrix u (\$1+0.5):(\$2+0.5):3 t'' w image pixels")
                              for r in result);
                       xtics="set xtics 10, 10, 50",
                       ytics="set ytics 10, 10, 50")
    end
end

function plot_long_times(output::String)
    init = rand(0:3, 50, 50)
    result = select(CA.run(Reflected{2}, init, 3), [10, 100, 1000])

    gnuplot() do g
        g |>
            send("set term svg size 900, 240") |>
            send("set output '$(output)'") |>
            send("load 'data/blue-to-red.pal'") |>
            send("set size square") |>
            send("set xrange [0:50]; set yrange [0:50]") |>
            send("unset colorbox") |>
            plot_array(3, 1, (send("data" => r) *
                              send("plot \$data matrix u (\$1+0.5):(\$2+0.5):3 t'' w image pixels")
                              for r in result);
                       xtics="set xtics 10, 10, 50",
                       ytics="set ytics 10, 10, 50")
    end
end
```

``` {.make file=figures.mk}
.RECIPEPREFIX = >
.PHONY: all _all

fig := docs/fig

all: _all

<<build>>

_all: $(targets)
```

``` {.make #build}
targets += $(fig)/burgess2013-fig3.svg
targets += $(fig)/burgess2013-long-times.svg

docs/fig/burgess2013-fig3.svg: src/figures/ca.jl
> julia --project=. -e 'include("$<"); plot("$@")'

docs/fig/burgess2013-long-times.svg: src/figures/ca.jl
> julia --project=. -e 'include("$<"); plot_long_times("$@")'
```

</details>

``` {.julia file=src/Burgess2013/Production.jl}
module Production

using ..Config: Species, Iₖ, k, gₘ

<<carbonate-production>>

function production_rate(I₀::Float64, s::Species, w::Float64)
    g(gₘ(s), I₀, Iₖ(s), k(s), w)
end

function production_rate(I₀::Float64, specs::Vector{Species}, spec_map::Matrix{Int}, w::Matrix{Int})
    production_rate.(I₀, Iterators.map(i -> specs[i], spec_map), w)
end

end
```

### Crowding
In crowded areas carbonate production rates are reduced. For cells where

$$n_{min} \le n \le n_{opt}$$

and

$$n_{opt} \le n \le n_{max}$$

($n_{min}$ and $n_{max}$ for living cells are 4 and 10)

we have a linear increase and linear decrease of production rate (i.e. a triangle function).

## Subsidence
Subsidence refers to gradual lowering or lifting of the underlying floor bed. This could be either sea level rise due to climate change, but also tectonic activity. Sea level also changes according to a given recipe with say three sinusoidals (e.g. Milankovich cycles). When a cell gets "subaerial exposure", i.e. the water level drops below the cell elevation (stupid jargon), accumulation stops and the cell becomes dormant. On reflooding, the cell resumes activity. From the text it is not entirely clear, but I think deactivated cells don't take part in the CA, so they count as dead neighbours.

- [reference on accommodation](http://strata.uga.edu/sequence/accommodation.html), also links to a model called [SedFlux](https://github.com/mcflugen/sedflux).

$$T + E = S + W$$

Saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.


## Steps

1. Update waterdepth, given subsidence
2. Update sea level elevation, given eustatics
3. Run CA
4. Compute thickness of carbonate production
5. Compute sediment transport

We need to keep a state with the following components:

- height map
- species
- global time, implying:
  - sea level and subsidence

Every cycle we may export a layer of sediment to archive.
