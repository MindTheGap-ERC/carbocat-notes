# ~/~ begin <<docs/carbocat.md#src/Burgess2013/Production.jl>>[init]
module Production

using ..Config: Species, Iₖ, k, gₘ

# ~/~ begin <<docs/bosscher-1992.md#carbonate-production>>[init]
g(gₘ, I₀, Iₖ, k, w) = gₘ * tanh(I₀/Iₖ * exp(-w * k))

struct Parameters
     I₀::Float64
     Iₖ::Float64
     k::Float64
     gₘ::Float64
end

g(p::Parameters, w) = g(p.gₘ, p.I₀, p.Iₖ, p.k, w)
# ~/~ end
# ~/~ begin <<docs/bosscher-1992.md#carbonate-production>>[1]
function model(p::Parameters, s, t_end::Float64, h₀::Float64)
     ∂h(h::Float64, _, t::Float64) = let w = h - s(t)
          w >= 0.0 ? -g(p, h - s(t)) : 0.0
     end
     ode = ODEProblem(∂h, h₀, (0.0, t_end), Nothing)
     solve(ode, Euler(), dt=10.0, reltol=1e-6, saveat=1000.0)
end
# ~/~ end

function production_rate(I₀::Float64, s::Species, w::Float64)
    g(gₘ(s), I₀, Iₖ(s), k(s), w)
end

function production_rate(I₀::Float64, specs::Vector{Species}, spec_map::Matrix{Int}, w::Matrix{Int})
    production_rate.(I₀, Iterators.map(i -> specs[i], spec_map), w)
end

end
# ~/~ end