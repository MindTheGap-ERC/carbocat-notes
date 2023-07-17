# ~/~ begin <<docs/bosscher-1992.md#src/BS92.jl>>[init]
module BS92

using DifferentialEquations
using CSV
using DataFrames
using Interpolations

function sealevel_curve()
     data = DataFrame(CSV.File("data/bs92-sealevel-curve.csv"))
     #scale(interpolate(data.depth, BSpline(Cubic())), LinRange(0.0, 80_000.0, 862))
     linear_interpolation(data.time, data.depth)
end

# ~/~ begin <<docs/bosscher-1992.md#carbonate-production>>[init]
g(gₘ, I₀, Iₖ, k, w) = if w > 0.0
     gₘ * tanh(I₀/Iₖ * exp(-w * k))
else
     0.0
end
# ~/~ end

struct Parameters
     I₀::Float64
     Iₖ::Float64
     k::Float64
     gₘ::Float64
end

g(p::Parameters, w) = g(p.gₘ, p.I₀, p.Iₖ, p.k, w)

function model(p::Parameters, s, t_end::Float64, h₀::Float64)
     ∂h(h::Float64, _, t::Float64) = -g(p, h - s(t))
     ode = ODEProblem(∂h, h₀, (0.0, t_end), Nothing)
     #solve(ode, AutoTsit5(Rosenbrock23()), reltol=1e-8, saveat=1000.0)
     solve(ode, Euler(), dt=10.0, reltol=1e-6, saveat=1000.0)
end

struct Scenario
     param::Parameters
     sealevel
     t_end::Float64
end

model(s::Scenario, h₀::Float64) = model(s.param, s.sealevel, s.t_end, h₀)

SCENARIO_A = Scenario(
     Parameters(2000.0, 250.0, 0.05, 0.005),
     sealevel_curve(),
     80_000.0)

SCENARIO_B = Scenario(
     Parameters(2000.0, 250.0, 0.05, 0.005),
     let sc = sealevel_curve()
          x -> sc(x) * 0.9
     end,
     80_000.0)

end
# ~/~ end