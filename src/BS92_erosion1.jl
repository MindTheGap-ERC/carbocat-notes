# ~/~ begin <<docs/bosscher-1992.md#src/BS92.jl>>[init]
# couple RULSE erosion model

using DifferentialEquations
using CSV
using DataFrames
using Interpolations

# ~/~ begin <<docs/bosscher-1992.md#b92-model>>[init]
# ~/~ begin <<docs/bosscher-1992.md#carbonate-production>>[init]
g(gₘ, I₀, Iₖ, k, w) = gₘ * tanh(I₀/Iₖ * exp(-w * k))

#he(R,K,LS,C) = R * K * LS * C
# ~/~ end

pCO2::Float64 = 1
T::Float64 = 300
p::Float64 = 3.17 * 10^(-8) # m/s
K1::Float64  = exp(-356.3-0.0609*T+21834/T+126.8*log(T)-1684915/T^2)
K2::Float64  = exp(-107.8-0.03252*T+5151/T+38.9*log(T)-563713/T^2)
KC::Float64  = exp(-171.9-0.007799*T+2839.3/T+71.59*log(T))
KH::Float64  = exp(108.38+0.01985*T-6919/T-40.45*log(T)+669365/T^2)
A::Float64  = -0.49 + 8.074 * 0.0001 * (T - 273)
B::Float64 = -0.32 + 1.6 * 0.0001 * (T - 273)
gama_Ca::Float64 = exp(-4A * sqrt(0.1) * (1 + 5 * 10^-8 * B * sqrt(0.1)))
gama_alk::Float64 =exp(-A * sqrt(0.1) * (1 + 5.4 * 10^-8 * B * sqrt(0.1)))
Caeq = pCO2 .* (K1*KC*KH/(4*K2*gama_Ca*(gama_alk.^2))).^(1/3)
dis_const = 40 * 1000 * Caeq ./ 2710
dis_const2 = 0.007
e = 0.001 * dis_const2 * 1000

struct Parameters
     I₀::Float64
     Iₖ::Float64
     k::Float64
     gₘ::Float64
end

#R = 12 * 1.735 * 10 ^ -0.8188 # assuming 1m/y precipitation woth no seasonality

struct Erosion
    R::Float16 # rainfall
    K::Float16 # erodbility
    LS::Float16 # topography
    C::Float16 # vegetation?
end

g(p::Parameters, w) = g(p.gₘ, p.I₀, p.Iₖ, p.k, w)
#he(e::Erosion) = he(e.R,e.K,e.LS,e.C)
# ~/~ end
# ~/~ begin <<docs/bosscher-1992.md#b92-model>>[1]
function model(p::Parameters, s, h, t_end::Float64, h₀::Float64)
    
    if h - s(t) >= 0.0
        ∂h(h::Float64, _, t::Float64) = -g(p, h - s(t)) 
    elseif h - e >= h₀
        ∂h(h::Float64, _, t::Float64) = e
    else
        ∂h(h::Float64, _, t::Float64) = h - s(t)
    end

     #∂h(h::Float64, _, t::Float64) = let w = h - s(t) # 0.05*t is subsidence+ 0.0005 * t
      #    w >= 0.0 ? -g(p, h - s(t)) : e
     #end
     ode = ODEProblem(∂h, h₀, (0.0, t_end), Nothing)
     solve(ode, Euler(), dt=10.0, reltol=1e-6, saveat=1000.0)
end
# ~/~ end

function sealevel_curve()
     data = DataFrame(CSV.File("data/bs92-sealevel-curve.csv"))
     linear_interpolation(data.time, data.depth)
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

# ~/~ end