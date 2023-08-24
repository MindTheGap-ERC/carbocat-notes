using DifferentialEquations
using CSV
using DataFrames
using Interpolations

function g(gₘ, I₀, Iₖ, k, w) 
    g = gₘ * tanh(I₀/Iₖ * exp(-w * k))
end

struct Parameters
    I₀::Float64
    Iₖ::Float64
    k::Float64
    gₘ::Float64
end

g(p::Parameters, w) = g(p.gₘ, p.I₀, p.Iₖ, p.k, w)


function sealevel_curve()
    data = DataFrame(CSV.File("data/bs92-sealevel-curve.csv"))
    linear_interpolation(data.time, data.depth)
end

function erosion()
    dis_const2 = 0.007
    e1 = 0.001 * dis_const2 * 1000 * 0
    #e2 = 
end


function model(p::Parameters,s,t_end::Float64, h₀::Float64)
    e = erosion()
    function ∂h(h, _, t)
        let w = h - s(t) + 0.005 * t
            if w > 0.0 
                -g(p, w) 
            else
                min(e,h₀ - h)
            end
        end
    end
    ode = ODEProblem(∂h, h₀, (0.0, t_end), Nothing)
    h = solve(ode, Euler(), dt=10.0, reltol=1e-6, saveat=1000.0)
end

struct Scenario
    param::Parameters
    sealevel
    t_end::Float64
end

model(s::Scenario, h₀::Float64) = model(s.param, s.sealevel, s.t_end, h₀::Float64)

SCENARIO_A = Scenario(
    Parameters(2000.0, 250.0, 0.05, 0.005),
    sealevel_curve(),
    80_000.0,
    )
 


using Plots

     h0 = LinRange(0, 200, 101)
     result = hcat([model(SCENARIO_A, h).u for h in h0]...) #BS92.BS92.
     t = LinRange(0, 80000, 81)

     plotlyjs()

     plot(h0, result',
          xaxis=("initial depth (m)"),
          yaxis=("depth (m)", :flip),
          legend_position=:none, lc=:steelblue,
          size=(700, 700), fontfamily="Merriweather,serif")

     plot!(t, SCENARIO_A.sealevel(t),#BS92.
          title="sea level curve", titlelocation=:left,
          titlefontsize=12,
          xaxis=("time (years)"),
          yaxis=("depth (m)", :flip),
          guidefontsize=10,
          legend_position=:none,
          lc=:steelblue,
          inset=(1, bbox(0.11, 0.60, 0.45, 0.28)),
          subplot=2,
          framestyle=:box)

          plot!(t,result[:,10],
          legend_position=:none,
          title="strati-ele with time",
          xaxis=("time (years)"),
          yaxis=("stratigraphy height (m)", :flip),
          inset=(1, bbox(0.61, 0.60, 0.45, 0.28)),
          subplot=3,
          framestyle=:box)

     #plot(result,t)
     #savefig("docs/fig/bs92-fig8.html")


#main()    