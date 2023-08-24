# ~/~ begin <<docs/bosscher-1992.md#src/BS92/fig8.jl>>[init]
#using MindTheGap.BS92
using Plots

function main()
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

          plot!(t,-result[:,10],
          legend_position=:none,
          title="strati-ele with time",
          xaxis=("time (years)"),
          yaxis=("stratigraphy height (m)", :right),
          inset=(1, bbox(0.61, 0.60, 0.45, 0.28)),
          subplot=3,
          framestyle=:box)

     #plot(result,t)
     #savefig("docs/fig/bs92-fig8.html")
end

main()
# ~/~ end
