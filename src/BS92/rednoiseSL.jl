# rednoise sealevel
module 
using SignalAnalysis
using Plots
using CSV, Tables, DataFrames

function redNoiseSL()
    sl = rand(RedGaussian(80000))
    df = DataFrame(
        time = LinRange(0,80000,80000), 
        depth = sl)
    CSV.write("data/rednoise-SLcurve.csv",df)
    plot(df.time,df.depth)
end

redNoiseSL()