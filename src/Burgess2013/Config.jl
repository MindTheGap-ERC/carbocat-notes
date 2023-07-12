# ~/~ begin <<docs/carbocat.md#src/Burgess2013/Config.jl>>[init]
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
# ~/~ end