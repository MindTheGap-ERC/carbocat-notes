# ~/~ begin <<docs/random-fields.md#src/RandomFields.jl>>[init]
module RandomFields

export Box, white_noise, fftfreq, mapk, gaussian_noise, power_law,
       discrete_scale_filter, sampled_gaussian_filter,
       fourier_power

using Random
using FFTW
using DataFrames
using Statistics

# ~/~ begin <<docs/random-fields.md#fourier>>[init]
wavenumber(n::Int, i::Int) = i - 1 > n ÷ 2 ? i - n - 1 : i - 1
fftfreq(n::Int, l::Float64, i::Int) = wavenumber(n, i) * 2π / l
fftfreq(n::Int, l::Float64) = (fftfreq(n, l, i) for i in 1:n)
# ~/~ end
# ~/~ begin <<docs/random-fields.md#box>>[init]
struct Box{dim}
    shape::Dims{dim}
    physical_size::NTuple{dim,Float64}
end

# ~/~ begin <<docs/random-fields.md#box-methods>>[init]
Box{dim}(n::Int, l::Float64) where {dim} = Box{dim}(
    Dims{dim}(Iterators.repeated(n, dim)),
    NTuple{dim}(Iterators.repeated(l, dim)))

Base.size(box::Box{dim}) where {dim} = *(box.shape...)
# ~/~ end
# ~/~ begin <<docs/random-fields.md#box-methods>>[1]
fftfreq(box::Box{dim}, idx::CartesianIndex{dim}) where {dim} =
    NTuple{dim,Float64}(fftfreq.(box.shape, box.physical_size, Tuple(idx)))
fftfreq(box::Box{dim}) where {dim} =
    Iterators.map(idx -> fftfreq(box, idx), CartesianIndices(box.shape))
# ~/~ end
# ~/~ end
# ~/~ begin <<docs/random-fields.md#filters>>[init]
struct Filter
    f::Function
end

(f::Filter)(k...) = f.f(k...)

Base.:*(f1::Filter, f2::Filter) =
    Filter((args...) -> f1.f(args...) * f2.f(args...))
# ~/~ end
# ~/~ begin <<docs/random-fields.md#filters>>[1]
sampled_gaussian_filter(sigma::Float64) = Filter(function (k...)
    absk2 = +((k .^ 2)...)
    exp(- sigma^2 * absk2 / 2)
end)
# ~/~ end
# ~/~ begin <<docs/random-fields.md#filters>>[2]
discrete_scale_filter(box::Box{2}, sigma::Float64) = Filter(function (k1, k2)
    γ = 1/3
    u = k1 * box.physical_size[1] / box.shape[1]
    v = k2 * box.physical_size[2] / box.shape[2]
    t = (sigma * box.shape[1] / box.physical_size[1])^2 / 2
    exp(t * ((1 - γ) * (cos(u) + cos(v)) + γ * (cos(u) * cos(v)) - 2 + γ))
end)
# ~/~ end
# ~/~ begin <<docs/random-fields.md#filters>>[3]
power_law(n::Float64) =
    Filter((k...) -> all(k .== 0.0) ? 1.0 : (+((k .^ 2)...))^(n/4))
# ~/~ end
# ~/~ begin <<docs/random-fields.md#noise>>[init]
white_noise(box::Box) = randn(box.shape...)

function mapk(f::Function, box::Box{dim}) where {dim}
    result = Array{Float64, dim}(undef, box.size...)
    for (i, k) in enumerate(fftfreq(box))
        result[i] = f(k)
    end
    result
end

function gaussian_noise(box::Box, P::Filter)
    noise = white_noise(box)
    f = fft(noise)
    for (i, k) in enumerate(fftfreq(box))
        f[i] *= P(k...)
    end
    real(ifft(f))
end
# ~/~ end

function fourier_power(box::Box{dim}, f::AbstractArray{Float64,dim}; n_bins = 64) where {dim}
    g = fft(f)
    power_table = DataFrame(
        ksqr = vec(collect(+((k .^ 2)...) for k in fftfreq(box))),
        power = vec(real(conj(g) .* g)) ./ size(box))
    sort!(power_table, [:ksqr])
    bin_size = size(box) ÷ n_bins
    power_table.bin = repeat(1:n_bins, inner=bin_size)
    agg = combine(
        groupby(power_table, :bin),
        :ksqr => mean, :power => mean)
    select(agg, :ksqr_mean, :power_mean)
end

end
# ~/~ end