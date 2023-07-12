# Initial conditions

``` {.julia file=src/RandomFields.jl}
module RandomFields

export Box, white_noise, fftfreq, mapk, gaussian_noise, power_law,
       discrete_scale_filter, sampled_gaussian_filter,
       fourier_power

using Random
using FFTW
using DataFrames
using Statistics

<<fourier>>
<<box>>
<<filters>>
<<noise>>

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
```

## The box

``` {.julia #box}
struct Box{dim}
    shape::Dims{dim}
    physical_size::NTuple{dim,Float64}
end

<<box-methods>>
```

``` {.julia #box-methods}
Box{dim}(n::Int, l::Float64) where {dim} = Box{dim}(
    Dims{dim}(Iterators.repeated(n, dim)),
    NTuple{dim}(Iterators.repeated(l, dim)))

Base.size(box::Box{dim}) where {dim} = *(box.shape...)
```

## Fourier basics
The Discrete Fourier Transform is defined as

$$F(\vec{y})_l = \sum_{j = 1}^N y_j \exp (i k_l x_j).$$

Our convention here is to include the usual factor $2\pi$ into the values of $k_l$, as this simplifies our bookkeeping.

### Frequencies
Discrete Fourier space is usually stored such that we start with positive ascending frequencies and the other half is negative ascending frequencies. For a grid with $N$ points and physical size $L$, we can give the grid coordinates

$$\vec{x} = {L \over N} \big[0, 1, 2, \dots, N-1\big].$$

The corresponding Fourier frequencies are,

$$\vec{k} = {2\pi \over L} \big[0, 1, 2, \dots, N_q-1, N_q, 1-N_q, 2-N_q, \dots, -2, -1\big],$$

where $N_q = N / 2$ is the Nyquist wave number, and $k_{N_q}$ is the Nyquist frequency (including the factor $2\pi$).

``` {.julia #fourier}
wavenumber(n::Int, i::Int) = i - 1 > n ÷ 2 ? i - n - 1 : i - 1
fftfreq(n::Int, l::Float64, i::Int) = wavenumber(n, i) * 2π / l
fftfreq(n::Int, l::Float64) = (fftfreq(n, l, i) for i in 1:n)
```

This easily extends to higher dimensions.

``` {.julia #box-methods}
fftfreq(box::Box{dim}, idx::CartesianIndex{dim}) where {dim} =
    NTuple{dim,Float64}(fftfreq.(box.shape, box.physical_size, Tuple(idx)))
fftfreq(box::Box{dim}) where {dim} =
    Iterators.map(idx -> fftfreq(box, idx), CartesianIndices(box.shape))
```

### Filters
We will make ample use of the Fourier convolution theorem. We'll define a set of filters in terms of their Fourier transforms.

``` {.julia #filters}
struct Filter
    f::Function
end

(f::Filter)(k...) = f.f(k...)

Base.:*(f1::Filter, f2::Filter) =
    Filter((args...) -> f1.f(args...) * f2.f(args...))
```

The Fourier transform of a Gaussian is a Gaussian. A Gaussian distribution with standard deviation $\sigma$,

$$G(x) = {1 \over {\sigma\sqrt{2 \pi}}} \exp\left(-{x^2 \over 2\sigma^2}\right),$$

has the Fourier transform,

$$\hat{G}(k) = \exp\left(- {{\sigma^2 k^2} \over 2}\right).$$

We can approximate a Gaussian convolution by sampling the Gaussian on a grid

``` {.julia #filters}
sampled_gaussian_filter(sigma::Float64) = Filter(function (k...)
    absk2 = +((k .^ 2)...)
    exp(- sigma^2 * absk2 / 2)
end)
```

For repeated filters with smaller width, it is better to use a *discrete scale filter*.

``` {.julia #filters}
discrete_scale_filter(box::Box{2}, sigma::Float64) = Filter(function (k1, k2)
    γ = 1/3
    u = k1 * box.physical_size[1] / box.shape[1]
    v = k2 * box.physical_size[2] / box.shape[2]
    t = (sigma * box.shape[1] / box.physical_size[1])^2 / 2
    exp(t * ((1 - γ) * (cos(u) + cos(v)) + γ * (cos(u) * cos(v)) - 2 + γ))
end)
```

## Gaussian Random Fields
We consider random functions $f : \mathbb{R}^n \to \mathbb{R}$. In general it is very hard to say anything meaningful about random functions without additional information. The first question would be: given any point $x \in \mathbb{R}^n$, what is the probability distribution for $y = f(x)$? This is called the *one-point* distribution function $P(y\ |\ y=f(x))$. If a random field is *homogeneous* or *translation invariant*, the one-point distribution does not depend on location, 

$$P(y\ |\ y=f(x_1)) = P(y\ |\ y=f(x_2))\ \forall\ x_1, x_2 \in \mathbb{R}^n, y \in \mathbb{R}.$$

If we want to know more about the random field, we can look at the $two-point$ distribution $P(y_1, y_2\ |\ y_1 = f(x_1), y_2 = f(x_2))$. For a homogeneous field, the two-point function should only depend on the distance $r = x_2 - x_1$, so we can write $P(y_1, y_2\ |\ y_1 = f(x), y_2 = f(x+r))$, or shorter $P(f(x_1), f(x_2)). Now we also have a tool to talk about isotropy: a random field is isotropic if all its $n$-point distributions depend only on the distance $|r_i|$.

We may consider the probability density function for an entire realisation of the field. Given the (possibly infinite) set of points $\vec{x}$, and all its values $\vec{y} = f.(\vec{x})$, we write the probability density $P(\vec{y}\ |\ \vec{y} = f.(\vec{x}))$.

A beautiful result that I will not derive here is that for Gaussian processes, all we need to know is the two-point distribution function, which fully determined by the two-point correlation function

$$\xi(x, x+r):=\langle f(x)^{\star}f(x+r) \rangle$$

When we talk about homogenous fields (as we usually do), $\xi$ does not depend on $x$, and by the *ergodic theorem* we can replace the ensemble average by an integral over space (and dividing by the volume).

$$\xi(r) = {1 \over V} \int f(x)^{\star} f(x+r) {\rm d} x.$$

This only works if $\xi(r)$ vanishes when $r$ approaches inifinity. The two-point probability density function is then givin by

$$P(\vec{y}) = {1 \over \sqrt{(2\pi)^n \det M}} \exp\left(
    -{1 \over 2} \vec{y}^{\dagger} M^{-1} \vec{y}\right),$$

where $M$ is the covariance matrix $M_{ij} = \langle y_i^{\star} y_j \rangle$.

### Power spectrum
Great things happen when looking at Gaussian random fields in Fourier space. We evaluate the two-point correlation function in Fourier space at the origin $x = 0$.

$$\begin{align}
\xi(r) &= \left\langle
        \int \hat{f}(k)^{\star} {{\rm d}k \over {(2\pi)^n}}
        \int \hat{f}(k') e^{-ik'r} {{\rm d}k \over {(2\pi)^n}}
        \right\rangle\\
       &= \int \int \left\langle \hat{f}(k)^{\star} \hat{f}(k') \right\rangle
       {{\rm d}k \over {(2\pi)^n}} e^{-ik'r} {{\rm d}k' \over {(2\pi)^n}}.
\end{align}$$

A translation in real space gives a phase shift in Fourier space, so

$$\mathcal{F}[f(x + r)] = \hat{f}(k) e^{-ikr}.$$

Since our integral needs to be invariant under translation, we may say,

$$\left\langle \hat{f}(k)^{\star}\hat{f}(k')\right\rangle
= \left\langle \hat{f}(k)^{\star}\hat{f}(k')\right\rangle e^{i(k-k')r},$$

for any $r \in \mathbb{R}^n$. Which is only true for $k = k'$. This means that *for any random field that is translation invariant, the Fourier coefficients are uncorrelated!* We can define the function $P(k)$ such that,

$$\left\langle \hat{f}(k)^{\star}\hat{f}(k') \right\rangle := (2 \pi)^n P(k) \delta(k - k'),$$

and,

$$\xi(r) = \int P(k) e^{ikr} {{\rm d}k \over {(2\pi)^n}}.$$

This gives us a way of imposing a chosen correlation function on some white noise.
Let's compute the power of white noise, for which $\xi(r) = \delta(r) \sigma^2$, then $P(k) = \sigma^2$. Let's see that it is so.

``` {.julia file=src/scripts/Power.jl}
using CosmicStructureFormation.InitialConditions: Box, fourier_power
using DataFrames
```



``` {.julia #filters}
power_law(n::Float64) =
    Filter((k...) -> all(k .== 0.0) ? 1.0 : (+((k .^ 2)...))^(n/4))
```

## Make some noise

``` {.julia #noise}
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
```
