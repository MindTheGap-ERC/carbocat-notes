# ~\~ language=Julia filename=src/Stencil.jl
# ~\~ begin <<docs/stencils.md|src/Stencil.jl>>[init]
module Stencil

export Boundary, Reflected, Periodic, Constant, stencil, convolution

# ~\~ begin <<docs/stencils.md|boundary-trait>>[init]
abstract type Boundary{dim} end
struct Reflected{dim}       <: Boundary{dim} end
struct Periodic{dim}        <: Boundary{dim} end
struct Constant{dim, value} <: Boundary{dim} end
# ~\~ end
# ~\~ begin <<docs/stencils.md|offset-indexing>>[init]
function offset_value(::Type{Periodic{dim}}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex) where {dim}
    z[mod1.(Tuple(i + Δi), size(z))...]
end

function offset_value(::Type{Reflected{dim}}, z::AbstractArray{T, dim}, i::CartesianIndex, Δi::CartesianIndex) where {T, dim}
    clip(i, a, b) = (i < a ? a + a - i : (i > b ? b + b - i : i))
    z[clip.(Tuple(i + Δi), ones(Int, dim), size(z))...]
end

function offset_value(::Type{Constant{dim, value}}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex) where {dim, value}
    j = i + Δi
    (checkbounds(Bool, z, j) ? z[j] : value)
end
# ~\~ end
# ~\~ begin <<docs/stencils.md|stencil-operation>>[init]
function stencil(::Type{T}, ::Type{BT}, n::NTuple{dim,Int}, f::Function) where {T, dim, BT <: Boundary{dim}}
    m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = zeros(T, n)

    function(z_in::AbstractArray{T, dim}, z_out::AbstractArray{T, dim}, args...)
        @assert (size(z_in) == size(z_out)) "sizes of arrays need to be equal"
        shape = size(z_in)
        for i in CartesianIndices(shape)
            for (k, Δi) in enumerate(CartesianIndices(stencil_shape))
                stencil[k] = offset_value(BT, z_in, i, Δi)
            end
            z_out[i] = f(stencil, args...)
        end
    end
end

function convolution(::Type{B}, kernel::Array{T, dim}) where { dim, T, B <: Boundary{dim} }
    stencil(T, B, size(kernel), s -> sum(s .* kernel))
end
# ~\~ end

end
# ~\~ end
