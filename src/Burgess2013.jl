# ~/~ begin <<docs/carbocat.md#src/Burgess2013.jl>>[init]
module Burgess2013

module Types
    # ~/~ begin <<docs/carbocat-transport.md#ck-types>>[init]
    struct Deposit{N}
        amount::NTuple{N, Float64}
    end

    Base.zero(::Type{Deposit{N}}) where {N} =
        Deposit{N}(ntuple(_ -> zero(Float64), N))
    # ~/~ end
    # ~/~ begin <<docs/carbocat.md#ck-types>>[0]
    export Product

    struct Product
        species::Int
        amount::Float64
    end

    Base.zero(::Type{Product}) = Product(0, 0.0)
    # ~/~ end
end

include("Burgess2013/Config.jl")
include("Burgess2013/CA.jl")
include("Burgess2013/Production.jl")
include("Burgess2013/Transport.jl")

end
# ~/~ end