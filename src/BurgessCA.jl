# ~\~ language=Julia filename=src/BurgessCA.jl
# ~\~ begin <<docs/carbocat.md|src/BurgessCA.jl>>[init]
module BurgessCA

rules(x) = let s = [sum(x == i) for i in 1:3]
end

end
# ~\~ end
