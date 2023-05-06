# CarboCAT
CarboCAT is primarily based on a very simple cellular automaton (CA). We may explore this CA as a first step in implementing the model in Julia.

``` {.julia file=src/BurgessCA.jl}
module BurgessCA

rules(x) = let s = [sum(x == i) for i in 1:3]
end

end
```
