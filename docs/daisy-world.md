---
title: Daisy World
subtitle: reading notes on Wood's 2008 review
---

# Daisy World
Reading the review paper on Daisy World by Wood et al. 2006.

We have two types of daisy: one black, the other white. Now we'll model the population growth of these two species. For a single species we have,

[$$\partial_t \alpha = \alpha(p - \alpha)\beta(x) - \alpha \gamma(x),$$]{#eq:daisy-single}

describing logistic growth with birthrate $\beta$ and deathrate $\gamma$. $\alpha$ is the fractional coverage by a species, $p$ the proportion of habitable ground. We may model the groth of this species by varying $\beta$ and $\gamma$.

Considering the meaning of the $\beta$ and $\gamma$ parameters, it seems reasonable to have $\beta > 0$ and $0 \le \gamma \le 1$. That means $\beta$ is probably best sampled from a Rayleigh distribution, which can be achieved by taking two independent Gaussian variables as a 2d vector and compute the length of that vector. For the $\gamma$ value, I propose taking the phase component of those same variables and take either the sine or cosine of that.

``` {.julia file=src/Daisy/SingleSpecies.jl}
using MindTheGap.RandomFields

function field_parameters(b::Box{2}, sigma, beta_mu)
end
```

Now for Daisy world, land is occupied by either one or the other species or bare, so $\alpha_a = p - \alpha_b$, and we can write two coupled equations:

[$$\begin{align}
\partial_t \alpha_a &= \alpha_a((p - \alpha_a - \alpha_b) \beta(x_a) - \gamma)\\
\partial_t \alpha_b &= \alpha_b((p - \alpha_a - \alpha_b) \beta(x_b) - \gamma)
\end{align}$$]{#eq:daisy-coupled}

This assumes that both species have the same death rate, but may prefer different locations.

``` {.julia file=src/Daisy.jl}
print("Hello Daisy")
```
