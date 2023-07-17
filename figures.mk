# ~/~ begin <<docs/carbocat.md#figures.mk>>[init]
.RECIPEPREFIX = >
.PHONY: all _all

fig := docs/fig

all: _all

# ~/~ begin <<docs/carbocat.md#build>>[init]
targets += $(fig)/burgess2013-fig3.svg
targets += $(fig)/burgess2013-long-times.svg

docs/fig/burgess2013-fig3.svg: src/figures/ca.jl
> julia --project=. -e 'include("$<"); plot("$@")'

docs/fig/burgess2013-long-times.svg: src/figures/ca.jl
> julia --project=. -e 'include("$<"); plot_long_times("$@")'
# ~/~ end
# ~/~ begin <<docs/bosscher-1992.md#build>>[0]
targets += $(fig)/bs92-fig2.svg

$(fig)/bs92-fig2.svg: src/figures/bs92-fig2.gnuplot
> @mkdir -p $(@D)
> gnuplot $< > $@
# ~/~ end

_all: $(targets)
# ~/~ end