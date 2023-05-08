# ~\~ language=Gnuplot filename=src/figures/plot-tanh.gnuplot
# ~\~ begin <<docs/carbocat.md|src/figures/plot-tanh.gnuplot>>[init]
set term svg
set xrange [-1:10]
plot tanh(exp(-x))
# ~\~ end
