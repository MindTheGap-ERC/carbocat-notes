# ~/~ begin <<docs/bosscher-1992.md#src/figures/bs92-fig2.gnuplot>>[init]
set term svg
set xrange [-100:2000]
set yrange [50:0]
set ylabel "depth (m)"
set xlabel "intensity"

set arrow from 20, 45 to 20, 0 nohead
set arrow from 20, 45 to 2000, 45 nohead

set parametric
set trange [0:50]
plot 2000*exp(-0.1 * t), t w l, \
     1400*tanh(2000*exp(-0.1*t) / 400), t w l
# ~/~ end