# plot_pngs.plt — writes PNG files
set term pngcairo size 1200,900 font ",12"
set datafile separator ","
set datafile columnheaders
set grid

# 1) Short orbit
set output "orbit_short.png"
set title "Mercury orbit (α=0, ~1 period)"
set xlabel "x [AU]"; set ylabel "y [AU]"
plot "orbit_short.csv" u 2:3 w l lw 2 t "trajectory"

# 2) Long stability
set output "orbit_long.png"
set title "Stable orbit (α=0, 100 periods)"
plot "orbit_long.csv" u 2:3 w l lw 1 t "trajectory"

# 3) Precession (α=0.01)
set output "orbit_precession.png"
set title "Precession of Mercury (α=0.01)"
plot "orbit_precession.csv" u 2:3 w l lw 2 t "trajectory"

# 4) ω(α) linear fit (requires alpha_omega.csv with >1 row)
set output "precession_fit.png"
set title "Angular velocity of precession vs α"
set xlabel "α [AU^2]"; set ylabel "ω [rad/year]"
f(x)=a*x
# If your gnuplot supports 'fit', this will work; otherwise comment the next two lines:
fit f(x) 'alpha_omega.csv' via a
plot 'alpha_omega.csv' u 1:2 w p pt 7 ps 1.5 t "data", \
     f(x) w l lw 2 lc "red" t sprintf("fit: ω=%.3e·α",a)

unset output
