# =============================
# gnuplot script for projectile
# =============================

set datafile separator ","

# ---- Task 2: Trajectories no drag ----
set title "Trajectories without drag (different dt)"
set xlabel "x [m]"
set ylabel "y [m]"
plot \
    "traj_n10.csv"   using 2:3 with lines title "n=10", \
    "traj_n20.csv"   using 2:3 with lines title "n=20", \
    "traj_n50.csv"   using 2:3 with lines title "n=50", \
    "traj_n100.csv"  using 2:3 with lines title "n=100", \
    "traj_n200.csv"  using 2:3 with lines title "n=200", \
    "traj_n500.csv"  using 2:3 with lines title "n=500"

pause -1 "Press Enter to continue"

# ---- Task 2: Global error vs dt ----
set title "Global error vs dt (log-log)"
set xlabel "dt [s]"
set ylabel "Global error [m]"
set logscale xy
plot "errors.csv" using 1:2 with linespoints title "Eglob"

pause -1 "Press Enter to continue"

# ---- Task 3: Drag force trajectories ----
set title "Trajectories with drag force (a=0)"
set xlabel "x [m]"
set ylabel "y [m]"
plot \
    "drag_traj_D0.csv"     using 2:3 with lines title "D=0", \
    "drag_traj_D100.csv"   using 2:3 with lines title "D=1e-4", \
    "drag_traj_D200.csv"   using 2:3 with lines title "D=2e-4", \
    "drag_traj_D500.csv"   using 2:3 with lines title "D=5e-4", \
    "drag_traj_D1000.csv"  using 2:3 with lines title "D=1e-3"

pause -1 "Press Enter to continue"

# ---- Task 3: Range vs firing angle ----
set title "Range vs firing angle"
set xlabel "Angle [deg]"
set ylabel "Range [m]"

plot \
    "range_D0.csv"   using 1:2 with linespoints title "D=0", \
    "range_D1.csv"   using 1:2 with linespoints title "D=1e-3", \
    "range_D2.csv"   using 1:2 with linespoints title "D=2e-3"

pause -1 "Press Enter to continue"

# ---- Task 4: Altitude correction ----
set title "Trajectories with/without altitude correction"
set xlabel "x [km]"
set ylabel "y [km]"
set xrange [0:*]
set yrange [0:*]
plot \
    "altcorr_angle35_a0.csv"   using ($2/1000):($3/1000) with lines title "35deg, a=0", \
    "altcorr_angle35_aPos.csv" using ($2/1000):($3/1000) with lines title "35deg, a>0", \
    "altcorr_angle45_a0.csv"   using ($2/1000):($3/1000) with lines title "45deg, a=0", \
    "altcorr_angle45_aPos.csv" using ($2/1000):($3/1000) with lines title "45deg, a>0"

pause -1 "Press Enter to finish"
