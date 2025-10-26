# plots.plt — Gnuplot script for RK4 oscillator (Tasks 2a–4)
# -> Adjusted for Task 3 filenames with 6 decimals from std::to_string
# Usage: gnuplot plots.plt

set datafile separator ","
set term pngcairo size 1200,800 font ",12"
set grid
set key outside top center horizontal

# ---------------------- Task 2a ----------------------
set output "energy_2a.png"
set title "Task 2a: Energies vs t (alpha=0, F0=0)"
set xlabel "t"; set ylabel "Energy"
plot "trajectory_2a.csv" u 1:4 w l lw 2 t "E_kin", \
     "" u 1:5 w l lw 2 t "E_pot", \
     "" u 1:6 w l lw 2 t "E_tot"

set output "timeseries_2a.png"
set title "Task 2a: x(t) and v(t)"
set xlabel "t"; set ylabel "x(t), v(t)"
plot "trajectory_2a.csv" u 1:2 w l lw 2 t "x(t)", \
     "" u 1:3 w l lw 2 t "v(t)"

set output "phase_2a.png"
set title "Task 2a: Phase space v(x)"
set xlabel "x"; set ylabel "v"
plot "trajectory_2a.csv" u 2:3 w l lw 2 t "v(x)"

# ---------------------- Task 2b ----------------------
set output "timeseries_compare_2b.png"
set title "Task 2b: x(t) numerical vs exact (alpha=0.1)"
set xlabel "t"; set ylabel "x(t)"
plot "compare_2b.csv" u 1:2 w l lw 2 t "x_num", \
     "" u 1:3 w l lw 2 dt (5,5) t "x_exact"

set output "energy_2b.png"
set title "Task 2b: Energies vs t (alpha=0.1)"
set xlabel "t"; set ylabel "Energy"
plot "trajectory_2b.csv" u 1:4 w l lw 2 t "E_kin", \
     "" u 1:5 w l lw 2 t "E_pot", \
     "" u 1:6 w l lw 2 t "E_tot"

set output "phase_2b.png"
set title "Task 2b: Phase space v(x) (alpha=0.1)"
set xlabel "x"; set ylabel "v"
plot "trajectory_2b.csv" u 2:3 w l lw 2 t "v(x)"

# ---------------------- Task 3 ----------------------
# NOTE: filenames use 6 decimals due to std::to_string in your C++.
#       Make sure you've run:
#       ./rk4_osc task=3 F0=0 Omega=0 tmax=50 N=10000 alphas3=0.0001,0.1,0.5,1.95
set output "dissipation_task3.png"
set title "Task 3: x(t) for various alpha (Fext=0)"
set xlabel "t"; set ylabel "x(t)"
plot \
  "trajectory_3_alpha_0.000100.csv" u 1:2 w l lw 2 t "alpha=1e-4", \
  "trajectory_3_alpha_0.100000.csv" u 1:2 w l lw 2 t "alpha=0.1", \
  "trajectory_3_alpha_0.500000.csv" u 1:2 w l lw 2 t "alpha=0.5", \
  "trajectory_3_alpha_1.950000.csv" u 1:2 w l lw 2 t "alpha=1.95"

# ---------------------- Task 4 ----------------------
# resonance_all.csv columns: 1:Omega, 2:alpha, 3:amp
alpha_list = "0.01 0.1 0.5 1.0"

set output "resonance_task4.png"
set title "Task 4: Resonance amplitude vs Ω (F0=1)"
set xlabel "Ω_ext"; set ylabel "Amplitude x_max"
set xrange [0.1:2.0]
set logscale y
set format y "10^{%L}"

# Analytic overlay for one alpha (choose any; assignment allows one)
alphaA = 0.1
m=1.0; F0=1.0; w0=1.0
tau(a) = 2.0*m/a
Aanal(O,a) = F0/(sqrt((w0*w0 - O*O)**2 + (2.0*O/tau(a))**2))

plot for [a in alpha_list] "resonance_all.csv" u (abs($2 - real(a))<1e-9 ? $1 : 1/0):3 w l lw 2 t sprintf("alpha=%s", a), \
     Aanal(x,alphaA) w l lw 2 dt (5,5) t sprintf("analytic (alpha=%.2f)", alphaA)

# ---------------------- Done ----------------------
set output
print "Plots written: energy_2a.png, timeseries_2a.png, phase_2a.png, timeseries_compare_2b.png, energy_2b.png, phase_2b.png, dissipation_task3.png, resonance_task4.png"
