#!/usr/bin/gnuplot

# set terminal epslatex size 8cm,5cm
# set terminal pngcairo size 1024,768
set terminal pngcairo size 1920,1080


set linestyle 1 lt rgb "steelblue" lw 1
set linestyle 2 lt rgb "red" lw 1
set linestyle 3 lt rgb "dark-spring-green" lw 1
set linestyle 4 lt rgb "steelblue" dashtype 4 lw 1
set linestyle 5 lt rgb "tan1" dashtype 4 lw 1

set autoscale

set samples 10000

set logscale y

set xtics auto
set ytics auto

set grid mytics ytics
set grid mxtics xtics

set output "FePt.png"

set title "Diffraction pattern of FePt at 300 K"
set xlabel "2θ (degrees)"
set ylabel "Intensity (arbitrary units)"

set xrange [20:120]
# set xrange [20:55]
# set xrange [22:26]
# set xrange [47:50]
# set yrange [0.01:10]

# plot "xrd_pattern.csv" using (2*$1):2 title "Simulated" with line ls 1, "AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta.txt" using 1:($2/10015.308134319914) title "Measurement" with line ls 2

plot "AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta_signal.txt" using 1:2 notitle with line ls 2

# plot "AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta.txt" using 1:2 title "Measurement" with line ls 1, "MgO_Background_XRD_20_120_2-Theta_Omega_processed.txt" using 1:2 title "Background" with line ls 2
# plot "MgO_Background_XRD_20_120_2-Theta_Omega.txt" using 1:(6.684593284*$2) notitle with line ls 1, "AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta.txt" using 1:2 notitle with line ls 2
# plot "MgO_Background_XRD_20_120_2-Theta_Omega.txt" using 1:(31.20317942417027*$2) notitle with line ls 1, "MgO_Background_XRD_20_120_2-Theta_Omega_processed.txt" using 1:2 title "Blurred" with line ls 2 


# set output "FePt_split.png"
# 
# set title "Diffraction pattern of tetragonal FePt (a=2.722, c=3.709) at 300 K"
# set xlabel "2θ (degrees)"
# set ylabel "Intensity (arbitrary units)"
# 
# # set xrange [20:120]
# set xrange [20:55]
# # set xrange [22:26]
# # set xrange [47:50]
# set yrange [1e-2:10]
# 
# set arrow from 23.9663966397,graph(0,0) to 23.9663966397,graph(1,1) nohead
# set label "001" at 23.9663966397,graph(1,1.015) center
# set arrow from 32.8622862286,graph(0,0) to 32.8622862286,graph(1,1) nohead
# set label "100" at 32.8622862286,graph(1,1.015) center
# set arrow from 40.5430543054,graph(0,0) to 40.5430543054,graph(1,1) nohead
# set label "110" at 40.5430543054,graph(1,1.015) center
# set arrow from 47.6537653765,graph(0,0) to 47.6537653765,graph(1,1) nohead
# set label "111" at 47.6537653765,graph(1,1.015) center
# set arrow from 49.0789078908,graph(0,0) to 49.0789078908,graph(1,1) nohead
# set label "002" at 49.0789078908,graph(1,1.015) center
# 
# plot "xrd_full_pattern.csv" using (2*$1):2 title "Simulated" with line ls 1, "AJA_1249_MgO-FePt-Pt_190s_XRD_Phil_Theta_2-Theta.txt" using 1:($2/10149) title "Measurement" with line ls 2
