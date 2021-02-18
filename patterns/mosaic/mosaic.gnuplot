#!/usr/bin/gnuplot

# set terminal epslatex size 8cm,5cm
set terminal png


set linestyle 1 lt rgb "steelblue" lw 1
set linestyle 2 lt rgb "red" lw 1

set autoscale

set samples 10000

# set logscale y

set xrange [0:180]

set xtics auto
set ytics auto

set grid mytics ytics
set grid mxtics xtics

set output "mosaic_KBr.png"

set title "(1 1 1) diffraction pattern of KBr at 300 K"
set xlabel "2θ (degrees)"
set ylabel "Intensity (arbitrary units)"

plot "xrd_pattern.csv" using (2*$1):2 notitle with line ls 1
