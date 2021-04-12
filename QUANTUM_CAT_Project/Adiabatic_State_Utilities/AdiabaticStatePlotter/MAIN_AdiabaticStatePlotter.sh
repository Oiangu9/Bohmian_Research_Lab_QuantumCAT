#!/bin/bash

echo "Introduce the adiabatic states for sections in x as a function of y, x and j"
read trash
read functionAdiabatic
echo "Input the potential field as a funciton of x and y"
read trash
read potential
echo "Introduce the maximum j for the states "
read trash
read jmax
echo "Introduce nx ny"
read trash trash
read nx ny
echo "Introduce the axis bounds for x: xmin xmax"
read trash trash
read xmin xmax
echo "Now for y dimension: ymin ymax"
read tras trash
read ymin ymax
echo "Introduce a name for the output animation file"
read trash
read expLabel

echo "EXECUTING CALCULATIONS....."

g++ -Wall -O CODE_GeneratorAdiabaticStatePlotter.cpp -o EXE_GeneratorAdiabaticStatePlotter

./EXE_GeneratorAdiabaticStatePlotter "$functionAdiabatic" "$potential" $jmax $nx $ny $xmin $xmax $ymin $ymax

g++ -Wall -O CODE_AdiabaticStatePlotter.cpp -o EXE_AdiabaticStatePlotter

./EXE_AdiabaticStatePlotter


echo "EXPORTING GNUPLOT ANIMATION..."

currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`

gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './OUTPUT_IMAGES/$expLabel $currentDateTime jmax($jmax).gif'; ix=0; ix_max=$nx; while(ix<ix_max){ xNow=$xmin+ix*($xmax-$xmin)/$nx; set multiplot; set origin 0.0, 0.0; set size 0.25, 1.0; clear; set key title; set xrange[-0.5:0.5]; set yrange [$ymin:$ymax]; set xlabel sprintf('Re{psi^{j}_{%f}(y)}', xNow); set ylabel 'Position y'; plot for [i=2:($jmax/2+1)] 'DATA_adiabaticStatePlot.txt' index ix using i:1 title sprintf('Re{psi^{%d}_{%f}(y)}',i-2, xNow) w l; set origin 0.25, 0.0; set size 0.25, 1.0; clear; set key title; set xrange[-0.5:0.5]; set yrange [$ymin:$ymax]; set xlabel sprintf('Re{psi^{j}_{%f}(y)}', xNow); set ylabel 'Position y'; plot for [i=($jmax/2+2):($jmax+2)] 'DATA_adiabaticStatePlot.txt' index ix using i:1 title sprintf('Re{psi^{%d}_{%f}(y)}',i-2, xNow) w l; set origin 0.5, 0.0; set size 0.5, 0.5; clear; set xrange [$xmin:$xmax]; set yrange [$ymin:$ymax]; set xlabel 'Position x'; set ylabel 'Position y'; set pm3d map; set palette rgbformulae 33,13,10; set colorbox; splot 'DATA_potentialToPlot.txt' title 'Potential Energy Map'; unset colorbox; set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$xmin:$xmax]; set yrange [$ymin:$ymax]; set xlabel 'Position x'; set ylabel 'Position y'; set pm3d map; set palette rgbformulae 33,13,10; splot exp(-2*(x-xNow)**2) title sprintf('Now plotting section in x=%f', xNow); unset colorbox; ix=ix+1; unset multiplot;}"

