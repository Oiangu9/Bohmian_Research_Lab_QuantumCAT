#!/bin/bash

pathFileWF=""
functionEigenstates=""
nx1=1
nx2=1
x1min=1.0
x1max=1.0
x2min=1.0
x2max=1.0
jmax=1

g++ -Wall -O CODE_GeneratorChiCalculator.cpp -o EXE_GeneratorChiCalculator

echo "Welcome to the Chi calculator user interface!"
echo "(1) INPUT DATA"
echo " Introduce the path of the file with the 2D wavefunction information:"
read trash
read pathFileWF
echo " Introduce a general C++ syntax set of operations that allow the obtention of the value of the adiabatic eigenstates for the sections in x, as a function of y, x and j (dont forget the return)"
read trash
read functionEigenstates
echo " Introduce the maximum j of the chis you want me to calculate"
read trash
read jmax
echo " Introduce the divisions in x and y separated by a space: nx1 nx2"
read trash trash
read nx1 nx2
echo " Introduce the minimum and maximum of considered x in the grid as : x1min x1max"
read trash trash
read x1min x1max
echo " Introduce the minimum and maximum of considered y in the grid as : x2min x2max"
read trash trash
read x2min x2max
echo "If you wish to output the animation as an image write G, if you want it to be animated live write L"
read trash    
read gif

if [[ $functionEigenstates != *"return"* ]]; then
        functionEigenstates="return ${functionEigenstates}"
fi

if [[ $gif == *"G"* ]]; then
    expLabel=""
    echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
    read expLabel
fi

echo "(2) COMPILING AND EXECUTING the ad-hoc CALCULATOR"

./EXE_GeneratorChiCalculator "$pathFileWF" "$functionEigenstates" $jmax $nx1 $nx2 $x1min $x1max $x2min $x2max

g++ -Wall -O CODE_ChiCalculator.cpp -o EXE_ChiCalculator

./EXE_ChiCalculator

echo "(3) GENERATING PLOT..."

if [[ $gif == *"G"* ]]; then
    currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
    echo " Generating GIF..." 

    gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './OUTPUT_IMAGES/$expLabel $currentDateTime jmax($jmax).gif'; set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|';plot for [i=2:($jmax/2+1)] 'DATA_chiInfo.txt' using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
    set origin 0, 0.5; set size 0.5, 0.5; set xrange [$x1min:$x1max]; set yrange [0:0.12]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|';plot for [i=($jmax/2+2):($jmax+1)] 'DATA_chiInfo.txt' using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
    set origin 0.49, 0; set size 0.5, 0.5; set xrange [-1:($jmax+1)]; set yrange [0:1.1]; set xticks 1; set xlabel 'jmax'; set ylabel 'integrate(\sum_{j=0}^{j=jmax}(|\chi^j(x)|^2))dx'; plot 'DATA_sumChiInfo.txt' using 1:2 with linespoint lw 3 pt 3 notitle, 'DATA_sumChiInfo.txt' using 1:2:(sprintf('%f', \$2)) with labels center offset -3.4,-.5 rotate by -45 notitle;
    set origin 0.5, 0.5; set size 0.5, 0.5; clear; set key title '|WF(x,y)|^2'; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; unset surface; splot 'DATA_plotWFInfo.txt' title columnheader(1); unset multiplot;"
    
    echo "Done!"

else

    gnuplot --persist -e "set terminal wxt size 1850,970; set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.5]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|';plot for [i=2:($jmax/2+1)] 'DATA_chiInfo.txt' using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
    set origin 0, 0.5; set size 0.5, 0.5; set xrange [$x1min:$x1max]; set yrange [0:0.12]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|';plot for [i=($jmax/2+2):($jmax+1)] 'DATA_chiInfo.txt' using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
    set origin 0.52, 0; set size 0.48, 0.5; set xtics 1; set xrange [-1:($jmax+1)]; set yrange [0:1.1]; set xlabel 'jmax'; set ylabel 'integrate(\sum_{j=0}^{j=jmax}(|\chi^j(x)|^2))dx'; plot 'DATA_sumChiInfo.txt' using 1:2 with linespoint lw 3 pt 3 notitle, 'DATA_sumChiInfo.txt' using 1:2:(sprintf('%f', \$2)) with labels center offset -3.4,-.5 rotate by -45 notitle; 
    set origin 0.5, 0.5; set size 0.5, 0.5; clear; set key title '|WF(x,y)|^2'; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; unset surface; splot 'DATA_plotWFInfo.txt' title columnheader(1); "
    
fi
        
