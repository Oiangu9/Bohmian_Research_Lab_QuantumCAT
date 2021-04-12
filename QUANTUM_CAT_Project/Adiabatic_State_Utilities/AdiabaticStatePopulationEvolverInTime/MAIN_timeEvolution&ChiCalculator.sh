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
numIt=1
option=1
customTrajsCode=""

cd ./SOURCE_FILES

echo "Welcome to the XO algorithm vs CN2 software user interface!"
echo "(1) INPUT DATA"
echo " Introduce the ks youll use"
read trash
read ks
echo "Introduce the number of ks"
read trash
read numk
echo "Introduce the kind of intial WF to use GG for gaussian(x)*gaussian(y), GS for gaussian(x)*firstEigenstate(y)"
read trash
read psi
if [[ $psi == *"GS"* ]]; then
    echo "Introduce the following parameters for the intial WF g a1 a2 Lmax Lmin mux sigmax"
    read trash trash trash trash trash trash trash
    read g a1 a2 Lmax Lmin mux sigmax
else
    echo "Introduce the following parameters for the initial WF: sigmax sigmay mux muy"
    read trash trash trash trash
    read sigmax sigmay mux muy
fi
echo "Introduce the potential energy profile as a function of x and y"
read trash
read potential
echo " Introduce a general C++ syntax set of operations that allow the obtention of the value of the adiabatic eigenstates for the sections in x, as a function of y, x and j (dont forget the return)"
read trash
read functionEigenstates
echo " Would you like to use Xabier's correction? (1 if yes, 0 if not)"
read trash
read correctionOiangu
echo " Would you like to use Xabier's correction? (1 if yes, 0 if not)"
read trash
read jmax
echo " Introduce the number of trajectories to compute, note that it is more costly to obtain many in the XO algorithm"
read trash
read numTrajs
echo " Introduce the divisions in x and y for CN separated by a space: nx ny"
read trash trash
read nx1CN nx2CN
echo " Introduce the minimum and maximum of considered x in the grid as : xmin xmax"
read trash trash
read x1min x1max
echo " Introduce the minimum and maximum of considered y in the grid as : ymin ymax"
read trash trash
read x2min x2max
echo "Introduce the number of time iterations to consider"
read trash
read numIt
echo " Introduce the time step"
read trash
read dt
echo " Introduce the mass at x1 and mass at x2 as m1 m2"
read trash trash
read mass1 mass2
echo " Introduce the value of hbar"
read trash
read hbar
echo "Input the xBound to consider the transmitance"
read trash
read xBound
echo "You want me to output only Data every n iterations, introduce n>1, else n=1"
read trash
read outputEvery


if [[ $functionEigenstates != *"return"* ]]; then
        functionEigenstates="return ${functionEigenstates}"
fi


if [[ $psiIni != *"return"* ]]; then
    psiIni="return ${psiIni}"
fi
if [[ $potential != *"return"* ]]; then
    potential="return ${potential}"
fi



tIts=$(echo "$numIt / $outputEvery" | bc)


for k0 in $ks
do

if [[ $psi == *"GS"* ]]; then
    psiIni="double g=${g},a1=${a1},a2=${a2},Lmax=${Lmax},Lmin=${Lmin}, kx=${k0}, mux=${mux}, sigmax=${sigmax}, o, L; o= -(Lmax-Lmin)/2.0*(1.0/(1.0+exp((x-a1)/g))+1.0/(1.0+exp((-x+a2)/g)))-Lmin/2.0; L=-2*o; if(y>(o+L) || y<o){return 0.0;} else{return sqrt(2.0/L)*sin(PI*(y-o)/L)*pow(1.0/(2*PI*sigmax*sigmax),0.25)*exp(J*kx*x-0.25*pow((x-mux)/sigmax,2));}"
else
    psiIni="double sigmax=${sigmax}, sigmay=${sigmay}, mux=${mux}, muy=${muy}, kx=${k0}, ky=0; return pow(1.0/(2*PI*sigmax*sigmax),0.25)*exp(J*kx*x-0.25*pow((x-mux)/sigmax,2))* pow(1.0/(2*PI*sigmay*sigmay),0.25)*exp(J*ky*y-0.25*pow((y-muy)/sigmay,2));"
fi

    expLabel="CN_k0_${k0}_${label}"
    # CRANCK NICOLSON!!!!---------------------------------------------------------------------------------------------
    trajOption="R"
    echo "-------CRANCK NICOLSON time!-------- "
    echo "(2) Executing calculations for the time evolution...."
    echo " Generating code and compiling..."

    ./EXE_codeFileGenerator_2D_CN_tINDEP "$psiIni" "$potential" $nx1CN $nx2CN $x1min $x1max $x2min $x2max $dt $numIt $mass1 $mass2 $hbar  $outputEvery $version_of_XOKA 1
    g++ -Wall -O CODE_simulator_2D_CN_tINDEP.cpp -o EXE_simulator_2D_CN_tINDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_2D_CN_tINDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""

    echo "(3) COMPILING AND EXECUTING the ad-hoc CHI CALCULATOR"

    START_TIME=$SECONDS

    tIts=$(echo "$numIt / $outputEvery" | bc)
    dt_berri=$(echo "$dt * $outputEvery" | bc)

    path=$(pwd)/DATA_rawSimulationData_2D_CN.txt
    ./EXE_codeFileGenerator_CN2D_ChiCalculator "$path" "$functionEigenstates" $jmax $nx1CN $nx2CN $x1min $x1max $x2min $x2max $tIts $xBound $numTrajs $k0 $hbar $mass1 $mass2 $dt_berri 1

    g++ -Wall -O CODE_CN2D_ChiCalculator.cpp -o EXE_CN2D_ChiCalculator

    ./EXE_CN2D_ChiCalculator

    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""

    cp DATA_chiInfo_CN.txt ./DATA/DATA_chiInfo_CN_k0_${k0}.txt
    cp DATA_sumChiInfo_CN.txt ./DATA/DATA_sumChiInfo_CN_k0_${k0}.txt
    echo "(2) GENERATING FINAL PLOT..."

    START_TIME=$SECONDS

        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating GIF..."

        gnuplot -e "set terminal gif size 1950, 675 font 'Helvetica,15' animate delay 1; set output '../ANIMATION_GIFS/CN $currentDateTime jmax($jmax).gif'; t=0; tCN=0; countTrajs=0; tmax=2*$numIt*$numTrajs/$outputEvery; dx1=($x1max-($x1min))/$nx1; dx2=($x2max-($x2min))/$nx2; dx1CN=($x1max-($x1min))/$nx1CN; dx2CN=($x2max-($x2min))/$nx2CN; array posx1[$nx1+1]; array posx2[$nx2+1]; array posx1CN[$nx1CN+1]; array posx2CN[$nx2CN+1]; do for [i=1:($nx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($nx2+1)] { posx2[i] = $x2min + i*dx2 }; do for [i=1:($nx1CN+1)] { posx1CN[i] = $x1min + i*dx1CN }; do for [i=1:($nx2CN+1)] { posx2CN[i] = $x2min + i*dx2CN }; while(t<tmax){ set multiplot;
        set origin -0.005, 0.0; set size 0.5, 1.0; clear; set key title 'CN |WF(x,y)|^2' at 1,30; set xtics auto; set palette rgbformulae 33,13,10; set pm3d map; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set cbrange[0:0.012]; set colorbox; splot 'DATA_probDensity_WF_CN.txt' index tCN title columnheader(1); unset pm3d; unset view; unset colorbox; unset key;
        set origin 0.5, 0.0; set size 0.5, 1.0; clear; set key title '' inside; set xrange [$x1min:$x1max]; set yrange [0:0.4]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot  for [i=2:($jmax+2)] 'DATA_chiInfo_CN.txt' index tCN using (posx1CN[\$0+1]):i title sprintf('CN |Chi^{%d}(x)|',i-1) w l;
        t=t+2; tCN=tCN+1; if(tCN>($tIts)){tCN=0;}; unset multiplot;}"

        echo "Done!"


    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""


done


g++ -Wall -O3 CODE_proportionUnifier_ErrorCalculator.cpp -o EXE_proportionUnifier_ErrorCalculator
    ./EXE_proportionUnifier_ErrorCalculator $numk $tIts $ks

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XOKA_CNArea_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = $numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 1.0, 1.0; set yrange [0:0.5]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('CN Probability Density kx_0= k_{%d}',i) w l; unset multiplot;"

    
