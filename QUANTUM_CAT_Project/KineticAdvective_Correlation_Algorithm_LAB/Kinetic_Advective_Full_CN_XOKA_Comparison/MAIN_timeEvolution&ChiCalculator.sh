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
echo "Introduce an expression for the first derivative in y of the eigenstates as a function of x and j - NOT NECESSARY FOR CN2"
read trash
read diffy_functionEigenstates
echo "Introduce an expression for the second derivative in y of the eigenstates as a function of x and j - NOT NECESSARY FOR CN2"
read trash
read diffyy_functionEigenstates
echo " Do you want to use Xabier's correction? If yes input 1, else 0"
read trash
read oianguCorrection
echo " Introduce the maximum j of the chis you want me to calculate"
read trash
read jmax
echo "Introduce the tolerance for the total summed and integrated chi -> the algorithm will use as many chis as the necessary ones to achieve this sum -try 0.90"
read trash
read chiSumTolerance
echo " Introduce the number of trajectories to compute, note that it is more costly to obtain many in the XO algorithm"
read trash
read numTrajs
echo " Introduce the divisions in x and y for CN separated by a space: nx ny"
read trash trash
read nx1CN nx2CN
echo " Introduce the divisions in x and y for XO separated by a space: nx ny"
read trash trash
read nx1 nx2
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
echo " Introduce whether to use or not (1 or 0) the following Correlation Potentials"
echo " Kinetic Advective G J"
read trash trash trash trash
read Kin Adv G J
echo "Introduce wheteher you want to use custom initial trajectories or random not 1, or 0"
read trash
read customTrajs
echo "If 1 Introduce code in C++ that saves the numTraj initial positions into already initialized and allocated arrays initialPosx and initialPosy, else enter blank"
read trash
read customTrajsCode
echo "Input the xBound to consider the transmitance"
read trash
read xBound
echo "You want me to output only Data every n iterations, introduce n>1, else n=1"
read trash
read outputEvery
echo "If you wish to output the animation as an image write G, if you want it to be animated live write L"
read trash
read gif

if [[ $functionEigenstates != *"return"* ]]; then
        functionEigenstates="return ${functionEigenstates}"
fi

if [[ $diffy_functionEigenstates != *"return"* ]]; then
        diffy_functionEigenstates="return ${diffy_functionEigenstates}"
fi

if [[ $functionEigenstates != *"return"* ]]; then
        diffyy_functionEigenstates="return ${diffyy_functionEigenstates}"
fi

if [[ $psiIni != *"return"* ]]; then
    psiIni="return ${psiIni}"
fi
if [[ $potential != *"return"* ]]; then
    potential="return ${potential}"
fi

if [[ $gif == *"G"* ]]; then
    expLabel=""
    echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
    read trash
    read label
fi

if [[ $customTrajs == *"0"* ]]; then
customTrajsCode=" "
fi

XOKA_version="1"
if [[ $oianguCorrection == *"1"* ]]; then
XOKA_version="3"
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

    # XO KINETIC ADVECTIVEEE!!!!!-------------------------------------------------------------------------------------------------------------
    echo "XO KA time!!------------------------"
    expLabel="XO_KA_k0_${k0}_${label}"
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."

    potentialPlotFineness=0.007
    eigenstatesForSectionsIny="return 0;"
    diffxEigenstatesForSectionsIny="return 0;"
    diffxxEigenstatesForSectionsIny="return 0;"
    yjmax=0
    b_y=0


    ./EXE_codeFileGenerator_2D_XO_KINADV_BornHuang_tINDEP "$psiIni" "$potential" $mass1 $mass2 $nx1 $nx2 $x1min $x1max $x2min $x2max $dt $numIt $numTrajs $potentialPlotFineness $hbar $outputEvery "$functionEigenstates" "$diffy_functionEigenstates" "$diffyy_functionEigenstates" "$eigenstatesForSectionsIny" "$diffxEigenstatesForSectionsIny" "$diffxxEigenstatesForSectionsIny" $jmax $yjmax $b_y $chiSumTolerance $xBound $k0 $Kin $Adv $G $J $customTrajs "$customTrajsCode" $XOKA_version

    g++ -Wall -O CODE_simulator_XO_KinAdv.cpp -o EXE_simulator_XO_KinAdv
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_XO_KinAdv
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""

    echo "(2) GENERATING FINAL PLOT..."

    START_TIME=$SECONDS

        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating GIF..."
        gnuplot -e "set terminal gif size 1800, 900 font 'Helvetica,10' animate delay 1; set output '../ANIMATION_GIFS/XOKA_$expLabel-{$currentDateTime}_numTrajs($numTrajs)_XOKA.gif'; t=0; tmax=$tIts+1; while(t<tmax){ set key title 'KINETIC ADVECTIVE Algorithm :' at 1,25; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set pm3d map; set palette rgbformulae 33,13,10; set cbrange[0:100]; set colorbox; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' title 'Potential Energy heatmap and aproximate KA Bohmian trajectories', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' title 'Animated trajectories in colors - current in red, resultant in black, current traced in colors' w points lt 5 pt 4 ps 0.18 lc rgb 'black', for [i=0:($numTrajs-1)] 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' every ::(i*tmax)::(i*tmax+t) notitle w points lt 5 pt 4 ps 0.16 linecolor palette frac (i+1.0)/($numTrajs+1.0), for [i=0:($numTrajs-1)] 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' every ::(i*tmax+t)::(i*tmax+t) notitle w points lt 5 pt 7 ps 1.0 linecolor rgb 'red';unset pm3d; unset view; unset colorbox; unset key; t=t+1;}"

        echo "Done!"


    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""

    # CRANCK NICOLSON!!!!---------------------------------------------------------------------------------------------
    trajOption="R"
    echo "-------CRANCK NICOLSON time!-------- "
    echo "(2) Executing calculations for the time evolution...."
    echo " Generating code and compiling..."

    ./EXE_codeFileGenerator_2D_CN_tINDEP "$psiIni" "$potential" $nx1CN $nx2CN $x1min $x1max $x2min $x2max $dt $numIt $mass1 $mass2 $hbar  $outputEvery 1
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
    currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`

    gnuplot -e "set terminal gif size 1800,900 font 'Helvetica,11' animate delay 1; set output '../ANIMATION_GIFS/$expLabel-{$currentDateTime}__numTrajs($numTrajs)_CN2D.gif'; t=0; tmax=$tIts; array zeros[$tIts+1]; do for [i=1:($tIsts+1)] { zeros[i] = 0.0 }; while(t<tmax){ set multiplot;
    set origin 0.0, 0.0; clear; set size 0.5, 1.0; clear; set key title '|WF(x,y)|^2' at 1,25; set cbrange [0:0.01]; set xtics auto; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set pm3d map; set palette rgbformulae 33,13,10; set colorbox; splot 'DATA_probDensity_WF_CN.txt' index t title columnheader(1), for [i=0:$numTrajs] 'DATA_CN_Trajs_k=$k0.txt' using 1:2:(zeros[\$0+1]) index t every ::i notitle w points lt 5 pt 7 ps 1 lc rgb 'black'; unset pm3d; unset view; unset colorbox;
    set origin 0.5, 0.0; set size 0.5, 1.0; clear; set key title 'Potential Energy heatmap and Bohmian trajectories' at 1,25; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set pm3d map; set palette rgbformulae 33,13,10; set cbrange[0:100]; set colorbox; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' notitle, for [i=0:$numTrajs] 'DATA_CN_Trajs_k=$k0.txt' using 1:2:(zeros[\$0+1]) every ::i::(i+1):$tIts notitle w points lt 5 pt 4 ps 0.2 lc rgb 'black', for [i=0:$numTrajs] 'DATA_CN_Trajs_k=$k0.txt' using 1:2:(zeros[\$0+1]) index t every ::i notitle w points lt 5 pt 7 ps 1 lc rgb 'red';unset pm3d; unset view; unset colorbox; unset key; t=t+1; unset multiplot;}"



    # XO NO G,J------------------------------------------------------------------------------
    echo "XO No GJ time!!---------------------------------"
    expLabel="XO_NoGJ_k0_${k0}_${label}"
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."

    potentialPlotFineness=0.007
    eigenstatesForSectionsIny="return 0;"
    diffxEigenstatesForSectionsIny="return 0;"
    diffxxEigenstatesForSectionsIny="return 0;"
    yjmax=0
    b_y=0

    ./EXE_codeFileGenerator_2D_XO_KINADV_BornHuang_tINDEP "$psiIni" "$potential" $mass1 $mass2 $nx1 $nx2 $x1min $x1max $x2min $x2max $dt $numIt $numTrajs $potentialPlotFineness $hbar $outputEvery $xBound $k0 $customTrajs "$customTrajsCode" 2

    g++ -Wall -O CODE_simulator_XO_NoGJ.cpp -o EXE_simulator_XO_NoGJ
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_XO_NoGJ
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""

    echo "(3) PLOTTING"
    echo " Generating GIF..."
    gnuplot -e "set terminal gif size 1800, 900 font 'Helvetica,10' animate delay 1; set output '../ANIMATION_GIFS/XOKA_$expLabel-{$currentDateTime}_numTrajs($numTrajs)_VanillaXO.gif'; t=0; tmax=$tIts+1; while(t<tmax){ set key title 'KINETIC ADVECTIVE Algorithm :' at 1,25; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set pm3d map; set palette rgbformulae 33,13,10; set cbrange[0:100]; set colorbox; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' title 'Potential Energy heatmap and aproximate KA Bohmian trajectories', 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' title 'Animated trajectories in colors - current in red, resultant in black, current traced in colors' w points lt 5 pt 4 ps 0.18 lc rgb 'black', for [i=0:($numTrajs-1)] 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' every ::(i*tmax)::(i*tmax+t) notitle w points lt 5 pt 4 ps 0.16 linecolor palette frac (i+1.0)/($numTrajs+1.0), for [i=0:($numTrajs-1)] 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' every ::(i*tmax+t)::(i*tmax+t) notitle w points lt 5 pt 7 ps 1.0 linecolor rgb 'red';unset pm3d; unset view; unset colorbox; unset key; t=t+1;}"

    echo ""
    echo "All the simulations for k=" $k0 " completed.-------------------------"
    echo "Thanks for trusting XOA engine"
    echo ""

done

g++ -Wall -O3 CODE_proportionUnifier_ErrorCalculator.cpp -o EXE_proportionUnifier_ErrorCalculator
    ./EXE_proportionUnifier_ErrorCalculator $numk $tIts $ks

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XOKA_CNTrajs_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:i title sprintf('CN Trajects kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs CN_Trajs'; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:i title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XOKA_CNArea_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = $numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('CN Probability Density kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs CN_ProbabDens'; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:a[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XOKA_XONoGJ_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:t[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XONoGJ_CNTrajs_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 4*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:i title sprintf('CN Trajects kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference CN Trajs vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:a[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XONoGJ_CNArea_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = $numk+i }; array b[$numk+1]; do for [i=1:($numk)] { b[i] = 4*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('CN Prob Density kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference CN Prob Denst vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:b[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"
