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



echo "Welcome to the XO algorithm vs CN2 software user interface!"
echo "(1) INPUT DATA"
echo " Introduce the ks youll use"
read trash
read ks
echo "Introduce the number of ks"
read trash
read numk
echo " Introduce the initial 2D wavefunction as a function of x and y in C++ syntax- can use Eigen and cmath functions, code lines separated by ; in a same actual line"
read trash
read psiIni
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

    tIts=$(echo "$numIt / $outputEvery" | bc)

for k0 in $ks
do

    psiIni="double sigmax=2.0, sigmay=6.0, mux=-8.0, muy=0.0, kx=${k0}, ky=0; return pow(1.0/(2*PI*sigmax*sigmax),0.25)*exp(J*kx*x-0.25*pow((x-mux)/sigmax,2))* pow(1.0/(2*PI*sigmay*sigmay),0.25)*exp(J*ky*y-0.25*pow((y-muy)/sigmay,2));
    "

    expLabel="CN_k0_${k0}_${label}"
    # CRANCK NICOLSON!!!!---------------------------------------------------------------------------------------------
    trajOption="R"
    echo "-------CRANCK NICOLSON CHOSEN!-------- "
    echo "(2) Executing calculations for the time evolution...."
    echo " Generating code and compiling..."
    :'
    ./EXE_CN2D_codeFileGenerator_2D_CN "$psiIni" "$potential" $nx1CN $nx2CN $x1min $x1max $x2min $x2max $dt $numIt $mass1 $mass2 $hbar $option $outputEvery
    g++ -Wall -O CODE_CN2D_simulator_2D_CN_tINDEP.cpp -o EXE_CN2D_simulator_2D_CN_tINDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_CN2D_simulator_2D_CN_tINDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    '
    echo "(3) COMPILING AND EXECUTING the ad-hoc CHI CALCULATOR"
    
    START_TIME=$SECONDS
    
    tIts=$(echo "$numIt / $outputEvery" | bc)
    dt_berri=$(echo "$dt * $outputEvery" | bc)
    :'
    path=$(pwd)/DATA_rawSimulationData_2D_CN.txt
    ./EXE_CN2D_GeneratorChiCalculator "$path" "$functionEigenstates" $jmax $nx1CN $nx2CN $x1min $x1max $x2min $x2max $tIts $xBound $numTrajs $k0 $hbar $mass1 $mass2 $dt_berri

    #OJO! SI QUIERES QUE SE PLOTTEN EL POTENCIAL Y UNAS TRAYECTORIAS TAMBIEN TENDRAS QUE AÑADIR ESO EN EL CHI CALCULATOR ESTE TB! Y ASI PODER PLOTEAR UN CONTOUR PLOT DEL POTKIAL SOBRE EL QUE HACER LA GAUSSIANITA Y LAS TRAJS
        
    g++ -Wall -O CODE_CN2D_ChiCalculator.cpp -o EXE_CN2D_ChiCalculator

    ./EXE_CN2D_ChiCalculator
    
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo "(4) GENERATING PLOT..."

    START_TIME=$SECONDS
    '
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating GIF..." 

        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './OUTPUT_IMAGES/$expLabel CN2D $currentDateTime jmax($jmax).gif'; t=0; tmax=$tIts; while(t<tmax){ set multiplot; set origin 0.0, 0.0; set size 0.33, 0.5; clear; set key title; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=2:($jmax/2+2)] 'DATA_chiInfo.txt' index t using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
        set origin 0, 0.5; set size 0.33, 0.5; clear; set key title; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=($jmax/2+2):($jmax+2)] 'DATA_chiInfo.txt' index t using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
        set origin 0.33, 0; set size 0.33, 0.5; clear; set xrange [-1:($jmax+2)]; set yrange [0:1.1]; set xtics 1; set xlabel 'jmax'; set ylabel 'integrate(\sum_{j=0}^{j=jmax}(|\chi^j(x)|^2))dx'; plot 'DATA_sumChiInfo.txt' index t using 1:2 with linespoint lw 3 pt 3 notitle, 'DATA_sumChiInfo.txt' index t using 1:2:(sprintf('%f', \$2)) with labels center offset -3.4,-.5 rotate by -45 notitle;
        set origin 0.33, 0.5; set size 0.33, 0.5; clear; set key title '|WF(x,y)|^2'; set xtics auto; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; set colorbox; unset surface; splot 'DATA_plotWFInfo.txt' index t title columnheader(1); unset key; unset colorbox; set origin 0.66, 0.5; set size 0.33, 0.5; clear; set palette rgbformulae 33,13,10; set pm3d map; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2';splot 'DATA_plotWFInfo.txt' index t title columnheader(1); unset pm3d; unset view; unset colorbox; set origin 0.66, 0.0; set size 0.33, 0.5; clear; set key title 'Bohmian Trajectories'; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; plot 'DATA_CN_Trajs_k=$k0.txt' index t w points lt 5 pt 7 ps 1 lc rgb 'black'; t=t+1; unset multiplot;}"
        
        echo "Done!"

    else

        gnuplot -e "set terminal wxt size 1850,970; t=0; tmax=$numIt; while(t<tmax){ set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=2:($jmax/2+2)] 'DATA_chiInfo.txt' index t using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
        set origin 0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=($jmax/2+2):($jmax+2)] 'DATA_chiInfo.txt' index t using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
        set origin 0.49, 0; set size 0.5, 0.5; clear; set xrange [-1:($jmax+2)]; set yrange [0:1.1]; set xtics 1; set xlabel 'jmax'; set ylabel 'integrate(\sum_{j=0}^{j=jmax}(|\chi^j(x)|^2))dx'; plot 'DATA_sumChiInfo.txt' index t using 1:2 with linespoint lw 3 pt 3 notitle, 'DATA_sumChiInfo.txt' index t using 1:2:(sprintf('%f', \$2)) with labels center offset -3.4,-.5 rotate by -45 notitle;
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set key title '|WF(x,y)|^2'; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; unset surface; splot 'DATA_plotWFInfo.txt' index t title columnheader(1); t=t+$speed; unset multiplot;}"
        
    fi    
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    # XO KINETIC ADVECTIVEEE!!!!!-------------------------------------------------------------------------------------------------------------
    echo "XO KA time!!------------------------"
    expLabel="XO_KA_k0_${k0}_${label}"
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    
    potentialPlotFineness=0.013
    eigenstatesForSectionsIny="return 0;"
    diffxEigenstatesForSectionsIny="return 0;"
    diffxxEigenstatesForSectionsIny="return 0;"
    yjmax=0
    b_y=0
    
    ./EXE_XO-KA_codeFileGenerator_2D_XO_KINADV_BornHuang "$psiIni" "$potential" $mass1 $mass2 $nx1 $nx2 $x1min $x1max $x2min $x2max $dt $numIt $numTrajs $potentialPlotFineness $hbar $outputEvery "$functionEigenstates" "$diffy_functionEigenstates" "$diffyy_functionEigenstates" "$eigenstatesForSectionsIny" "$diffxEigenstatesForSectionsIny" "$diffxxEigenstatesForSectionsIny" $jmax $yjmax $b_y $chiSumTolerance $xBound $k0 1
    
    g++ -Wall -O CODE_XO-KA_simulator_2D_XO_KINADV_BornHuang_tINDEP.cpp -o EXE_XO-KA_simulator_2D_XO_KINADV_BornHuang_tINDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_XO-KA_simulator_2D_XO_KINADV_BornHuang_tINDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo "(3) GENERATING PLOT..."

    START_TIME=$SECONDS
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating GIF..." 

        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './OUTPUT_IMAGES/$expLabel XOKA $currentDateTime jmax($jmax).gif'; t=0; countTrajs=0; tmax=2*$numIt*$numTrajs/$outputEvery; dx1=($x1max-($x1min))/$nx1; dx2=($x2max-($x2min))/$nx2; array posx1[$nx1+1]; array posx2[$nx2+1]; do for [i=1:($nx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($nx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){ set multiplot; set origin 0.0, 0.0; set size 0.33, 0.5; clear; set key title; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=1:($jmax/2+2)] 'DATA_chiInfo.txt' index t/2 using (posx1[\$0+1]):i title sprintf('|Chi^{%d}(x)|',i-1) w l;
        set origin 0, 0.5; set size 0.33, 0.5; clear; set key title; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=($jmax/2+1):($jmax)] 'DATA_chiInfo.txt' index t/2 using (posx1[\$0+1]):i title sprintf('|Chi^{%d}(x)|',i-1) w l;
        set origin 0.33, 0; set size 0.33, 0.5; clear; set xrange [-1:($jmax+2)]; set yrange [0:1.1]; set xtics 1; set xlabel 'jmax'; set ylabel 'integrate(\sum_{j=0}^{j=jmax}(|\chi^j(x)|^2))dx'; plot 'DATA_sumChiInfo.txt' index t/2 using 1:2 with linespoint lw 3 pt 3 notitle, 'DATA_sumChiInfo.txt' index t/2 using 1:2:(sprintf('%f', \$2)) with labels center offset -3.4,-.5 rotate by -45 notitle;
        set origin 0.33, 0.5; set size 0.33, 0.5; clear; set xtics auto; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; if ((t-(2*$numIt*countTrajs)) >= 2*$numIt){ countTrajs=countTrajs+1; }; set pm3d map; set palette rgbformulae 33,13,10; set colorbox; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index 1 title 'Potential Energy Map', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' every ::1::(t/2+2) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red'; unset colorbox;
        set origin 0.66, 0; set size 0.33, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.04]; set key title 'x Conditional Wave Function Probability Density'; set xlabel 'Position x'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index t using (posx1[\$0+1]):1 title columnheader(1) w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' every ::(t/2)::(t/2+2) using 1:3 title 'x component of the Trajectory' w points lt 5 pt 6 ps 1 lc rgb 'red';
        set origin 0.66, 0.5; set size 0.33, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.04];set key title 'y Conditional Wave Function Probability Density';  set xlabel 'Position y'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index (t+1) using (posx2[\$0+1]):1 title columnheader(1) w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' every ::(t/2)::(t/2+2) using 2:3 title 'y component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        t=t+20; unset multiplot;}"
        #Ikusi ia index en el potential si realmente hace falta, a mi me da que nop
        #plotie nx ta nyren balixoak klar
        echo "Done!"

    else

        gnuplot -e "set terminal wxt size 1850,970; t=0; tmax=$numIt; while(t<tmax){ set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=2:($jmax/2+2)] 'DATA_chiInfo.txt' index t using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
        set origin 0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=($jmax/2+2):($jmax+2)] 'DATA_chiInfo.txt' index t using 1:i title sprintf('|Chi^{%d}(x)|',i-2) w l;
        set origin 0.49, 0; set size 0.5, 0.5; clear; set xrange [-1:($jmax+2)]; set yrange [0:1.1]; set xtics 1; set xlabel 'jmax'; set ylabel 'integrate(\sum_{j=0}^{j=jmax}(|\chi^j(x)|^2))dx'; plot 'DATA_sumChiInfo.txt' index t using 1:2 with linespoint lw 3 pt 3 notitle, 'DATA_sumChiInfo.txt' index t using 1:2:(sprintf('%f', \$2)) with labels center offset -3.4,-.5 rotate by -45 notitle;
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set key title '|WF(x,y)|^2'; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; unset surface; splot 'DATA_plotWFInfo.txt' index t title columnheader(1); t=t+$speed; unset multiplot;}"
        
    fi    
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
        
    
    
    # XO NO G,J------------------------------------------------------------------------------
    echo "XO No GJ time!!---------------------------------"
    expLabel="XO_NoGJ_k0_${k0}_${label}"
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    
    potentialPlotFineness=0.013
    eigenstatesForSectionsIny="return 0;"
    diffxEigenstatesForSectionsIny="return 0;"
    diffxxEigenstatesForSectionsIny="return 0;"
    yjmax=0
    b_y=0
    
    ./EXE_XO-KA_codeFileGenerator_2D_XO_KINADV_BornHuang "$psiIni" "$potential" $mass1 $mass2 $nx1 $nx2 $x1min $x1max $x2min $x2max $dt $numIt $numTrajs $potentialPlotFineness $hbar $outputEvery "$functionEigenstates" "$diffy_functionEigenstates" "$diffyy_functionEigenstates" "$eigenstatesForSectionsIny" "$diffxEigenstatesForSectionsIny" "$diffxxEigenstatesForSectionsIny" $jmax $yjmax $b_y $chiSumTolerance $xBound $k0 0
    
    g++ -Wall -O CODE_XO-KA_simulator_2D_XO_KINADV_BornHuang_tINDEP.cpp -o EXE_XO-KA_simulator_2D_XO_KINADV_BornHuang_tINDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_XO-KA_simulator_2D_XO_KINADV_BornHuang_tINDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo "(3) PLOTTING"
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  
                
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './OUTPUT_IMAGES/$expLabel XO_NoGJ_$currentDateTime x,yDivs($nx1,$nx2) x,ymin($x1min,$x2min) x,ymax($x1max,$x2max) nTrjs($numTrajs) dt($dt) tIt($numIt) mx,y($mass1,$mass2) outEvry($outputEvery).gif'; t=0; countTrajs=0; tmax=2*$numIt*$numTrajs/$outputEvery; dx1=($x1max-($x1min))/$nx1; dx2=($x2max-($x2min))/$nx2; array posx1[$nx1+1]; array posx2[$nx2+1]; do for [i=1:($nx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($nx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){  set multiplot; 
        set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.04]; set key title 'x Conditional Wave Function Probability Density'; set xlabel 'Position x'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index t using (posx1[\$0+1]):1 title columnheader(1) w l, 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' every ::(t/2)::(t/2+2) using 1:3 title 'x component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.04];set key title 'y Conditional Wave Function Probability Density';  set xlabel 'Position y'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index (t+1) using (posx2[\$0+1]):1 title columnheader(1) w l, 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' every ::(t/2)::(t/2+2) using 2:3 title 'y component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red'; set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; if ((t-(2*$numIt*countTrajs)) >= 2*$numIt){ countTrajs=countTrajs+1; }; set pm3d map; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index 1 title 'Potential Energy Map', 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' every ::1::(t/2+2) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
        t=t+20; unset multiplot; }"
        
    else
        
        
        gnuplot -e "set terminal WTX size 1800,900; t=0; countTrajs=0; tmax=2*$tIt*$numTrajs; dx1=($xmax-($xmin))/$xDivs; dx2=($ymax-($ymin))/$yDivs; array posx1[$xDivs+1]; array posx2[$yDivs+1]; do for [i=1:($xDivs+1)] { posx1[i] = $xmin + i*dx1 }; do for [i=1:($yDivs+1)] { posx2[i] = $ymin + i*dx2 }; while(t<tmax){  set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$xmin:$xmax]; set yrange [0:0.1]; set xlabel 'Position x'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index t using (posx1[\$0+1]):1 title 'x Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' every ::(t/2)::(t/2+2) using 1:3 title 'x component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$ymin:$ymax]; set yrange [0:0.1]; set xlabel 'Position y'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'y Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' every ::(t/2)::(t/2+2) using 2:3 title 'y component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$xmin:$xmax]; set yrange [$ymin:$ymax]; set xlabel 'Position x'; set ylabel 'Position y'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set pm3d map; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' every ::1::(t/2+2) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
        t=t+2; unset multiplot; }"

	fi

    echo ""
    echo "All the simulations for k=" $k0 " completed.-------------------------"
    echo "Thanks for trusting XOA engine"
    echo ""
    
done

g++ -Wall -O3 CODE_proportionUnifier_ErrorCalculator.cpp -o EXE_proportionUnifier_ErrorCalculator
    ./EXE_proportionUnifier_ErrorCalculator $numk $tIts $ks
        
    gnuplot -e "set terminal png size 1800,500; set output '$label Transmitace_XOKA_CNTrajs_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:i title sprintf('CN Trajects kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs CN_Trajs'; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:i title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"
    
    gnuplot -e "set terminal png size 1800,500; set output '$label Transmitace_XOKA_CNArea_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = $numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('CN Probability Density kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs CN_ProbabDens'; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:a[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '$label Transmitace_XOKA_XONoGJ_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:t[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '$label Transmitace_XONoGJ_CNTrajs_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 4*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:i title sprintf('CN Trajects kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference CN Trajs vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:a[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"
    
    gnuplot -e "set terminal png size 1800,500; set output '$label Transmitace_XONoGJ_CNArea_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = $numk+i }; array b[$numk+1]; do for [i=1:($numk)] { b[i] = 4*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('CN Prob Density kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference CN Prob Denst vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:b[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"
    
