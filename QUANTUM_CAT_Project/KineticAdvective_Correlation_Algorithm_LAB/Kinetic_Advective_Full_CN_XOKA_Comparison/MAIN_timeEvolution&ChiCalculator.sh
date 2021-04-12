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
echo " Would you like to use Xabier's correction? (1 if yes, 0 if not)"
read trash
read correctionOiangu
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

tIts=$(echo "$numIt / $outputEvery" | bc)

version_of_XOKA="1"
if [[ $correctionOiangu == *"1"* ]]; then
version_of_XOKA="3"
fi

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


    ./EXE_codeFileGenerator_2D_XO_KINADV_BornHuang_tINDEP "$psiIni" "$potential" $mass1 $mass2 $nx1 $nx2 $x1min $x1max $x2min $x2max $dt $numIt $numTrajs $potentialPlotFineness $hbar $outputEvery "$functionEigenstates" "$diffy_functionEigenstates" "$diffyy_functionEigenstates" "$eigenstatesForSectionsIny" "$diffxEigenstatesForSectionsIny" "$diffxxEigenstatesForSectionsIny" $jmax $yjmax $b_y $chiSumTolerance $xBound $k0 $Kin $Adv $G $J $customTrajs "$customTrajsCode" $version_of_XOKA

    g++ -Wall -O CODE_simulator_XO_KinAdv.cpp -o EXE_simulator_XO_KinAdv
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_XO_KinAdv
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""

    # Combined information calculation!
    echo "Combine XO trajectories with CN!!---------------"
    echo " (1) CODE COMPILATION and EXECUTION"

    DATA_traj_XO_File="DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt"
    DATA_wf_CN_File="DATA_rawSimulationData_2D_CN.txt"

    ./EXE_codeFileGenerator_CN_Conditional_WFs_with_XO_Trajectories "$DATA_traj_XO_File" "$DATA_wf_CN_File" $nx1CN $nx2CN $jmax $x1min $x2min $x1max $x2max $tIts $numTrajs "$functionEigenstates" $chiSumTolerance $mass1 $mass2 $hbar 1

    g++ -Wall -O CODE_Simulator_CN_Conditional_WFs_with_XO_Trajectories.cpp -o EXE_Simulator_CN_Conditional_WFs_with_XO_Trajectories

    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_Simulator_CN_Conditional_WFs_with_XO_Trajectories

    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""

    echo "(2) GENERATING FINAL PLOT..."

    START_TIME=$SECONDS

        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating GIF..."

        gnuplot -e "set terminal gif size 1950, 1350 font 'Helvetica,10' animate delay 1; set output '../ANIMATION_GIFS/$expLabel CN vs XOKA $currentDateTime jmax($jmax) numTrajs($numTrajs).gif'; t=0; tCN=0; countTrajs=0; tmax=2*$numIt*$numTrajs/$outputEvery; dx1=($x1max-($x1min))/$nx1; dx2=($x2max-($x2min))/$nx2; dx1CN=($x1max-($x1min))/$nx1CN; dx2CN=($x2max-($x2min))/$nx2CN; array posx1[$nx1+1]; array posx2[$nx2+1]; array posx1CN[$nx1CN+1]; array posx2CN[$nx2CN+1]; do for [i=1:($nx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($nx2+1)] { posx2[i] = $x2min + i*dx2 }; do for [i=1:($nx1CN+1)] { posx1CN[i] = $x1min + i*dx1CN }; do for [i=1:($nx2CN+1)] { posx2CN[i] = $x2min + i*dx2CN }; while(t<tmax){ set multiplot;
        set origin -0.005, 0.66; set size 0.25, 0.33; clear; set key title 'CN |WF(x,y)|^2' at 1,30; set xtics auto; set palette rgbformulae 33,13,10; set pm3d map; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set cbrange[0:0.012]; set colorbox user origin 0.215,0.72 size 0.008,0.22; splot 'DATA_probDensity_WF_CN.txt' index tCN title columnheader(1); unset pm3d; unset view; unset colorbox; unset key;
        set origin 0, 0.33; set size 0.25, 0.33; clear; set key title 'Potential Energy heatmap and KA trajectories'; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set pm3d map; set palette rgbformulae 33,13,10; set cbrange[0:100]; set colorbox; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' title 'Potential Energy Map', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' every ::1::(t/2+2) title 'traced KA Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';unset pm3d; unset view; unset colorbox; unset key;
        set origin -0.007,0.0; set size 0.25, 0.33; clear; set key title 'R^4 G_x(x,y,t) for this trajectory KA y(t)'; set palette rgbformulae 33,13,10; set pm3d map; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set cbrange[-2e-7<*:*<1e-6]; set colorbox; splot 'DATA_Full_WF_G_J_CN.txt' index t/2 using 1:2:3 title 'Full WF R^4 G_x(x,y,t)'; unset pm3d; unset view; unset colorbox; unset key;
        set origin 0.25, 0.66; set size 0.25, 0.33; clear; set xrange [-1:($jmax+2)]; set yrange [0:1.7]; set xtics 1; set xlabel 'jmax'; set ylabel 'integrate(\sum_{j=0}^{j=jmax}(|\chi^j(x)|^2))dx'; set key title ' ' inside; plot 'DATA_sumChiInfo_XO.txt' index t/2 using 1:2 with linespoint lw 3 pt 3 title 'KA', 'DATA_sumChiInfo_XO.txt' index t/2 using 1:2:(sprintf('%f', \$2)) with labels center offset -5.4,-.5 rotate by -45 notitle, 'DATA_sumChiInfo_CN.txt' index tCN using 1:2 with linespoint lw 3 pt 3 title 'exact CN', 'DATA_sumChiInfo_CN.txt' index tCN using 1:2:(sprintf('%f', \$2)) with labels center offset -3.4,-.5 rotate by -45 notitle, 'DATA_sumChiInfo_CN_with_KA_traj.txt' index t/2 using 1:2 with linespoint lw 3 pt 3 title 'Cn w KA trajs', 'DATA_sumChiInfo_CN_with_KA_traj.txt' index t/2 using 1:2:(sprintf('%f', \$2)) with labels center offset -1.4,-.5 rotate by -45 notitle; unset key;
        set origin 0.25, 0.33; set size 0.25, 0.33; clear; set xtics auto; set key title ' ' inside; set xrange [$x1min:$x1max]; set yrange [-0.7:0.7]; set xlabel 'Position x'; set ylabel 'Kin'; set y2tics nomirror; set y2label 'Adv'; set y2range [-0.7:0.7]; plot 'DATA_KinAdv_x.txt' index t/2 using (posx1[\$0+1]):1 title 'Re{Kin(x)}' axes x1y1 w l, 'DATA_KinAdv_x.txt' index t/2 using (posx1[\$0+1]):2 title 'Im{Kin(x)}' axes x1y1 w l, 'DATA_KinAdv_x.txt' index t/2 using (posx1[\$0+1]):3 title 'Re{Adv(x)}' axes x1y2 w l, 'DATA_KinAdv_x.txt' index t/2 using (posx1[\$0+1]):4 title 'Im{Adv(x)}' axes x1y2 w l; unset y2tics; unset y2label;
        set origin 0.25,0.0; set size 0.25, 0.33; clear; set xtics auto; set key title 'R^4 J_x(x,y,t) for this trajectory KA y(t)' at 1,30; set palette rgbformulae 33,13,10; set pm3d map; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; set cbrange[-2e-7<*:*<1.5e-7]; set colorbox user origin 0.48,0.06 size 0.008,0.22; splot 'DATA_Full_WF_G_J_CN.txt' index t/2 using 1:2:4 title 'Full WF R^4 J_x(x,y,t)'; unset pm3d; unset view; unset colorbox; unset key;
        set origin 0.5, 0.66; set size 0.25, 0.33; clear; set key title '' inside; set xrange [$x1min:$x1max]; set yrange [0:0.4]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; plot for [i=2:(($jmax+1)/2+1)] 'DATA_chiInfo_XO.txt' index t/2 using (posx1[\$0+1]):i title sprintf('KA |Chi^{%d}(x)|',i-1) w l, for [i=2:(($jmax+1)/2+1)] 'DATA_chiInfo_CN.txt' index tCN using (posx1CN[\$0+1]):i title sprintf('CN |Chi^{%d}(x)|',i-1) w l, for [i=2:(($jmax+1)/2+1)] 'DATA_chiInfo_CN_w_KA_trajs.txt' index t/2 using (posx1CN[\$0+1]):i title sprintf('CN with KA traj condit wfs |Chi^{%d}(x)|',i-1) w l;
        set origin 0.5, 0.33; set size 0.25, 0.33; clear; set xtics auto; set xrange [$x1min:$x1max]; set yrange [0:0.013]; set key title 'x Conditional Wave Function Probability Density' inside; set xlabel 'Position x'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_2D_XO_KinAdv_BornHuang_tINDEP.txt' index t using (posx1[\$0+1]):1 title columnheader(1) w l, 'DATA_probDensity_conditional_WF_CN.txt' index t using (posx1CN[\$0+1]):1 title columnheader(1) w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' using 1:3 every ::(t/2)::(t/2+2) title 'x component of the Trajectory' w points lt 5 pt 6 ps 1 lc rgb 'red';
        set origin 0.51, 0.0; set size 0.26, 0.33; clear; set key title ' ' inside; set xrange [$x1min:$x1max]; set xlabel 'Position x'; set ylabel 'G(x)'; set y2tics nomirror; set y2label 'J(x)'; set yrange [-1e-5<*:*<1e-5]; set y2range [-1e-5<*:*<1e-5]; plot 'DATA_G_J_with_CN_and_KA_trajs.txt' index t/2 using (posx1CN[\$0+1]):1 title 'R^4 G(x) CN with KA trajs' axes x1y1 w l lc rgb 'red', 'DATA_G_J_with_CN_and_KA_trajs.txt' index t/2 using (posx1CN[\$0+1]):2 title 'R^4 J(x) CN with KA trajs' axes x1y2 w l lc rgb 'blue'; unset y2tics; unset y2label; unset key;
        set origin 0.75, 0.66; set size 0.25, 0.33; clear; set key title '' inside; set xrange [$x1min:$x1max]; set yrange [0:0.4]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; if($jmax+1>2){plot for [i=(($jmax+1)/2+2):($jmax+1)] 'DATA_chiInfo_XO.txt' index t/2 using (posx1[\$0+1]):i title sprintf('KA |Chi^{%d}(x)|',i-1) w l, for [i=(($jmax+1)/2+2):($jmax+1)] 'DATA_chiInfo_CN.txt' index tCN using (posx1CN[\$0+1]):i title sprintf('CN |Chi^{%d}(x)|',i-1) w l, for [i=(($jmax+1)/2+2):($jmax+1)] 'DATA_chiInfo_CN_w_KA_trajs.txt' index t/2 using (posx1CN[\$0+1]):i title sprintf('CN with KA traj condit wfs |Chi^{%d}(x)|',i-1) w l;};
        set origin 0.75, 0.33; set size 0.25, 0.33; clear; set xrange [$x2min:$x2max]; set yrange [0:0.013]; set key title 'y Conditional Wave Function Probability Density';  set xlabel 'Position y'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_2D_XO_KinAdv_BornHuang_tINDEP.txt' index (t+1) using (posx2[\$0+1]):1 title columnheader(1) w l, 'DATA_probDensity_conditional_WF_CN.txt' index (t+1) using (posx2CN[\$0+1]):1 title columnheader(1) w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$k0.txt' using 2:3 every ::(t/2)::(t/2+2) title 'y component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.76, 0.0; set size 0.25, 0.33; clear; set key title ' ' inside; set xrange [$x1min:$x1max]; set yrange [-0.1:0.2]; set xlabel 'Position x'; set ylabel 'G(x)'; set y2tics nomirror; set y2label 'J(x)'; set y2range [-0.1:0.1]; plot 'DATA_G_J_x_KA.txt' index t/2 using (posx1[\$0+1]):1 title 'G(x) KA aprox' axes x1y1 w l lc rgb 'red', 'DATA_G_J_x_KA.txt' index t/2 using (posx1[\$0+1]):2 title 'J(x) KA aprox' axes x1y2 w l lc rgb 'blue'; unset y2tics; unset y2label; unset key;
        t=t+2; tCN=tCN+1; if(tCN>($tIts)){tCN=0;}; unset multiplot;}"

        echo "Done!"

        # ta goikoan plot for [i=1:($jmax/2+2)]

        #set origin 0, 0.5; set size 0.35, 0.5; clear; set key title; set xrange [$x1min:$x1max]; set yrange [0:0.3]; set xlabel 'Position x'; set ylabel '|Chi^j(x)|'; set y2tics nomirror; set y2label 'J(x)'; set y2range [-0.45:0.2]; plot for [i=($jmax/2+1):($jmax)] 'DATA_chiInfo.txt' index t/2 using (posx1[\$0+1]):i title sprintf('|Chi^{%d}(x)|',i-1) axes x1y1 w l, 'DATA_G_J_x.txt' index t/2 using (posx1[\$0+1]):2 title 'J(x)' axes x1y2 w l; unset y2tics; unset y2label;

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
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."

        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output '../ANIMATION_GIFS/$expLabel XO_NoGJ_$currentDateTime x,yDivs($nx1,$nx2) x,ymin($x1min,$x2min) x,ymax($x1max,$x2max) nTrjs($numTrajs) dt($dt) tIt($numIt) mx,y($mass1,$mass2) outEvry($outputEvery).gif'; t=0; countTrajs=0; tmax=2*$numIt*$numTrajs/$outputEvery; dx1=($x1max-($x1min))/$nx1; dx2=($x2max-($x2min))/$nx2; array posx1[$nx1+1]; array posx2[$nx2+1]; do for [i=1:($nx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($nx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){  set multiplot;
        set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.025]; set key title 'x Conditional Wave Function Probability Density'; set xlabel 'Position x'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index t using (posx1[\$0+1]):1 title columnheader(1) w l, 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' every ::(t/2)::(t/2+2) using 1:3 title 'x component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.025];set key title 'y Conditional Wave Function Probability Density';  set xlabel 'Position y'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' index (t+1) using (posx2[\$0+1]):1 title columnheader(1) w l, 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' every ::(t/2)::(t/2+2) using 2:3 title 'y component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red'; set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x'; set ylabel 'Position y'; if ((t-(2*$numIt*countTrajs)) >= 2*$numIt){ countTrajs=countTrajs+1; }; set pm3d map; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt' title 'Potential Energy Map', 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$k0.txt' every ::1::(t/2+2) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
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

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XOKA_CNTrajs_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:i title sprintf('CN Trajects kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs CN_Trajs'; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:i title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XOKA_CNArea_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = $numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('CN Probability Density kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs CN_ProbabDens'; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:a[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XOKA_XONoGJ_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 2*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('SPCWF using KinAdv kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference SPCWF KinAdv vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:t[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XONoGJ_CNTrajs_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = 4*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:i title sprintf('CN Trajects kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference CN Trajs vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:a[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"

    gnuplot -e "set terminal png size 1800,500; set output '../ANIMATION_GIFS/$label Transmitace_XONoGJ_CNArea_gc_time.png'; array t[$numk+1]; do for [i=1:($numk)] { t[i] = $numk+i }; array b[$numk+1]; do for [i=1:($numk)] { b[i] = 4*$numk+i }; array a[$numk+1]; do for [i=1:($numk)] { a[i] = 3*$numk+i }; tmax1=$numk*4+1; tmax2=$numk*5+1; set multiplot; set origin 0.0,0.0; set size 0.5, 1.0; set yrange [0:0.4]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Aproximated Transmited Density'; plot for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:t[i] title sprintf('CN Prob Density kx_0= k_{%d}',i) w l, for [i=1:($numk)] 'DATA_porpsToPlot.txt' using tmax1:a[i] title sprintf('SPCWF No G,J kx_0= k_{%d}',i) w l; set origin 0.5, 0.0; set size 0.5, 1.0; clear; set yrange [0:0.1]; set xrange [0:$tIts]; set xlabel 'time (a.u.)'; set ylabel 'Absolute Difference CN Prob Denst vs SPCWF No GJ '; plot for [i=1:($numk)] 'DATA_errorsToPlot.txt' using tmax2:b[i] title sprintf('Initial k = k_{%d}',i) w l; unset multiplot;"
