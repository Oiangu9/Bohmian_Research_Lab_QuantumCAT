#!/bin/bash

psiIni=""
potential=""
nx1=1
nx2=1
x1min=1.0
x1max=1.0
x2min=1.0
x2max=1.0
dt=1.0
tIt=1
plotOpt=1
numTrajs=0
mass1=1.0
mass2=1.0
hbar=1.0

echo ""
echo " (1) DATA INPUT"
echo "Introduce in C++ syntax the initial WaveFunction as a function of space variables x1 and x2 please"
read trash
read psiIni
echo "Now introduce the Potential Energy field as a function of x1 and x2"
read trash
read potential
echo "Introduce the number of space divisions for the grid in x1 and x2"
read trash
read nx1 nx2
echo "Introduce x1min x1max x2min x2max"
read trash
read x1min x1max x2min x2max
echo " Introduce the mass at x1 and mass at x2 as m1 m2"
read trash
read mass1 mass2
echo " Introduce the value of hbar"
read trash
read hbar

if [[ $psiIni != *"return"* ]]; then
    psiIni="return ${psiIni}"
fi
if [[ $potential != *"return"* ]]; then
    potential="return ${potential}"
fi

echo " (2) CODE COMPILATION and EXECUTION"
echo " Generating code and compiling..."

g++ -Wall -O3 CODE_Generator_EnergyCalculator.cpp -o EXE_Generator_EnergyCalculator

./EXE_Generator_EnergyCalculator "$psiIni" "$potential" $nx1 $nx2 $x1min $x1max $x2min $x2max $mass1 $mass2 $hbar $dt
g++ -Wall -O3 CODE_2D_WF_EnergyCalculator.cpp -o EXE_2D_WF_EnergyCalculator
echo " Done!"
echo ""
echo " Executing calculations..."
START_TIME=$SECONDS
./EXE_2D_WF_EnergyCalculator
CALCULATION_TIME=$SECONDS
echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
echo ""
rm EXE_Generator_EnergyCalculator
rm EXE_2D_WF_EnergyCalculator

