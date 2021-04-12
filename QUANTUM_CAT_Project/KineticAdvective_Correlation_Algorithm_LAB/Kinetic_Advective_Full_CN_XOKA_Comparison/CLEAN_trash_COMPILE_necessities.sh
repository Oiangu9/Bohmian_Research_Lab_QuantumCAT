#!/bin/bash

cd SOURCE_FILES
rm DATA_chiInfo_CN_w_KA_trajs.txt
rm DATA_chiInfo_XO.txt
rm DATA_chiInfo_CN.txt
#rm DATA_CN_Im_Uj.txt
#rm DATA_CN_Re_Uj.txt
#rm DATA_XO_Im_Uj.txt
#rm DATA_XO_Re_Uj.txt
rm DATA_G_J_x_KA.txt
rm DATA_G_J_with_CN_and_KA_trajs.txt
rm DATA_KinAdv_x.txt
rm DATA_probDensity_conditional_WF_CN.txt
rm DATA_probDensity_WF_CN.txt
rm DATA_sumChiInfo_CN_with_KA_traj.txt
rm DATA_sumChiInfo_XO.txt
rm DATA_sumChiInfo_CN.txt
rm DATA_rawSimulationData_2D_CN.txt
rm DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt
rm DATA_probabilityToPlot_2D_XO_KinAdv_BornHuang_tINDEP.txt
rm DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt
rm DATA_trajectoriesToPlot*
rm DATA_Full_WF_G_J_CN.txt
rm DATA_derivatives_S_R.txt

rm CODE_CN2D_ChiCalculator.cpp
rm CODE_simulator_*
rm CODE_Simulator_*
rm *.png

rm EXE_*

g++ -Wall -O CODE_codeFileGenerator_2D_CN_tINDEP.cpp -o EXE_codeFileGenerator_2D_CN_tINDEP

g++ -Wall -O CODE_codeFileGenerator_CN2D_ChiCalculator.cpp -o EXE_codeFileGenerator_CN2D_ChiCalculator

g++ -Wall -O CODE_codeFileGenerator_2D_XO_KINADV_BornHuang_tINDEP.cpp -o EXE_codeFileGenerator_2D_XO_KINADV_BornHuang_tINDEP

g++ -Wall -O3 CODE_proportionUnifier_ErrorCalculator.cpp -o EXE_proportionUnifier_ErrorCalculator

g++ -Wall -O CODE_codeFileGenerator_CN_Conditional_WFs_with_XO_Trajectories.cpp -o EXE_codeFileGenerator_CN_Conditional_WFs_with_XO_Trajectories
