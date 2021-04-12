#!/bin/bash

rm DATA_chiInfo.txt
rm DATA_plotWFInfo.txt
rm DATA_sumChiInfo.txt
rm DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt
rm DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt
rm DATA_rawSimulationData_2D_CN.txt
rm DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt

rm CODE_CN2D_ChiCalculator.cpp
rm CODE_CN2D_simulator_2D_CN_tINDEP.cpp
rm CODE_XO-KA_simulator_2D_XO_KINADV_BornHeun_tINDEP.cpp


rm EXE_CN2D_ChiCalculator
rm EXE_CN2D_codeFileGenerator_2D_CN
rm EXE_CN2D_GeneratorChiCalculator
rm EXE_CN2D_simulator_2D_CN_tINDEP
rm EXE_XO-KA_codeFileGenerator_2D_XO_KINADV_BornHeun
rm EXE_XO-KA_simulator_2D_XO_KINADV_BornHeun_tINDEP

g++ -Wall -O CODE_CN2D_codeFileGenerator_2D_CN.cpp -o EXE_CN2D_codeFileGenerator_2D_CN

g++ -Wall -O CODE_CN2D_GeneratorChiCalculator.cpp -o EXE_CN2D_GeneratorChiCalculator

g++ -Wall -O CODE_XO-KA_codeFileGenerator_2D_XO_KINADV_BornHeun.cpp -o EXE_XO-KA_codeFileGenerator_2D_XO_KINADV_BornHeun

g++ -Wall -O3 CODE_proportionUnifier_ErrorCalculator.cpp -o EXE_proportionUnifier_ErrorCalculator

