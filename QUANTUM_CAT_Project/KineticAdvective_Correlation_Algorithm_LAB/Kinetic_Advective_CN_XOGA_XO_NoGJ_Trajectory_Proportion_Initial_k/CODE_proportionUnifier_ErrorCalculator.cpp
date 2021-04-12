// ERRORS IN PROPORTIONS!!!
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <complex>
#include "LIB_dcomplex.h" // Macro for using dcomplex as std::complex<double> and J as the complex 0.0+i*1.0
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <cmath>
using namespace std::complex_literals;
using namespace Eigen;
using namespace std;
const double PI = 3.141592653589793238463; const double EULER=2.718281828459045;

int main(int argNum, char **argVec){
    // First argument is the number of different k-s you are inputting
    // next the number of time its
    // and next the sequence of k-s
    int numberOfks, numIts;
    sscanf(argVec[1], "%d", &numberOfks);
    sscanf(argVec[2], "%d", &numIts);
    ifstream propFiles;
    
    ArrayXXd trajProps(numIts+1, numberOfks*4+1);
    ArrayXXd CN_trajProps(numIts+1, numberOfks), CN_areaProps(numIts+1,numberOfks), XO_KinAdv_trajProps(numIts+1, numberOfks), XO_NoGJ_trajProps(numIts+1, numberOfks);
    
    for(int it=0; it<=numIts; ++it){trajProps(it, numberOfks*4)=it;}
    
    string line;
    double arrayEl;
    
    std::ostringstream s;

    
    for(int k=0; k<numberOfks; ++k){
        
        s<< "DATA_CN_trajProps_k="<<argVec[k+3]<<".txt";
        //std::cout << s.str();
        propFiles.open(s.str());
        
        for(int it=0; it<=numIts; ++it){
            getline(propFiles,line);
            istringstream ss(line);
            ss>>arrayEl;
            trajProps(it, k) = arrayEl;
        }
        
        propFiles.close();
        s.str("");
        s.clear();
        
        s << "DATA_CN_areaProps_k="<<argVec[k+3]<<".txt";
        propFiles.open(s.str());
        
        for(int it=0; it<=numIts; ++it){
            getline(propFiles,line);
            istringstream ss(line);
            ss>>arrayEl;
            trajProps(it, numberOfks+k) = arrayEl;
        }
        
        propFiles.close();
        s.str("");
        s.clear();
        
        s << "DATA_XO_KA_trajProps_k="<<argVec[k+3]<<".txt";
        propFiles.open(s.str());
        
        for(int it=0; it<=numIts; ++it){
            getline(propFiles,line);
            istringstream ss(line);
            ss>>arrayEl;
            trajProps(it, numberOfks*2+k) = arrayEl;
        }
        
        propFiles.close();
        s.str("");
        s.clear();
        
        s << "DATA_XO_NoGJ_trajProps_k="<<argVec[k+3]<<".txt";
        propFiles.open(s.str());
        
        for(int it=0; it<=numIts; ++it){
            getline(propFiles,line);
            istringstream ss(line);
            ss>>arrayEl;
            trajProps(it, numberOfks*3+k) = arrayEl;
        }
        
        propFiles.close();
        s.str("");
        s.clear();
        
    }
    CN_trajProps = trajProps.block(0,0, numIts, numberOfks);
    CN_areaProps = trajProps.block(0, numberOfks, numIts, numberOfks);
    
    XO_KinAdv_trajProps = trajProps.block(0, 2*numberOfks, numIts, numberOfks);
    
    XO_NoGJ_trajProps = trajProps.block(0, 3*numberOfks, numIts, numberOfks);
    
    ofstream plotOutputFile, errorOutputFile, errorPlotfile;
    plotOutputFile.open("DATA_porpsToPlot.txt");
    errorOutputFile.open("DATA_errorsInfo.txt");
    errorPlotfile.open("DATA_errorsToPlot.txt");
    
    plotOutputFile << trajProps<<endl;
    
    ArrayXXd errors(numIts, 5*numberOfks+1);
    for(int it=0; it<numIts; ++it){errors(it, 5*numberOfks)=it;}
    
    for(int k=0; k<numberOfks; ++k){
        errors.col(k)=abs(XO_KinAdv_trajProps.col(k)-CN_trajProps.col(k));
        errors.col(numberOfks+k) = abs(XO_KinAdv_trajProps.col(k)-CN_areaProps.col(k));
        errors.col(2*numberOfks+k) = abs(XO_KinAdv_trajProps.col(k)-XO_NoGJ_trajProps.col(k));
        
        errors.col(3*numberOfks+k) = abs(XO_NoGJ_trajProps.col(k) - CN_trajProps.col(k));
        errors.col(4*numberOfks+k) = abs(XO_NoGJ_trajProps.col(k) - CN_areaProps.col(k));
    }
    errorPlotfile << errors<<endl;
    
    errorOutputFile << errors.colwise().maxCoeff();
    
    plotOutputFile.close();
    errorOutputFile.close();
    errorPlotfile.close();
    return 0;
}
    
    

    
