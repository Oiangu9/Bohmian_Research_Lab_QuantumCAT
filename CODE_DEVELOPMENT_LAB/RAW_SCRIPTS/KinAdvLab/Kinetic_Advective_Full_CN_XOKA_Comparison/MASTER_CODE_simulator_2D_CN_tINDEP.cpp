$$14$1$CODE_simulator_2D_CN_tINDEP.cpp$$

// SCHRODINGER EQUATION SOLVER 2D
// for TIME INDEPENDENT POTENTIALS
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
//We define the initial WF
cdouble psi0(double x, double y){
 $1$;
}
//We define the Potential energy field TIME INDEPENDENT
double V(double x, double y){
 $2$;
}
int main(){
    //defining the particle parameters
    double m1=$11$, m2=$12$, hbar=$13$;
    //define the SPACE GRID
    double x1min=$5$, x1max=$6$, x2min=$7$, x2max=$8$;
    int nx1=$3$, nx2=$4$, outputDataEvery=$14$;
    double dx1=(x1max-x1min)/nx1, dx2=(x2max-x2min)/nx2;
    int gridPoints=(nx1+1)*(nx2+1);
    //define the TIME GRID
    //double tmin=0.0;
    double dt=$9$;
    int numIt=$10$;
    //Build the planed initial WavePacket (which should be seen as a 2d array
    VectorXcd psi_ini(gridPoints); //nx, ny intervalos, then nx+1, ny+1 ptos por cada dim
    for(int ix=0; ix<=nx1; ix++){
        for (int iy=0; iy<=nx2;iy++){
            psi_ini(ix*(nx2+1)+iy)=psi0(x1min+ix*dx1, x2min+iy*dx2);
        }
    }
    //Build the propagators U1 and U2
    SparseMatrix<cdouble> UL(gridPoints, gridPoints);
    UL.reserve(VectorXi::Constant(gridPoints,5));
    cdouble a=J*dt*hbar/(4*m1*dx1*dx1), d=J*dt*hbar/(4*m2*dx2*dx2);
    int j;
    for(int ix=0;ix<=nx1;ix++){
        for(int iy=0; iy<=nx2;iy++){
            j=ix*(nx2+1)+iy;
            UL.insert(j,j)=1.0+J*dt*(hbar*hbar*(1/(m1*dx1*dx1)+1/(m2*dx2*dx2))+V(x1min+ix*dx1,x2min+iy*dx2))/((cdouble)2*hbar);
            if(iy!=nx2) UL.insert(j,j+1)=-d;
            if(iy!=0) UL.insert(j,j-1)=-d;
            if(ix>0) UL.insert(j,j-nx2-1)=-a;
            if(ix<nx1) UL.insert(j,j+nx2+1)=-a;
        }
    }
    UL.makeCompressed();
    // UR is exactly de complex conjugate of UL
    SparseMatrix<cdouble> UR(gridPoints, gridPoints);
    UR.reserve(VectorXi::Constant(gridPoints,5));
    UR = UL.conjugate();
    UR.makeCompressed();
    // Make the LU decomposition of UL
    SparseLU<SparseMatrix<cdouble>> LUsolver;
    LUsolver.compute(UL);
    if(LUsolver.info()!=Success) {
        cout << "LU decomposition FAILED!" << endl;
        return 1;
    }
    //we prepare everything for the loop
    VectorXcd URpsi(gridPoints), psiNext(gridPoints);
    psiNext = psi_ini;
    ofstream outsideFile;
    outsideFile.open("DATA_rawSimulationData_2D_CN.txt");
    outsideFile << nx1 << " " << nx2 << " " << numIt << " " << dt << " " << x1min << " " << x1max << " " << x2min << " " << x2max << " " << m1 << " " << m2 << " "<<hbar<< " (Spatial divisions in x1, x2, time iterations, x1min, x1max, x2min, x2max) CN 2D" << endl;
    outsideFile << std::setprecision(17);
    for(int it=0; it<=numIt; it++){
        URpsi= UR*psiNext;
        psiNext = LUsolver.solve(URpsi);
        if(it%outputDataEvery == 0){
            outsideFile << "it"<<it<<endl<< psiNext<<endl;
        }
    }
    outsideFile.close();
    return 0;
}
