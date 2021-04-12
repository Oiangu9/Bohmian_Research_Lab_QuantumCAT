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
#include <random>
using namespace std::complex_literals;
using namespace Eigen;
using namespace std;
#define PI 3.141592653589793238463
#define INF 1000000.0
int nx=500, ny=700, jmax=10;
double xmin=-50.0, xmax=70.0, ymin=-20.0, ymax=20.0, dx, dy;
cdouble EigenstatesForSectionsInx(double y, double x, int j){ double g=3.0,a1=5.0,a2=40.0,Lmax=40.0,Lmin=5.0,o, L; o= -(Lmax-Lmin)/2.0*(1.0/(1.0+exp((x-a1)/g))+1.0/(1.0+exp((-x+a2)/g)))-Lmin/2.0; L=-2*o; if(y>(o+L) || y<o){return 0.0;} else{return sqrt(2.0/L)*sin(j*PI*(y-o)/L);}}
double potential(double x, double y){double g=3.0,a1=5.0,a2=40.0,Lmax=40.0,Lmin=5.0,L1; L1= (Lmax-Lmin)/2.0*(1.0/(1.0+exp((x-a1)/g))+1.0/(1.0+exp((-x+a2)/g)))+Lmin/2.0; if(y>L1 || y<-L1){return 100.0;} else{return 0.0;}}
ArrayXd posx(nx+1), posy(ny+1);
int main(){
ofstream plotFile, potentialFile;
plotFile.open("DATA_adiabaticStatePlot.txt");
//Prepare the grid in x y
dx=(xmax-xmin)/nx;
dy=(ymax-ymin)/ny;
for(int ix=0; ix<=nx; ix++){
posx(ix) = xmin+ix*dx;
}
for(int iy=0; iy<=ny; iy++){
posy(iy) = ymin+iy*dy;
}
for(int ix=0; ix<=nx; ++ix){
for(int iy=0; iy<=ny; ++iy){
plotFile << posy(iy)<<" ";
for(int j=0; j<=jmax; ++j){
plotFile << real(EigenstatesForSectionsInx(posy(iy), posx(ix), j))<<" ";
}
plotFile<<endl;
}
plotFile<<endl<<endl; }
plotFile.close();
potentialFile.open("DATA_potentialToPlot.txt");
for(int i=0; i<=nx; i+=1){
for(int j=0; j<=ny; j+=1){
potentialFile << posx(i) << " " << posy(j) << " " << potential(posx(i), posy(j))<< endl;
}potentialFile<<endl;}
potentialFile.close();
return 0;
}
