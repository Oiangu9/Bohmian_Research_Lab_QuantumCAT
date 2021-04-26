$$17$1$CODE_CN2D_ChiCalculator.cpp$$

// This is a programm intended to calculate the chi coefficients given a 2D input waveFunction and the adiabatic states for the sections as analytic expressions
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
//USER INPUT------
//The information for the input 2D wavefunction psi(x,y)

ifstream inputFileWF;

int nx1=$4$,nx2=$5$, numIt=$10$, numTrajs=$12$, gridPoints, aux, wholex1, wholex2;
double x1min=$6$,x1max=$7$,x2min=$8$,x2max=$9$,xBound=$11$,dx1,dx2, normWF, sumCorners, sumBorders, sumInterior, sumaParaChisx, hbar=$14$, m1=$15$, m2=$16$, fractionalx1, fractionalx2, wholex1f, wholex2f, dt=$17$;
double xmin=$6$,xmax=$7$,ymin=$8$,ymax=$9$;
double * trajectoriesx1, * trajectoriesx2;
//The adiabatic states for sections in x and the maximum desired j to calculate
int jmax=$3$;
cdouble integral;
ArrayXXcd Chijx_container(nx1+1, jmax+1); // Each column k contains for the energy level j=k the value of chi as a function of the section in x (the position in x)
ArrayXd sumaChisHasta(jmax+1);
cdouble EigenstatesForSectionsInx(double y, double x, int j){$2$}

//We initialize the transmitance counter vector

ArrayXd trajProportionCrossed = ArrayXd::Zero(numIt+1);
ArrayXd areaProportionCrossed = ArrayXd::Zero(numIt+1);
//WF EXTRACTION FROM FILE-------
int main(){
inputFileWF.open("$1$"); //the path of the file where the WF is saved in a sequentialized matrix vector fashion
string line;
cdouble arrayEl;
getline(inputFileWF,line);
//istringstream ss(line);
//ss >> nx1 >> nx2 >> nt >> dt >>x1min >> x1max >> x2min >> x2max>>m1>>m2>>hbar; //nx is the number of space division, which is one less than the number of spatial points

trajectoriesx1=(double*) malloc(numTrajs*sizeof(double));
trajectoriesx2=(double*) malloc(numTrajs*sizeof(double));
gridPoints=(nx1+1)*(nx2+1);
ArrayXcd WF(gridPoints);

ArrayXcd auxComplexVectorx1(gridPoints), auxComplexVectorx2(gridPoints), conjPsi(gridPoints);
ArrayXd vFieldx1(gridPoints), vFieldx2(gridPoints);
ArrayXd posx1(gridPoints), posx2(gridPoints), absWF(gridPoints);
MatrixXd results(nx2+1, 3);
//Prepare the grid in x y
dx1=(x1max-x1min)/nx1;
dx2=(x2max-x2min)/nx2;
for(int ix=0; ix<=nx1; ix++){
    for(int iy=0; iy<=nx2; iy++){
        posx1(ix*(nx2+1)+iy) = x1min+ix*dx1;
        posx2(ix*(nx2+1)+iy) = x2min+iy*dx2;
    }
}
//We define the initial trajectories randomly according to the probability density of the initial wf
for(int xIt=0; xIt<gridPoints; xIt++){
    getline(inputFileWF,line);
    istringstream ss(line);
    ss>>arrayEl;
    WF(xIt)=arrayEl;
}
absWF=100*abs2(WF);
double* probabClist = absWF.data();
std::default_random_engine generator;
std::discrete_distribution<int> distribution (probabClist,probabClist+gridPoints-1);
for(int i=0; i<numTrajs; i++){
    aux=distribution(generator);
    trajectoriesx1[i] = posx1(aux);
    trajectoriesx2[i] = posx2(aux);
}
inputFileWF.seekg(0);
//Output files initialized
ofstream chiInfo, sumChiInfo, wfPlot, trajFile, trajProps, areaProps;
chiInfo.open("DATA_chiInfo_CN.txt");
sumChiInfo.open("DATA_sumChiInfo_CN.txt");
wfPlot.open("DATA_probDensity_WF_CN.txt");
trajFile.open("DATA_CN_Trajs_k=$13$.txt");
//BEGINING OF TIME ITERATIONS!--
getline(inputFileWF,line);
istringstream ss(line);
for(int tIt=0; tIt<=numIt; ++tIt){
    //Actually extract the wf
    getline(inputFileWF,line);
    for(int xIt=0; xIt<gridPoints; xIt++){
        getline(inputFileWF,line);
        istringstream ss(line);
        ss>>arrayEl;
        WF(xIt)=arrayEl;
    }
    absWF=abs2(WF);
    //CAlculate the velocity fields in x and y
    conjPsi= conj(WF);


    // psi^-1* diff(psi,x)/ absWF is computed
    for(int i=0; i<=nx1; i++){
        for(int j=0; j<nx2;j++){
            if(i==0){
                auxComplexVectorx1(j)=conjPsi(j)*(WF(nx2+1+j)-WF(j))/(absWF(j)*dx1);
            } else if(i==nx1){
                auxComplexVectorx1(nx1*(nx2+1)+j)=conjPsi(nx1*(nx2+1)+j)*(WF(nx1*(nx2+1)+j)-WF((nx1-1)*(nx2+1)+j))/(absWF(nx1*(nx2+1)+j)*dx1);
            } else if(i==1 || i==(nx1-1)){
                auxComplexVectorx1(i*(nx2+1)+j)=conjPsi(i*(nx2+1)+j)*(WF((i+1)*(nx2+1)+j)-WF((i-1)*(nx2+1)+j))/(2.0*absWF(i*(nx2+1)+j)*dx1);
            } else{
                auxComplexVectorx1(i*(nx2+1)+j)=conjPsi(i*(nx2+1)+j)*(-WF((i+2)*(nx2+1)+j)+8.0*WF((i+1)*(nx2+1)+j)-8.0*WF((i-1)*(nx2+1)+j)+WF((i-2)*(nx2+1)+j))/(12.0*absWF(i*(nx2+1)+j)*dx1);
            }
            if(j==0){
                auxComplexVectorx2(i*(nx2+1))=conjPsi(i*(nx2+1))*(WF(i*(nx2+1)+1)-WF(i*(nx2+1)))/(absWF(i*(nx2+1))*dx2);
            } else if(j==nx2){
                auxComplexVectorx2(i*(nx2+1)+nx2)=conjPsi(i*(nx2+1)+nx2)*(WF(i*(nx2+1)+nx2)-WF(i*(nx2+1)+nx2-1))/(absWF(i*(nx2+1)+nx2)*dx2);
            } else if(j==1 || j==(nx2-1)){
                auxComplexVectorx2(i*(nx2+1)+j)=conjPsi(i*(nx2+1)+j)*(WF(i*(nx2+1)+j+1)-WF(i*(nx2+1)+j-1))/(2.0*absWF(i*(nx2+1)+j)*dx2);
            }else{
                auxComplexVectorx2(i*(nx2+1)+j)=conjPsi(i*(nx2+1)+j)*(-WF(i*(nx2+1)+j+2)+8.0*WF(i*(nx2+1)+j+1)-8.0*WF(i*(nx2+1)+j-1)+WF(i*(nx2+1)+j-2))/(12.0*absWF(i*(nx2+1)+j)*dx2);
            }
        }
    }
    // imaginary part is extracted and Jk obtained
    vFieldx1 = (hbar/m1)*imag(auxComplexVectorx1);
    vFieldx2 = (hbar/m2)*imag(auxComplexVectorx2);
    for(int i=0; i<numTrajs;i++){
        //we apply the discretisation of the grid to the traj positions
        fractionalx1 = std::modf((trajectoriesx1[i]-x1min)/dx1, &wholex1f);
        fractionalx2 = std::modf((trajectoriesx2[i]-x2min)/dx2, &wholex2f);
        wholex1 = wholex1f;
        wholex2 = wholex2f;
        if(wholex1>=nx1-1){wholex1=nx1-2;}else if(wholex1<0){wholex1=0;}
        if(wholex2>=nx2-1){wholex2=nx2-2;}else if(wholex2<0){wholex2=0;}
        trajectoriesx1[i] = trajectoriesx1[i]+( (1-fractionalx1)*vFieldx1( wholex1*(nx2+1) + wholex2)+fractionalx1*vFieldx1( (wholex1+1)*(nx2+1) + wholex2 ))*dt;
        trajectoriesx2[i] = trajectoriesx2[i]+( (1-fractionalx2)*vFieldx2( wholex1*(nx2+1) + wholex2)+fractionalx2*vFieldx2( wholex1*(nx2+1) + wholex2+1 ))*dt;
        trajFile << trajectoriesx1[i] << " " << trajectoriesx2[i] << endl;
         if (trajectoriesx1[i]>=xBound){trajProportionCrossed(tIt)+=1;}
    }
    trajFile<<endl<<endl;
    // calculate the normWF of the WF using the 2d trapezium rule
    sumCorners = 0;
    sumBorders = 0;
    sumInterior = 0;
    sumCorners = absWF(0) + absWF(nx2) + absWF(gridPoints-1) + absWF(nx1*(nx2+1));
    for(int ix=1; ix<nx1; ++ix){
        for(int iy=1; iy<nx2; ++iy){
            sumInterior += absWF(ix*(nx2+1)+iy);
        }
    }
    for(int iy=1; iy<nx2; ++iy){
        sumBorders += absWF(iy) + absWF(nx1*(nx2+1) + iy); // ix=0 row and ix=nx1 row
    }
    for(int ix=1; ix<nx1; ++ix){
        sumBorders += absWF(ix*(nx2+1)) + absWF(ix*(nx2+1)+nx2); //iy=0 row and iy=nx2 row
    }
    normWF=0.25*dx1*dx2*(sumCorners + 2*sumBorders +4* sumInterior);
    for(int ix=0; ix<=nx1; ++ix){
        for(int iy=0; iy<=nx2; ++iy){
            if(posx1(ix*(nx2+1)+iy)>=xBound){
                areaProportionCrossed(tIt)+=absWF(ix*(nx2+1)+iy);
            }
        }
    }
    areaProportionCrossed(tIt)=areaProportionCrossed(tIt)/normWF;
    //CALCULATE PROBABILITY AREA PASSING xBOUND inside the for loop with j and i in the next section
    /*
    //CACLCULATE CHI_j(x,t)-----------
    sumaParaChisx=0;
    for(int j=0; j<=jmax; ++j){
        for(int ix=0; ix<=nx1; ++ix){
            // for each section in x we integrate the WF in y restricted to this x
            integral=0.5*(WF(ix*(nx2+1))+WF(ix*(nx2+1)+nx2));
            for(int iy=1; iy<nx2; ++iy){

            integral += EigenstatesForSectionsInx(x2min+iy*dx2, x1min+ix*dx1, j)*WF(ix*(nx2+1)+iy);
            }
            Chijx_container(ix, j)=integral*dx2;
        }
        //we calculate the incremental sum of the moduluous of chi integrated
        sumaParaChisx+=abs2(Chijx_container.col(j)).sum();
        sumaChisHasta(j)=sumaParaChisx*dx1;
    } //OUTPUT Obtained CHI information for plots--------

    //chiInfo << "x_position";
    for(int j=0; j<=jmax; ++j){
        //chiInfo << " "<< sumaChisHasta(j); //header like the sum of the chis will be printed until the given j
        sumChiInfo<<j<<" "<<sumaChisHasta(j)<<endl;
    }
    sumChiInfo<<endl<<endl;

    chiInfo << abs(Chijx_container) <<endl <<endl <<endl;
    */
    //OUTPUT WF information for plots-------------

    wfPlot <<"Norm="<< normWF << endl;
    //output the position together with the probability
    for(int ix=0; ix<=nx1; ix++){
        results << posx1.block(ix*nx2+ix,0, nx2+1, 1), posx2.block(ix*nx2+ix,0, nx2+1, 1), absWF.block(ix*nx2+ix,0, nx2+1, 1);
        wfPlot <<results<<endl<< endl;
    }
    wfPlot <<endl;

}
wfPlot.close();
chiInfo.close();
sumChiInfo.close();
trajFile.close();
inputFileWF.close();
areaProportionCrossed= dx1*dx2*areaProportionCrossed;

trajProps.open("DATA_CN_trajProps_k=$13$.txt");
areaProps.open("DATA_CN_areaProps_k=$13$.txt");
trajProportionCrossed = (1/(double)numTrajs)*trajProportionCrossed;
trajProps<<trajProportionCrossed<<endl;
areaProps<<areaProportionCrossed<<endl;
trajProps.close();
areaProps.close();

return 0;
}
