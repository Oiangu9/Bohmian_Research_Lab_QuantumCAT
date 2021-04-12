$$34$1$CODE_simulator_XO_KinAdv.cpp$$

// 2D SCHRODINGER EQUATION SOLVER - XO ALGORITHM Kinetic and Advective Correlation Potential approximation of G and J:
// The Conditional Single Particle Wave Function (CSPWF) of the dimensions (x,y) will be evolved for each initial conditions using a 1D Cranck Nicolson method for the Pseudo Schrodinger Equations
//WARNING! The eigenstates for the sections are not set cdouble because in the box case they are fully real. However, this should be adjusted for a generalized algorithm. As such, complex conjugate should be done in the integrals for the Ujx!!!
//We include the necessary libraries
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
#define customTrajs $33$
//USER INPUT DECLARATION----------------------------------------------------------------------
//We declare the Spatial and Time grids - names are self-explaning
double xmin = $7$, xmax = $8$, ymin = $9$, ymax = $10$, t0=0.0, dt=$11$, posx, posy, Nx, Ny;
cdouble Uj;
int xDivs=$5$, yDivs = $6$, timeIts=$12$, xBound=$27$, aux;
double dx=(xmax-xmin)/xDivs;
double dy=(ymax-ymin)/yDivs;
int outputDataEvery = $16$; // data will be outputed to the external file only every x time iterations, if 1, every time iteration will be outputted
//The constants: hbar and the masses of each dimension are declared
double hbar=$15$, mx=$3$, my=$4$;
//We declare the initial full Wave Function (x,y,t=t0). It will be used in order to obtain the initial CSPWF by evaluating at each dimension the initial position of the other trajectory
cdouble initialFullWF(double x, double y){
$1$
}
//We declare the static external Potential Field
double W(double x,double y){
$2$
}
//The highest used adiabatic states' index in the Born Huang expansion's truncated version
int xjmax = $23$, yjmax = $24$, lastjUsedInItx=0, lastjUsedInIty=0;
//The tolerance for the sum squared of the chis, the algorithm will use in each time iteration as many chis as required to achieve this squared sum tolerance (or rather they will arrive the maximum tolerated j)
double chiSumTolerance=$26$;
double sumaParaChisx, sumaParaChisy;
ArrayXd sumaChisx(xjmax+1), sumaChisy(yjmax+1);
cdouble correlPot, advectiveCor, kineticCor;
//The "analytical" expression for the adiabatic states for each x as a function of y and j are expressed
double eigenstatesForSectionsInx(double y, double x, int j){ //the so called psi_x^j(y)
$17$
}
//The "analytical" expression for the adiabatic states for each y as a function of x and j are expressed
double eigenstatesForSectionsIny(double x, double y, int j){ //the so called psi_y^j(x)
$20$
}
//The analytical expression for the above functions' first and second derivatives
double diffyEigenstatesForSectionsInx(double y, double x, int j){ //the so called d psi_x^j(y)/dy
$18$
}
double diffyyEigenstatesForSectionsInx(double y, double x, int j){ //the so called d**2 psi_x^j(y)/dy**2
$19$
}
double diffxEigenstatesForSectionsIny(double x, double y, int j){ //the so called d psi_y^j(x)/dx
$21$
}
double diffxxEigenstatesForSectionsIny(double x, double y, int j){ //the so called d**2 psi_y^j(x)/dx**2
$22$
}
//We choose if we want the Uj for y dimension to be calculated as well
double b_y=$25$;
//We declare the matrix that will save the Uj(x,t) values for each x and j, in order to allow the calculation of the normalized Uj
ArrayXXcd Ujx_container(xDivs+1, xjmax+1), Ujy_container(yDivs+1, yjmax +1), Chijx_container(xDivs+1, xjmax+1), Chijy_container(yDivs+1, yjmax+1);
//ArrayXd Ujx_normFactors(xjmax+1), Ujy_normFactors(yjmax+1);
//We declare a matrix to store the computed correlation potential G(x)+i*J(x) for a certain time iteration in order to allow its plot. Another one for G(y)+i*J(y)
//We also declare a matrix to store the real and imaginary parts of the Kinetic and Advective correlation potentials approximated for x and for y. (4 columns -> Re{Kin}, Im{Kin}, Re{Adv}, Im{Adv}
ArrayXXd G_J_x(xDivs+1, 2), G_J_y(yDivs+1, 2), KinAdv_x(xDivs+1, 4), KinAdv_y(yDivs+1, 4);

//We initialize a vector that will number the trajectories that have crossed a certain boundary in each time
ArrayXd trajProportionCrossed = ArrayXd::Zero(timeIts+1);
//The grid positions in x and y are saved into a vector as they are heavily used in the computation of Uj
ArrayXd xgrid(xDivs+1), ygrid(yDivs+1);
int main(){
for(int k=0; k<=xDivs; ++k){xgrid(k)=xmin+k*dx;}
for(int k=0; k<=yDivs;++k){ygrid(k)=ymin+dy*k;}


int numTrajs=$13$; // we choose the number of trajectories that will be evolved
int gridPointsFullWF=(xDivs+1)*(yDivs+1), whole;
double fractional, wholef;
double* initialPosx = (double*) malloc(numTrajs*sizeof(double));
double* initialPosy = (double*) malloc(numTrajs*sizeof(double));

if(customTrajs==0){ //The initial positions of each trajectory that will be evolved using the algorithm are chosen according to the probability distribution given by the modulus squared of the initial wave function
    ArrayXcd initialFullPsi(gridPointsFullWF);
    ArrayXd probabDensity(gridPointsFullWF);
    // the initial state of the full wavefunction is generated in order to obtain its modulus squared in each point
    for(int i=0; i<=xDivs; ++i){
      for(int j=0; j<=yDivs; ++j){
        initialFullPsi(i*(yDivs+1) + j) = initialFullWF(xgrid(i), ygrid(j));
      }
    }
    // the probability associated with each space point is generated
    probabDensity =100*abs2(initialFullPsi);
    // the random number generator initialised
    double* probabClist = probabDensity.data();
    std::default_random_engine generator;
    std::discrete_distribution<int> distribution (probabClist,probabClist+gridPointsFullWF-1);
    // and the initial positions chosen according to the associated probabilities
    for(int i=0; i<numTrajs; i++){
      aux=distribution(generator); //returns the winner index of the prob vector -> we must revert the indexing to 2d indexes
      initialPosx[i] = xmin + ((int) aux/(yDivs+1))*dx;
      initialPosy[i] = ymin + (aux%(yDivs+1))*dy;
    }
}else{ // custom trajectory selection by user. Initial positions should be charged in arrays initialPosx, initialPosy

    $34$

}

// begin the time iterations for each evolved trajectory - UL, UR must be renamed in every iteration - as if it was a time dependant potential algorithm for a 1D particle
//we declare and prepare the propagator matrices for the Cranck Nicolson (CN) evolution
SparseMatrix<cdouble> U1x(xDivs+1, xDivs+1), U2x(xDivs+1, xDivs+1);
SparseMatrix<cdouble> U1y(yDivs+1, yDivs+1), U2y(yDivs+1, yDivs+1);
U1x.reserve(VectorXi::Constant(xDivs+1,3));
U2x.reserve(VectorXi::Constant(xDivs+1,3));
U1y.reserve(VectorXi::Constant(yDivs+1,3));
U2y.reserve(VectorXi::Constant(yDivs+1,3));
cdouble ax=J*dt*hbar/(4.0*mx*dx*dx);
cdouble ay=J*dt*hbar/(4.0*my*dy*dy);
for(int i=1;i<xDivs;i++){
    U1x.insert(i,i)= 1.0*J; //just initialise to some random variable
    U1x.insert(i,i+1)= -ax;
    U1x.insert(i,i-1)= -ax;
}
U1x.insert(0,0)= 1.0*J;
U1x.insert(0,1)= -ax;
U1x.insert(xDivs,xDivs)= 1.0*J;
U1x.insert(xDivs,xDivs-1)= -ax;
U1x.makeCompressed();
U2x = U1x.conjugate();
U2x.makeCompressed();
for(int i=1;i<yDivs;i++){
    U1y.insert(i,i)= 1.0*J; //just initialise to some random variable
    U1y.insert(i,i+1)= -ay;
    U1y.insert(i,i-1)= -ay;
}
U1y.insert(0,0)= 1.0*J;
U1y.insert(0,1)= -ay;
U1y.insert(yDivs,yDivs)= 1.0*J;
U1y.insert(yDivs,yDivs-1)= -ay;
U1y.makeCompressed();
U2y = U1y.conjugate();
U2y.makeCompressed();
//We initialise the LU solvers
SparseLU<SparseMatrix<cdouble>> LUsolverx;
SparseLU<SparseMatrix<cdouble>> LUsolvery;
//We declare the SPCWFs and the auxiliar vectors for velocity field computation
VectorXcd psiX(xDivs+1), psiY(yDivs+1), U2psix(xDivs+1), U2psiy(yDivs+1);
VectorXcd conjPsix(xDivs+1), conjPsiy(yDivs+1), auxX(xDivs+1), auxY(yDivs+1);
ArrayXd probDensityx(xDivs+1), probDensityy(yDivs+1), velocityFieldx(xDivs+1), velocityFieldy(yDivs+1), auxArrayx(xDivs+1), auxArrayy(yDivs+1);
//we define the trajectory matrix
double** traj=new double*[timeIts+2];
for (int i=0; i<=(timeIts+1); ++i){ traj[i]= new double[4];} // the trajectory is saved in an array of timeIts arrays of 4 doubles (xi, yi, vxi, vyi)
// each of the timeIts arrays contains the value for the trajectory in each of the x,y at that iteration
traj[timeIts+1][2]=0.0; //timeIts+2 traj points are computed, while only timeIts+1 general iterations are made -thus that many speeds we will have
traj[timeIts+1][3]=0.0;

double vx, vy;
//We open the output streams
ofstream probabDataFile, trajDataFile, DATA_chiInfo, DATA_sumChiInfo, DATA_G_J_x, DATA_G_J_y, DATA_KinAdv_x, DATA_KinAdv_y, DATA_XO_Re_Uj_x, DATA_XO_Im_Uj_x;
//psiDataFile.open("DATA_rawSimulationData_nD_XO_ZERO_CN_ABC_tDEP.txt");
probabDataFile.open("DATA_probabilityToPlot_2D_XO_KinAdv_BornHuang_tINDEP.txt");
trajDataFile.open("DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$28$.txt");
DATA_chiInfo.open("DATA_chiInfo_XO.txt");
DATA_sumChiInfo.open("DATA_sumChiInfo_XO.txt");
DATA_G_J_x.open("DATA_G_J_x_KA.txt");
//DATA_G_J_y.open("DATA_G_J_y_KA.txt");
//DATA_XO_Re_Uj_x.open("DATA_XO_Re_Uj.txt");
//DATA_XO_Im_Uj_x.open("DATA_XO_Im_Uj.txt");

DATA_KinAdv_x.open("DATA_KinAdv_x.txt");
//DATA_KinAdv_y.open("DATA_KinAdv_y.txt");

//psiDataFile << std::setprecision(17);
probabDataFile << std::setprecision(17);
trajDataFile << std::setprecision(17);
//BEGINNING OF THE ALGORITHM FOR EACH OF THE INITIAL CONDITIONS----------------------------------------------------------------
for(int trajNum=0; trajNum<numTrajs; ++trajNum){ //this is a potential multithreading branching point
  //We initialise the SPCWF conditioning the full WF to the intial values of the trajectories of this iteration
  for(int i=0; i<=xDivs; ++i){
    psiX(i) = initialFullWF(xgrid(i), initialPosy[trajNum]);
  }
  for(int i=0; i<=yDivs; ++i){
    psiY(i) = initialFullWF(initialPosx[trajNum], ygrid(i));
  }
  traj[0][0]=initialPosx[trajNum];
  traj[0][1]=initialPosy[trajNum];
  // TIME ITERATIONS BEGIN -----------------------------------------------------------------------------------------------------
  for(int it=0; it<=timeIts; ++it){
    //NEXT POSITION OF THE TRAJECTORY OBTAINED -----------------------------
    //Using the current SPCWF, the velocity field of each dimension at this time is obtained and the next position of the trajectory calculated
    //first for the x dimension-------------------------------------
    //we first get the probability density function and the inverse psi
    probDensityx =abs2(psiX.array());
    conjPsix = conj(psiX.array());
    // psi_qi^-1* diff(psi_qi,qi) is computed:
    //the borders are get with an Euler difference o(dq) and the immediate divisions with a central difference o(dq^2)
    auxX(0) = conjPsix(0)*(psiX(1)-psiX(0))/(dx*probDensityx(0));
    auxX(1) = conjPsix(1)*(psiX(2)-psiX(0))/(2.0*dx*probDensityx(1));
    //the rest of points are got with a o(dq^4) difference
    for(int i=2; i<=xDivs-2; ++i){
      auxX(i) = conjPsix(i)*(-psiX(i+2) + 8.0*psiX(i+1) - 8.0*psiX(i-1) +psiX(i-2))/(12*dx*probDensityx(i));
    }
    auxX(xDivs-1) = conjPsix(xDivs-1)*(psiX(xDivs)-psiX(xDivs-2))/(2.0*dx*probDensityx(xDivs-1));
    auxX(xDivs) = conjPsix(xDivs)*(psiX(xDivs)-psiX(xDivs-1))/(dx*probDensityx(xDivs));
    // imaginary part is extracted and the velocity field obtained
    velocityFieldx = (hbar/mx)*imag(auxX.array());
    //now the y dimension------------------------------------
    //we first get the probability density function and the inverse psi
    probDensityy =abs2(psiY.array());
    conjPsiy = conj(psiY.array());
    // psi_qi^-1* diff(psi_qi,qi) is computed:
    //the borders are get with an Euler difference o(dq) and the immediate divisions with a central difference o(dq^2)
    auxY(0) = conjPsiy(0)*(psiY(1)-psiY(0))/(dy*probDensityy(0));
    auxY(1) = conjPsiy(1)*(psiY(2)-psiY(0))/(2.0*dy*probDensityy(1));
    //the rest of points are got with a o(dq^4) difference
    for(int i=2; i<=yDivs-2; ++i){
      auxY(i) = conjPsiy(i)*(-psiY(i+2) + 8.0*psiY(i+1) - 8.0*psiY(i-1) +psiY(i-2))/(12*dy*probDensityy(i));
    }
    auxY(yDivs-1) = conjPsiy(yDivs-1)*(psiY(yDivs)-psiY(yDivs-2))/(2.0*dy*probDensityy(yDivs-1));
    auxY(yDivs) = conjPsiy(yDivs)*(psiY(yDivs)-psiY(yDivs-1))/(dy*probDensityy(yDivs));
    // imaginary part is extracted and the velocity field obtained
    velocityFieldy = (hbar/my)*imag(auxY.array());
    //we apply the discretisation of the grid to the traj positions
    fractional = std::modf((traj[it][0]-xmin)/dx, &wholef);
    whole = wholef;
    if(whole>=xDivs){whole=xDivs-2;}else if(whole<0){whole=0;}
    vx=(1-fractional)*velocityFieldx(whole)+fractional*velocityFieldx(whole+1);
    traj[it+1][0] = traj[it][0]+vx*dt;
    traj[it][2] = vx;
    fractional = std::modf((traj[it][1]-ymin)/dy, &wholef);
    whole = wholef;
    if(whole>=yDivs){whole=yDivs-2;}else if(whole<0){whole=0;}
    vy= (1-fractional)*velocityFieldy(whole)+fractional*velocityFieldy(whole+1);
    traj[it+1][1] = traj[it][1]+vy*dt;
    traj[it][3] = vy;

    //The norms of the SPCWFs Nx and Ny for Uj term calculation are obtained with a composed trapezium rule--------------------------------------------------------------------
    //for Nx
    Nx=0.5*(probDensityx(0)+probDensityx(xDivs));
    for(int i=1; i<xDivs; ++i){Nx+=probDensityx(i);}
    Nx*=dx;
    Ny=0.5*(probDensityy(0)+probDensityy(yDivs));
    for(int i=1; i<yDivs; ++i){Ny+=probDensityy(i);}
    Ny*=dy;
    //Uj(x,t) values of the XO algorithm are calculated so they can be ----------------------------------------------------------------------------
    lastjUsedInItx=-1.0;
    sumaParaChisx=0.0;
    for(int j=0; j<=xjmax; ++j){
      for(int i=0; i<=xDivs; ++i){
        posx = xgrid(i);
        //we get the Uj for this x and this j

        Uj=0.5*(eigenstatesForSectionsInx(ymin,posx,j)*psiY(0) + eigenstatesForSectionsInx(ymax,posx,j)*psiY(yDivs));
        for(int k=1; k<yDivs; ++k){ Uj=Uj+eigenstatesForSectionsInx(ygrid(k),posx,j)*psiY(k);}
        Ujx_container(i, j)=Uj*dy/((cdouble) sqrt(Nx*Ny));
        Chijx_container(i,j)=Ujx_container(i, j)*psiX(i);
      }
      lastjUsedInItx=j;
      sumaParaChisx+=abs2(Chijx_container.col(j)).sum();
      if((sumaChisx(j)=sumaParaChisx*dx)>=chiSumTolerance){break;}
      //Ujx_normFactors(j) = sqrt((abs2(psiX.array()*Ujx_container.col(j)).sum())*(xjmax+1)); //sure we miss the 1/2 for the first and last elements in the trapezium rule, but due to the number of entries in the vector this will turn out negligible
    }
    if(b_y!=0){
    lastjUsedInIty=-1.0;
    sumaParaChisy=0.0;
    for(int j=0; j<=yjmax; ++j){
      for(int i=0; i<=yDivs; ++i){
        posy = ygrid(i);
        //we get the Uj for this y and this j
        Uj=0.5*(eigenstatesForSectionsIny(xmin,posy,j)*psiX(0) + eigenstatesForSectionsIny(xmax,posy,j)*psiX(xDivs));
        for(int k=1; k<xDivs; ++k){ Uj=Uj+eigenstatesForSectionsIny(xgrid(k),posy,j)*psiX(k);}
         Ujy_container(i, j)=Uj*dx/((cdouble) sqrt(Nx*Ny));
        Chijy_container(i,j)=Ujy_container(i, j)*psiY(i);
      }
      lastjUsedInIty=j;
      sumaParaChisy+=abs2(Chijy_container.col(j)).sum();
      if((sumaChisy(j)=sumaParaChisy*dy)>=chiSumTolerance){break;}
      //Ujy_normFactors(j) = sqrt((abs2(psiY.array()*Ujy_container.col(j)).sum())*(yjmax+1));
    }
    }
    //The Ui propagator matrices of each dimension x,y are updated for this time iteration-------------------------------------------------------------------------------------
    posy = traj[it][1];
    for(int i=0; i<=xDivs; ++i){
      posx = xgrid(i);
      kineticCor = 0.0;
      advectiveCor = 0.0;
      //correlPot =0.0;
      for(int j=0; j<=lastjUsedInItx; ++j){ //generate the kinetic and advective correlation potentials for this spatial grid point posx
        kineticCor = kineticCor - Ujx_container(i,j)* 0.5*hbar*hbar*diffyyEigenstatesForSectionsInx(posy, posx, j)/my;
        advectiveCor = advectiveCor + Ujx_container(i,j)* vy*hbar*diffyEigenstatesForSectionsInx(posy, posx, j);

        //correlPot = correlPot + ( Ujx_container(i,j) )*(kineticCor + J*advectiveCor);
      }
      correlPot=$29$*kineticCor+J*advectiveCor*$30$;
      U1x.coeffRef(i,i) = 1.0+J*dt*(hbar*hbar/(mx*dx*dx)+ W(posx, posy) + $31$*correlPot.real()+J*correlPot.imag()*$32$ )/((cdouble)2.0*hbar);
      U2x.coeffRef(i,i) = 1.0-J*dt*(hbar*hbar/(mx*dx*dx)+ W(posx, posy) + $31$*correlPot.real()+J*correlPot.imag()*$32$ )/((cdouble)2.0*hbar);
      G_J_x(i,0) = $31$*correlPot.real();
      G_J_x(i,1) = $32$*correlPot.imag();
      KinAdv_x(i,0)=$29$*kineticCor.real();
      KinAdv_x(i,1)=$29$*kineticCor.imag();
      KinAdv_x(i,2)=$30$*advectiveCor.real();
      KinAdv_x(i,3)=$30$*advectiveCor.imag();
    }
    posx = traj[it][0];
    for(int i=0; i<=yDivs; ++i){
      posy = ygrid(i);
      correlPot =0.0;
      if(b_y!=0){
        for(int j=0; j<=lastjUsedInIty; ++j){ //generate the kinetic and advective correlation potentials for this spatial grid point posx
          kineticCor = -0.5*hbar*hbar*diffxxEigenstatesForSectionsIny(posx, posy, j)/mx;
          advectiveCor = vx*hbar*diffxEigenstatesForSectionsIny(posx, posy, j);

          correlPot = correlPot + ( Ujy_container(i,j) )*(kineticCor + J*advectiveCor);
        }
      }
      U1y.coeffRef(i,i)= 1.0+J*dt*(hbar*hbar/(my*dy*dy)+ W(posx,posy) + b_y*correlPot)/((cdouble)2.0*hbar);
      U2y.coeffRef(i,i)= 1.0-J*dt*(hbar*hbar/(my*dy*dy)+ W(posx,posy) + b_y*correlPot)/((cdouble)2.0*hbar);

      //G_J_y(i,0) = correlPot.real();//nG_J_y(i,1) = correlPot.imag();

    }
    U1x.makeCompressed();
    U1y.makeCompressed();
    U2x.makeCompressed();
    U2y.makeCompressed();
    //LU decomposition done
    LUsolverx.compute(U1x);
    if(LUsolverx.info()!=Success) {
      cout << "LUx decomposition FAILED!" << endl;
      return 1;
    }
    U2psix= U2x*psiX;
    psiX = LUsolverx.solve(U2psix); //the wavefunction of the next time iteration is generated
    LUsolvery.compute(U1y);
    if(LUsolvery.info()!=Success) {
      cout << "LUy decomposition FAILED!" << endl;
      return 1;
    }
    U2psiy= U2y*psiY;
    psiY = LUsolvery.solve(U2psiy);
    if( it%outputDataEvery == 0){ //then we output the data
      probabDataFile <<"KA-Norm_x=" << Nx<<endl<<probDensityx << endl << endl<<endl;
      probabDataFile<<"KA-Norm_y=" << Ny<<endl << probDensityy << endl << endl<<endl;
      for(int j=0; j<=lastjUsedInItx; ++j){
        DATA_sumChiInfo<<j<<" "<<sumaChisx(j)<<endl;
      }
      DATA_sumChiInfo<<endl<<endl;
      DATA_chiInfo <<abs(Chijx_container.leftCols(lastjUsedInItx+1))<<endl<<endl<<endl;
      DATA_G_J_x<<G_J_x<<endl<<endl<<endl;
      //DATA_G_J_y<<G_J_y<<endl<<endl<<endl;
      DATA_KinAdv_x << KinAdv_x<<endl<<endl<<endl;
      //DATA_KinAdv_y << KinAdv_y<<endl<<endl<<endl;
      //DATA_XO_Re_Uj_x << Ujx_container.real() << endl<<endl<<endl;
      //DATA_XO_Im_Uj_x << Ujx_container.imag() << endl << endl << endl;
    }
  } //end time iteration loop
  for(int it=0; it<=timeIts; ++it){
    if( it%outputDataEvery == 0){ //then we output the data
      trajDataFile << traj[it][0] << " " << traj[it][1] << " 0 "<<traj[it][2]<< " "<< traj[it][3]<< endl;
      if(traj[it][0]>=xBound){trajProportionCrossed(it)+=1;}
    }
  }
} //end TrajectoryNumber loop
probabDataFile.close();
trajDataFile.close();
DATA_chiInfo.close();
DATA_sumChiInfo.close();
//DATA_XO_Re_Uj_x.close();
//DATA_XO_Im_Uj_x.close();

//We output the shape of the potential in order to be able to plot it
ofstream potentialToPlot, trajProps;
//we define some output finnes parameters in case it is not necessary to plot the potential to full accuracy (make the algorithm faster)
double potentialPlotFinness=$14$;
int enoughStepx=xDivs*potentialPlotFinness;
int enoughStepy=yDivs*potentialPlotFinness;
double* posArx=new double[xDivs+1];
double* posAry=new double[yDivs+1];
for(int i=0; i<=xDivs; i+=enoughStepx){
  posArx[i]=xgrid(i);
}
for(int j=0; j<=yDivs; j+=enoughStepy){
  posAry[j]=ygrid(j);
}
potentialToPlot.open("DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt");
trajProps.open("DATA_XO_KA_trajProps_k=$28$.txt");
trajProportionCrossed = (1/(double)numTrajs)*trajProportionCrossed;
// we also output the trajectory proportion that crossed x=0 at each time - in order to compute the transmission -
for(int it=0; it<=timeIts; ++it){
    if(it%outputDataEvery ==0){trajProps<<trajProportionCrossed(it)<<endl;
  }
}
for(int i=0; i<=xDivs; i+=enoughStepx){
  for(int j=0; j<=yDivs; j+=enoughStepy){
      potentialToPlot << posArx[i] << " " << posAry[j] << " " << W(posArx[i], posAry[j])<< endl;
  }potentialToPlot<<endl;
}potentialToPlot<<endl;

potentialToPlot.close();
trajProps.close();

for(int i=0; i<=(timeIts+1); i++){ delete[] traj[i]; }
delete[] traj;

return 0;
}


// #################################################################################################
// #################################################################################################

$$2$CODE_simulator_XO_NoGJ.cpp$$


// 2D SCHRODINGER EQUATION SOLVER - XO ALGORITHM NO G, J:
// The Conditional Single Particle Wave Function (CSPWF) of the dimensions (x,y) will be evolved for each initial conditions using a 1D Cranck Nicolson method for the Pseudo Schrodinger Equations
//We include the necessary libraries
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
#define customTrajs $19$

//USER INPUT DECLARATION----------------------------------------------------------------------
//We declare the Spatial and Time grids - names are self-explaning
double xmin = $7$, xmax = $8$, ymin = $9$, ymax = $10$, t0=0.0, dt=$11$, posx, posy, Nx, Ny;
int xDivs=$5$, yDivs = $6$, timeIts=$12$, xBound=$17$, aux;
double dx=(xmax-xmin)/xDivs;
double dy=(ymax-ymin)/yDivs;
int outputDataEvery = $16$; // data will be outputed to the external file only every x time iterations, if 1, every time iteration will be outputted
//The constants: hbar and the masses of each dimension are declared
double hbar=$15$, mx=$3$, my=$4$;
//We declare the initial full Wave Function (x,y,t=t0). It will be used in order to obtain the initial CSPWF by evaluating at each dimension the initial position of the other trajectory
cdouble initialFullWF(double x, double y){
    $1$
}
//We declare the static external Potential Field
double W(double x,double y){
    $2$
}

//We initialize a vector that will number the trajectories that have crossed a certain boundary in each time
ArrayXd trajProportionCrossed = ArrayXd::Zero(timeIts+1);
//The grid positions in x and y are saved into a vector as they are heavily used in the computation of Uj
ArrayXd xgrid(xDivs+1), ygrid(yDivs+1);
int main(){
for(int k=0; k<=xDivs; ++k){xgrid(k)=xmin+k*dx;}
for(int k=0; k<=yDivs;++k){ygrid(k)=ymin+dy*k;}

//The initial positions of each trajectory that will be evolved using the algorithm are chosen according to the probability distribution given by the modulus squared of the initial wave function
int numTrajs=$13$; // we choose the number of trajectories that will be evolved
int gridPointsFullWF=(xDivs+1)*(yDivs+1), whole;
double fractional, wholef;
ArrayXd probabDensity(gridPointsFullWF);
ArrayXcd initialFullPsi(gridPointsFullWF);
double* initialPosx = (double*) malloc(numTrajs*sizeof(double));
double* initialPosy = (double*) malloc(numTrajs*sizeof(double));

if(customTrajs==0){
// the initial state of the full wavefunction is generated in order to obtain its modulus squared in each point
    for(int i=0; i<=xDivs; ++i){
        for(int j=0; j<=yDivs; ++j){
            initialFullPsi(i*(yDivs+1) + j) = initialFullWF(xgrid(i), ygrid(j));
        }
    }
    // the probability associated with each space point is generated
    probabDensity =100*abs2(initialFullPsi);
    // the random number generator initialised
    double* probabClist = probabDensity.data();
    std::default_random_engine generator;
    std::discrete_distribution<int> distribution (probabClist,probabClist+gridPointsFullWF-1);
    // and the initial positions chosen according to the associated probabilities
    for(int i=0; i<numTrajs; i++){
        aux=distribution(generator); //returns the winner index of the prob vector -> we must revert the indexing to 2d indexes
        initialPosx[i] = xmin + ((int) aux/(yDivs+1))*dx;
        initialPosy[i] = ymin + (aux%(yDivs+1))*dy;
    }

}else{ // custom trajectory selection by user. Initial positions should be charged in arrays initialPosx, initialPosy

    $20$

}
// begin the time iterations for each evolved trajectory - UL, UR must be renamed in every iteration - as if it was a time dependant potential algorithm for a 1D particle
//we declare and prepare the propagator matrices for the Cranck Nicolson (CN) evolution
SparseMatrix<cdouble> U1x(xDivs+1, xDivs+1), U2x(xDivs+1, xDivs+1);
SparseMatrix<cdouble> U1y(yDivs+1, yDivs+1), U2y(yDivs+1, yDivs+1);
U1x.reserve(VectorXi::Constant(xDivs+1,3));
U2x.reserve(VectorXi::Constant(xDivs+1,3));
U1y.reserve(VectorXi::Constant(yDivs+1,3));
U2y.reserve(VectorXi::Constant(yDivs+1,3));
cdouble ax=J*dt*hbar/(4.0*mx*dx*dx);
cdouble ay=J*dt*hbar/(4.0*my*dy*dy);
for(int i=1;i<xDivs;i++){
    U1x.insert(i,i)= 1.0*J; //just initialise to some random variable
    U1x.insert(i,i+1)= -ax;
    U1x.insert(i,i-1)= -ax;
}
U1x.insert(0,0)= 1.0*J;
U1x.insert(0,1)= -ax;
U1x.insert(xDivs,xDivs)= 1.0*J;
U1x.insert(xDivs,xDivs-1)= -ax;
U1x.makeCompressed();
U2x = U1x.conjugate();
U2x.makeCompressed();
for(int i=1;i<yDivs;i++){
    U1y.insert(i,i)= 1.0*J; //just initialise to some random variable
    U1y.insert(i,i+1)= -ay;
    U1y.insert(i,i-1)= -ay;
}
U1y.insert(0,0)= 1.0*J;
U1y.insert(0,1)= -ay;
U1y.insert(yDivs,yDivs)= 1.0*J;
U1y.insert(yDivs,yDivs-1)= -ay;
U1y.makeCompressed();
U2y = U1y.conjugate();
U2y.makeCompressed();
//We initialise the LU solvers
SparseLU<SparseMatrix<cdouble>> LUsolverx;
SparseLU<SparseMatrix<cdouble>> LUsolvery;
//We declare the SPCWFs and the auxiliar vectors for velocity field computation
VectorXcd psiX(xDivs+1), psiY(yDivs+1), U2psix(xDivs+1), U2psiy(yDivs+1);
VectorXcd conjPsix(xDivs+1), conjPsiy(yDivs+1), auxX(xDivs+1), auxY(yDivs+1);
ArrayXd probDensityx(xDivs+1), probDensityy(yDivs+1), velocityFieldx(xDivs+1), velocityFieldy(yDivs+1), auxArrayx(xDivs+1), auxArrayy(yDivs+1);
//we define the trajectory matrix
double** traj=new double*[timeIts+2];
for (int i=0; i<=(timeIts+1); ++i){ traj[i]= new double[2];} // the trajectory is saved in an array of timeIts arrays of 2 doubles (xi, yi)
// each of the timeIts arrays contains the value for the trajectory in each of the x,y at that iteration
double vx, vy;
//We open the output streams
ofstream probabDataFile, trajDataFile;
//psiDataFile.open("DATA_rawSimulationData_nD_XO_ZERO_CN_ABC_tDEP.txt");
probabDataFile.open("DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt");
trajDataFile.open("DATA_trajectoriesToPlot_2D_XO_CN_NoGJ_BornHuang_tINDEP_k=$18$.txt");

//psiDataFile << std::setprecision(17);
probabDataFile << std::setprecision(17);
trajDataFile << std::setprecision(17);

//BEGINNING OF THE ALGORITHM FOR EACH OF THE INITIAL CONDITIONS----------------------------------------------------------------
for(int trajNum=0; trajNum<numTrajs; ++trajNum){ //this is a potential multithreading branching point
    //We initialise the SPCWF conditioning the full WF to the intial values of the trajectories of this iteration
    for(int i=0; i<=xDivs; ++i){
        psiX(i) = initialFullWF(xgrid(i), initialPosy[trajNum]);
    }
    for(int i=0; i<=yDivs; ++i){
        psiY(i) = initialFullWF(initialPosx[trajNum], ygrid(i));
    }
    traj[0][0]=initialPosx[trajNum];
    traj[0][1]=initialPosy[trajNum];
    // TIME ITERATIONS BEGIN -----------------------------------------------------------------------------------------------------
    for(int it=0; it<=timeIts; ++it){
    //NEXT POSITION OF THE TRAJECTORY OBTAINED -----------------------------
    //Using the current SPCWF, the velocity field of each dimension at this time is obtained and the next position of the trajectory calculated
    //first for the x dimension-------------------------------------
    //we first get the probability density function and the inverse psi
        probDensityx =abs2(psiX.array());
        conjPsix = conj(psiX.array());
        // psi_qi^-1* diff(psi_qi,qi) is computed:
        //the borders are get with an Euler difference o(dq) and the immediate divisions with a central difference o(dq^2)
        auxX(0) = conjPsix(0)*(psiX(1)-psiX(0))/(dx*probDensityx(0));
        auxX(1) = conjPsix(1)*(psiX(2)-psiX(0))/(2.0*dx*probDensityx(1));
        //the rest of points are got with a o(dq^4) difference
        for(int i=2; i<=xDivs-2; ++i){
            auxX(i) = conjPsix(i)*(-psiX(i+2) + 8.0*psiX(i+1) - 8.0*psiX(i-1) +psiX(i-2))/(12*dx*probDensityx(i));
        }
        auxX(xDivs-1) = conjPsix(xDivs-1)*(psiX(xDivs)-psiX(xDivs-2))/(2.0*dx*probDensityx(xDivs-1));
        auxX(xDivs) = conjPsix(xDivs)*(psiX(xDivs)-psiX(xDivs-1))/(dx*probDensityx(xDivs));
        // imaginary part is extracted and the velocity field obtained
        velocityFieldx = (hbar/mx)*imag(auxX.array());
        //now the y dimension------------------------------------
        //we first get the probability density function and the inverse psi
        probDensityy =abs2(psiY.array());
        conjPsiy = conj(psiY.array());
        // psi_qi^-1* diff(psi_qi,qi) is computed:
        //the borders are get with an Euler difference o(dq) and the immediate divisions with a central difference o(dq^2)
        auxY(0) = conjPsiy(0)*(psiY(1)-psiY(0))/(dy*probDensityy(0));
        auxY(1) = conjPsiy(1)*(psiY(2)-psiY(0))/(2.0*dy*probDensityy(1));
        //the rest of points are got with a o(dq^4) difference
        for(int i=2; i<=yDivs-2; ++i){
            auxY(i) = conjPsiy(i)*(-psiY(i+2) + 8.0*psiY(i+1) - 8.0*psiY(i-1) +psiY(i-2))/(12*dy*probDensityy(i));
        }
        auxY(yDivs-1) = conjPsiy(yDivs-1)*(psiY(yDivs)-psiY(yDivs-2))/(2.0*dy*probDensityy(yDivs-1));
        auxY(yDivs) = conjPsiy(yDivs)*(psiY(yDivs)-psiY(yDivs-1))/(dy*probDensityy(yDivs));
        // imaginary part is extracted and the velocity field obtained
        velocityFieldy = (hbar/my)*imag(auxY.array());
        //we apply the discretisation of the grid to the traj positions
        fractional = std::modf((traj[it][0]-xmin)/dx, &wholef);
        whole = wholef;
        if(whole>=xDivs){whole=xDivs-2;}else if(whole<0){whole=0;}
        vx=(1-fractional)*velocityFieldx(whole)+fractional*velocityFieldx(whole+1);
        traj[it+1][0] = traj[it][0]+vx*dt;
        fractional = std::modf((traj[it][1]-ymin)/dy, &wholef);
        whole = wholef;
        if(whole>=yDivs){whole=yDivs-2;}else if(whole<0){whole=0;}
        vy= (1-fractional)*velocityFieldy(whole)+fractional*velocityFieldy(whole+1);
        traj[it+1][1] = traj[it][1]+vy*dt;
        //The norms of the SPCWFs Nx and Ny for Uj term calculation are obtained with a composed trapezium rule--------------------------------------------------------------------
        //for Nx
        Nx=0.5*(probDensityx(0)+probDensityx(xDivs));
        for(int i=1; i<xDivs; ++i){Nx+=probDensityx(i);}
        Nx*=dx;
        Ny=0.5*(probDensityy(0)+probDensityy(yDivs));
        for(int i=1; i<yDivs; ++i){Ny+=probDensityy(i);}
        Ny*=dy;
        //The Ui propagator matrices of each dimension x,y are updated for this time iteration-------------------------------------------------------------------------------------
        posy = traj[it][1];
        for(int i=0; i<=xDivs; ++i){
            posx = xgrid(i);
            U1x.coeffRef(i,i) = 1.0+J*dt*(hbar*hbar/(mx*dx*dx)+ W(posx, posy))/((cdouble)2.0*hbar);
            U2x.coeffRef(i,i) = 1.0-J*dt*(hbar*hbar/(mx*dx*dx)+ W(posx, posy) )/((cdouble)2.0*hbar);
        }
        posx = traj[it][0];
        for(int i=0; i<=yDivs; ++i){
        posy = ygrid(i);
            U1y.coeffRef(i,i)= 1.0+J*dt*(hbar*hbar/(my*dy*dy)+ W(posx,posy) )/((cdouble)2.0*hbar);
            U2y.coeffRef(i,i)= 1.0-J*dt*(hbar*hbar/(my*dy*dy)+ W(posx,posy))/((cdouble)2.0*hbar);
        }
        U1x.makeCompressed();
        U1y.makeCompressed();
        U2x.makeCompressed();
        U2y.makeCompressed();
        //LU decomposition done
        LUsolverx.compute(U1x);
        if(LUsolverx.info()!=Success) {
            cout << "LUx decomposition FAILED!" << endl;
            return 1;
        }
        U2psix= U2x*psiX;
        psiX = LUsolverx.solve(U2psix); //the wavefunction of the next time iteration is generated
        LUsolvery.compute(U1y);
        if(LUsolvery.info()!=Success) {
            cout << "LUy decomposition FAILED!" << endl;
            return 1;
        }
        U2psiy= U2y*psiY;
        psiY = LUsolvery.solve(U2psiy);

        if( it%outputDataEvery == 0){ //then we output the data
            probabDataFile <<"Norm_x=" << Nx<<endl<<probDensityx << endl << endl<<endl;
            probabDataFile<<"Norm_y=" << Ny<<endl << probDensityy << endl << endl<<endl;
        }

     } // END TIME ITERATIONS
    for(int it=0; it<=timeIts; ++it){
        if( it%outputDataEvery == 0){ //then we output the data
            trajDataFile << traj[it][0] << " " << traj[it][1] << " ";
            trajDataFile <<" 0"<<endl;
            if(traj[it][0]>=xBound){trajProportionCrossed(it)+=1;}
            }
    }
}// END TARJECTORY ITERATIONS
probabDataFile.close();
trajDataFile.close();
//We output the shape of the potential in order to be able to plot it
ofstream potentialToPlot, trajProps;
//we define some output finnes parameters in case it is not necessary to plot the potential to full accuracy (make the algorithm faster)
double potentialPlotFinness=$14$;
int enoughStepx=xDivs*potentialPlotFinness;
int enoughStepy=yDivs*potentialPlotFinness;
double* posArx=new double[xDivs+1];
double* posAry=new double[yDivs+1];
for(int i=0; i<=xDivs; i+=enoughStepx){
posArx[i]=xgrid(i);}
for(int j=0; j<=yDivs; j+=enoughStepy){
posAry[j]=ygrid(j);}
potentialToPlot.open("DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt");
for(int i=0; i<=xDivs; i+=enoughStepx){
    for(int j=0; j<=yDivs; j+=enoughStepy){
        potentialToPlot << posArx[i] << " " << posAry[j] << " " << W(posArx[i], posAry[j])<< endl;
    }
    potentialToPlot<<endl;
}

trajProps.open("DATA_XO_NoGJ_trajProps_k=$18$.txt");
trajProportionCrossed = (1/(double)numTrajs)*trajProportionCrossed;
for(int it=0; it<=timeIts; ++it){
    if(it%outputDataEvery ==0){
        trajProps<<trajProportionCrossed(it)<<endl;
    }
}
potentialToPlot.close();
trajProps.close();

for(int i=0; i<=(timeIts+1); i++){ delete[] traj[i]; }
delete[] traj;
delete[] posArx;
delete[] posAry;
return 0;
}



// ###############################################################################################
// ###############################################################################################

$$3$CODE_simulator_XO_KinAdv.cpp$$

// 2D SCHRODINGER EQUATION SOLVER - XO ALGORITHM Kinetic and Advective Correlation Potential approximation of G and J WITH XABIER'S CORRECTION!
// The Conditional Single Particle Wave Function (CSPWF) of the dimensions (x,y) will be evolved for each initial conditions using a 1D Cranck Nicolson method for the Pseudo Schrodinger Equations
//WARNING! The eigenstates for the sections are not set cdouble because in the box case they are fully real. However, this should be adjusted for a generalized algorithm. As such, complex conjugate should be done in the integrals for the Ujx!!!
//We include the necessary libraries
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
#define customTrajs $33$
//USER INPUT DECLARATION----------------------------------------------------------------------
//We declare the Spatial and Time grids - names are self-explaning
double xmin = $7$, xmax = $8$, ymin = $9$, ymax = $10$, t0=0.0, dt=$11$, posx, posy, Nx, Ny, prob_in_traj_y, prob_in_traj_x;
cdouble Uj, psiY_in_y_traj_pos, psiX_in_x_traj_pos;
int xDivs=$5$, yDivs = $6$, timeIts=$12$, xBound=$27$, aux;
double dx=(xmax-xmin)/xDivs;
double dy=(ymax-ymin)/yDivs;
int outputDataEvery = $16$; // data will be outputed to the external file only every x time iterations, if 1, every time iteration will be outputted
//The constants: hbar and the masses of each dimension are declared
double hbar=$15$, mx=$3$, my=$4$;
//We declare the initial full Wave Function (x,y,t=t0). It will be used in order to obtain the initial CSPWF by evaluating at each dimension the initial position of the other trajectory
cdouble initialFullWF(double x, double y){
$1$
}
//We declare the static external Potential Field
double W(double x,double y){
$2$
}
//The highest used adiabatic states' index in the Born Huang expansion's truncated version
int xjmax = $23$, yjmax = $24$, lastjUsedInItx=0, lastjUsedInIty=0;
//The tolerance for the sum squared of the chis, the algorithm will use in each time iteration as many chis as required to achieve this squared sum tolerance (or rather they will arrive the maximum tolerated j)
double chiSumTolerance=$26$;
double sumaParaChisx, sumaParaChisy;
ArrayXd sumaChisx(xjmax+1), sumaChisy(yjmax+1);
cdouble correlPot, advectiveCor, kineticCor;
//The "analytical" expression for the adiabatic states for each x as a function of y and j are expressed
double eigenstatesForSectionsInx(double y, double x, int j){ //the so called psi_x^j(y)
$17$
}
//The "analytical" expression for the adiabatic states for each y as a function of x and j are expressed
double eigenstatesForSectionsIny(double x, double y, int j){ //the so called psi_y^j(x)
$20$
}
//The analytical expression for the above functions' first and second derivatives
double diffyEigenstatesForSectionsInx(double y, double x, int j){ //the so called d psi_x^j(y)/dy
$18$
}
double diffyyEigenstatesForSectionsInx(double y, double x, int j){ //the so called d**2 psi_x^j(y)/dy**2
$19$
}
double diffxEigenstatesForSectionsIny(double x, double y, int j){ //the so called d psi_y^j(x)/dx
$21$
}
double diffxxEigenstatesForSectionsIny(double x, double y, int j){ //the so called d**2 psi_y^j(x)/dx**2
$22$
}
//We choose if we want the Uj for y dimension to be calculated as well
double b_y=$25$;
//We declare the matrix that will save the Uj(x,t) values for each x and j, in order to allow the calculation of the normalized Uj
ArrayXXcd Ujx_container(xDivs+1, xjmax+1), Ujy_container(yDivs+1, yjmax +1), Chijx_container(xDivs+1, xjmax+1), Chijy_container(yDivs+1, yjmax+1);
//ArrayXd Ujx_normFactors(xjmax+1), Ujy_normFactors(yjmax+1);
//We declare a matrix to store the computed correlation potential G(x)+i*J(x) for a certain time iteration in order to allow its plot. Another one for G(y)+i*J(y)
//We also declare a matrix to store the real and imaginary parts of the Kinetic and Advective correlation potentials approximated for x and for y. (4 columns -> Re{Kin}, Im{Kin}, Re{Adv}, Im{Adv}
ArrayXXd G_J_x(xDivs+1, 2), G_J_y(yDivs+1, 2), KinAdv_x(xDivs+1, 4), KinAdv_y(yDivs+1, 4);

//We initialize a vector that will number the trajectories that have crossed a certain boundary in each time
ArrayXd trajProportionCrossed = ArrayXd::Zero(timeIts+1);
//The grid positions in x and y are saved into a vector as they are heavily used in the computation of Uj
ArrayXd xgrid(xDivs+1), ygrid(yDivs+1);
int main(){
for(int k=0; k<=xDivs; ++k){xgrid(k)=xmin+k*dx;}
for(int k=0; k<=yDivs;++k){ygrid(k)=ymin+dy*k;}


int numTrajs=$13$; // we choose the number of trajectories that will be evolved
int gridPointsFullWF=(xDivs+1)*(yDivs+1), whole;
double fractional, wholef;
double* initialPosx = (double*) malloc(numTrajs*sizeof(double));
double* initialPosy = (double*) malloc(numTrajs*sizeof(double));

if(customTrajs==0){ //The initial positions of each trajectory that will be evolved using the algorithm are chosen according to the probability distribution given by the modulus squared of the initial wave function
    ArrayXcd initialFullPsi(gridPointsFullWF);
    ArrayXd probabDensity(gridPointsFullWF);
    // the initial state of the full wavefunction is generated in order to obtain its modulus squared in each point
    for(int i=0; i<=xDivs; ++i){
      for(int j=0; j<=yDivs; ++j){
        initialFullPsi(i*(yDivs+1) + j) = initialFullWF(xgrid(i), ygrid(j));
      }
    }
    // the probability associated with each space point is generated
    probabDensity =100*abs2(initialFullPsi);
    // the random number generator initialised
    double* probabClist = probabDensity.data();
    std::default_random_engine generator;
    std::discrete_distribution<int> distribution (probabClist,probabClist+gridPointsFullWF-1);
    // and the initial positions chosen according to the associated probabilities
    for(int i=0; i<numTrajs; i++){
      aux=distribution(generator); //returns the winner index of the prob vector -> we must revert the indexing to 2d indexes
      initialPosx[i] = xmin + ((int) aux/(yDivs+1))*dx;
      initialPosy[i] = ymin + (aux%(yDivs+1))*dy;
    }
}else{ // custom trajectory selection by user. Initial positions should be charged in arrays initialPosx, initialPosy

    $34$

}

// begin the time iterations for each evolved trajectory - UL, UR must be renamed in every iteration - as if it was a time dependant potential algorithm for a 1D particle
//we declare and prepare the propagator matrices for the Cranck Nicolson (CN) evolution
SparseMatrix<cdouble> U1x(xDivs+1, xDivs+1), U2x(xDivs+1, xDivs+1);
SparseMatrix<cdouble> U1y(yDivs+1, yDivs+1), U2y(yDivs+1, yDivs+1);
U1x.reserve(VectorXi::Constant(xDivs+1,3));
U2x.reserve(VectorXi::Constant(xDivs+1,3));
U1y.reserve(VectorXi::Constant(yDivs+1,3));
U2y.reserve(VectorXi::Constant(yDivs+1,3));
cdouble ax=J*dt*hbar/(4.0*mx*dx*dx);
cdouble ay=J*dt*hbar/(4.0*my*dy*dy);
for(int i=1;i<xDivs;i++){
    U1x.insert(i,i)= 1.0*J; //just initialise to some random variable
    U1x.insert(i,i+1)= -ax;
    U1x.insert(i,i-1)= -ax;
}
U1x.insert(0,0)= 1.0*J;
U1x.insert(0,1)= -ax;
U1x.insert(xDivs,xDivs)= 1.0*J;
U1x.insert(xDivs,xDivs-1)= -ax;
U1x.makeCompressed();
U2x = U1x.conjugate();
U2x.makeCompressed();
for(int i=1;i<yDivs;i++){
    U1y.insert(i,i)= 1.0*J; //just initialise to some random variable
    U1y.insert(i,i+1)= -ay;
    U1y.insert(i,i-1)= -ay;
}
U1y.insert(0,0)= 1.0*J;
U1y.insert(0,1)= -ay;
U1y.insert(yDivs,yDivs)= 1.0*J;
U1y.insert(yDivs,yDivs-1)= -ay;
U1y.makeCompressed();
U2y = U1y.conjugate();
U2y.makeCompressed();
//We initialise the LU solvers
SparseLU<SparseMatrix<cdouble>> LUsolverx;
SparseLU<SparseMatrix<cdouble>> LUsolvery;
//We declare the SPCWFs and the auxiliar vectors for velocity field computation
VectorXcd psiX(xDivs+1), psiY(yDivs+1), U2psix(xDivs+1), U2psiy(yDivs+1);
VectorXcd conjPsix(xDivs+1), conjPsiy(yDivs+1), auxX(xDivs+1), auxY(yDivs+1);
ArrayXd probDensityx(xDivs+1), probDensityy(yDivs+1), velocityFieldx(xDivs+1), velocityFieldy(yDivs+1), auxArrayx(xDivs+1), auxArrayy(yDivs+1);
//we define the trajectory matrix
double** traj=new double*[timeIts+2];
for (int i=0; i<=(timeIts+1); ++i){ traj[i]= new double[4];} // the trajectory is saved in an array of timeIts arrays of 4 doubles (xi, yi, vxi, vyi)
// each of the timeIts arrays contains the value for the trajectory in each of the x,y at that iteration
traj[timeIts+1][2]=0.0; //timeIts+2 traj points are computed, while only timeIts+1 general iterations are made -thus that many speeds we will have
traj[timeIts+1][3]=0.0;

double vx, vy;
//We open the output streams
ofstream probabDataFile, trajDataFile, DATA_chiInfo, DATA_sumChiInfo, DATA_G_J_x, DATA_G_J_y, DATA_KinAdv_x, DATA_KinAdv_y, DATA_XO_Re_Uj_x, DATA_XO_Im_Uj_x;
//psiDataFile.open("DATA_rawSimulationData_nD_XO_ZERO_CN_ABC_tDEP.txt");
probabDataFile.open("DATA_probabilityToPlot_2D_XO_KinAdv_BornHuang_tINDEP.txt");
trajDataFile.open("DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP_k=$28$.txt");
DATA_chiInfo.open("DATA_chiInfo_XO.txt");
DATA_sumChiInfo.open("DATA_sumChiInfo_XO.txt");
DATA_G_J_x.open("DATA_G_J_x_KA.txt");
//DATA_G_J_y.open("DATA_G_J_y_KA.txt");
//DATA_XO_Re_Uj_x.open("DATA_XO_Re_Uj.txt");
//DATA_XO_Im_Uj_x.open("DATA_XO_Im_Uj.txt");

DATA_KinAdv_x.open("DATA_KinAdv_x.txt");
//DATA_KinAdv_y.open("DATA_KinAdv_y.txt");

//psiDataFile << std::setprecision(17);
probabDataFile << std::setprecision(7);
trajDataFile << std::setprecision(17);
//BEGINNING OF THE ALGORITHM FOR EACH OF THE INITIAL CONDITIONS----------------------------------------------------------------
for(int trajNum=0; trajNum<numTrajs; ++trajNum){ //this is a potential multithreading branching point
  //We initialise the SPCWF conditioning the full WF to the intial values of the trajectories of this iteration
  for(int i=0; i<=xDivs; ++i){
    psiX(i) = initialFullWF(xgrid(i), initialPosy[trajNum]);
  }
  for(int i=0; i<=yDivs; ++i){
    psiY(i) = initialFullWF(initialPosx[trajNum], ygrid(i));
  }
  traj[0][0]=initialPosx[trajNum];
  traj[0][1]=initialPosy[trajNum];
  // TIME ITERATIONS BEGIN -----------------------------------------------------------------------------------------------------
  for(int it=0; it<=timeIts; ++it){
    //NEXT POSITION OF THE TRAJECTORY OBTAINED -----------------------------
    //Using the current SPCWF, the velocity field of each dimension at this time is obtained and the next position of the trajectory calculated
    //first for the x dimension-------------------------------------
    //we first get the probability density function and the inverse psi
    probDensityx =abs2(psiX.array());
    conjPsix = conj(psiX.array());
    // psi_qi^-1* diff(psi_qi,qi) is computed:
    //the borders are get with an Euler difference o(dq) and the immediate divisions with a central difference o(dq^2)
    auxX(0) = conjPsix(0)*(psiX(1)-psiX(0))/(dx*probDensityx(0));
    auxX(1) = conjPsix(1)*(psiX(2)-psiX(0))/(2.0*dx*probDensityx(1));
    //the rest of points are got with a o(dq^4) difference
    for(int i=2; i<=xDivs-2; ++i){
      auxX(i) = conjPsix(i)*(-psiX(i+2) + 8.0*psiX(i+1) - 8.0*psiX(i-1) +psiX(i-2))/(12*dx*probDensityx(i));
    }
    auxX(xDivs-1) = conjPsix(xDivs-1)*(psiX(xDivs)-psiX(xDivs-2))/(2.0*dx*probDensityx(xDivs-1));
    auxX(xDivs) = conjPsix(xDivs)*(psiX(xDivs)-psiX(xDivs-1))/(dx*probDensityx(xDivs));
    // imaginary part is extracted and the velocity field obtained
    velocityFieldx = (hbar/mx)*imag(auxX.array());
    //now the y dimension------------------------------------
    //we first get the probability density function and the inverse psi
    probDensityy =abs2(psiY.array());
    conjPsiy = conj(psiY.array());
    // psi_qi^-1* diff(psi_qi,qi) is computed:
    //the borders are get with an Euler difference o(dq) and the immediate divisions with a central difference o(dq^2)
    auxY(0) = conjPsiy(0)*(psiY(1)-psiY(0))/(dy*probDensityy(0));
    auxY(1) = conjPsiy(1)*(psiY(2)-psiY(0))/(2.0*dy*probDensityy(1));
    //the rest of points are got with a o(dq^4) difference
    for(int i=2; i<=yDivs-2; ++i){
      auxY(i) = conjPsiy(i)*(-psiY(i+2) + 8.0*psiY(i+1) - 8.0*psiY(i-1) +psiY(i-2))/(12*dy*probDensityy(i));
    }
    auxY(yDivs-1) = conjPsiy(yDivs-1)*(psiY(yDivs)-psiY(yDivs-2))/(2.0*dy*probDensityy(yDivs-1));
    auxY(yDivs) = conjPsiy(yDivs)*(psiY(yDivs)-psiY(yDivs-1))/(dy*probDensityy(yDivs));
    // imaginary part is extracted and the velocity field obtained
    velocityFieldy = (hbar/my)*imag(auxY.array());
    //we apply the discretisation of the grid to the traj positions
    fractional = std::modf((traj[it][0]-xmin)/dx, &wholef);
    whole = wholef;
    if(whole>=xDivs){whole=xDivs-2;}else if(whole<0){whole=0;}
    vx=(1-fractional)*velocityFieldx(whole)+fractional*velocityFieldx(whole+1);
    traj[it+1][0] = traj[it][0]+vx*dt;
    traj[it][2] = vx;
    psiX_in_x_traj_pos = (1-fractional)*psiX(whole)+ fractional*psiX(whole+1);
    prob_in_traj_x = probDensityx(whole);

    fractional = std::modf((traj[it][1]-ymin)/dy, &wholef);
    whole = wholef;
    if(whole>=yDivs){whole=yDivs-2;}else if(whole<0){whole=0;}
    vy= (1-fractional)*velocityFieldy(whole)+fractional*velocityFieldy(whole+1);
    traj[it+1][1] = traj[it][1]+vy*dt;
    traj[it][3] = vy;
    psiY_in_y_traj_pos = (1-fractional)*psiY(whole)+ fractional*psiY(whole+1);
    prob_in_traj_y = probDensityy(whole);

    //cout << psiX_in_x_traj_pos << " " << psiY_in_y_traj_pos << endl;

    //The norms of the SPCWFs Nx and Ny for Uj term calculation are obtained with a composed trapezium rule--------------------------------------------------------------------
    //for Nx
    Nx=0.5*(probDensityx(0)+probDensityx(xDivs));
    for(int i=1; i<xDivs; ++i){Nx+=probDensityx(i);}
    Nx*=dx;
    Ny=0.5*(probDensityy(0)+probDensityy(yDivs));
    for(int i=1; i<yDivs; ++i){Ny+=probDensityy(i);}
    Ny*=dy;
    //Uj(x,t) values of the XO algorithm are calculated so they can be ----------------------------------------------------------------------------
    lastjUsedInItx=-1.0;
    sumaParaChisx=0.0;
    for(int j=0; j<=xjmax; ++j){
      for(int i=0; i<=xDivs; ++i){
        posx = xgrid(i);
        //we get the Uj for this x and this j

        Uj=0.5*(eigenstatesForSectionsInx(ymin,posx,j)*psiY(0) + eigenstatesForSectionsInx(ymax,posx,j)*psiY(yDivs));
        for(int k=1; k<yDivs; ++k){ Uj=Uj+eigenstatesForSectionsInx(ygrid(k),posx,j)*psiY(k);}
        Ujx_container(i, j)=Uj*dy/(0.5*(psiY_in_y_traj_pos+psiX_in_x_traj_pos));
        Chijx_container(i,j)=Ujx_container(i, j)*psiX(i);
      }
      lastjUsedInItx=j;
      sumaParaChisx+=abs2(Chijx_container.col(j)).sum();
      if((sumaChisx(j)=sumaParaChisx*dx)>=chiSumTolerance){break;}
      //Ujx_normFactors(j) = sqrt((abs2(psiX.array()*Ujx_container.col(j)).sum())*(xjmax+1)); //sure we miss the 1/2 for the first and last elements in the trapezium rule, but due to the number of entries in the vector this will turn out negligible
    }
    if(b_y!=0){
    lastjUsedInIty=-1.0;
    sumaParaChisy=0.0;
    for(int j=0; j<=yjmax; ++j){
      for(int i=0; i<=yDivs; ++i){
        posy = ygrid(i);
        //we get the Uj for this y and this j
        Uj=0.5*(eigenstatesForSectionsIny(xmin,posy,j)*psiX(0) + eigenstatesForSectionsIny(xmax,posy,j)*psiX(xDivs));
        for(int k=1; k<xDivs; ++k){ Uj=Uj+eigenstatesForSectionsIny(xgrid(k),posy,j)*psiX(k);}
         Ujy_container(i, j)=Uj*dx/(psiX_in_x_traj_pos);
        Chijy_container(i,j)=Ujy_container(i, j)*psiY(i);
      }
      lastjUsedInIty=j;
      sumaParaChisy+=abs2(Chijy_container.col(j)).sum();
      if((sumaChisy(j)=sumaParaChisy*dy)>=chiSumTolerance){break;}
      //Ujy_normFactors(j) = sqrt((abs2(psiY.array()*Ujy_container.col(j)).sum())*(yjmax+1));
    }
    }
    //The Ui propagator matrices of each dimension x,y are updated for this time iteration-------------------------------------------------------------------------------------
    posy = traj[it][1];
    for(int i=0; i<=xDivs; ++i){
      posx = xgrid(i);
      kineticCor = 0.0;
      advectiveCor = 0.0;
      //correlPot =0.0;
      for(int j=0; j<=lastjUsedInItx; ++j){ //generate the kinetic and advective correlation potentials for this spatial grid point posx
        kineticCor = kineticCor - Ujx_container(i,j)* 0.5*hbar*hbar*diffyyEigenstatesForSectionsInx(posy, posx, j)/my;
        advectiveCor = advectiveCor + Ujx_container(i,j)* vy*hbar*diffyEigenstatesForSectionsInx(posy, posx, j);

        //correlPot = correlPot + ( Ujx_container(i,j) )*(kineticCor + J*advectiveCor);
      }
      correlPot=$29$*kineticCor+J*advectiveCor*$30$;
      U1x.coeffRef(i,i) = 1.0+J*dt*(hbar*hbar/(mx*dx*dx)+ W(posx, posy) + $31$*correlPot.real()+J*correlPot.imag()*$32$ )/((cdouble)2.0*hbar);
      U2x.coeffRef(i,i) = 1.0-J*dt*(hbar*hbar/(mx*dx*dx)+ W(posx, posy) + $31$*correlPot.real()+J*correlPot.imag()*$32$ )/((cdouble)2.0*hbar);
      G_J_x(i,0) = $31$*correlPot.real();
      G_J_x(i,1) = $32$*correlPot.imag();
      KinAdv_x(i,0)=$29$*kineticCor.real();
      KinAdv_x(i,1)=$29$*kineticCor.imag();
      KinAdv_x(i,2)=$30$*advectiveCor.real();
      KinAdv_x(i,3)=$30$*advectiveCor.imag();
    }
    posx = traj[it][0];
    for(int i=0; i<=yDivs; ++i){
      posy = ygrid(i);
      correlPot =0.0;
      if(b_y!=0){
        for(int j=0; j<=lastjUsedInIty; ++j){ //generate the kinetic and advective correlation potentials for this spatial grid point posx
          kineticCor = -0.5*hbar*hbar*diffxxEigenstatesForSectionsIny(posx, posy, j)/mx;
          advectiveCor = vx*hbar*diffxEigenstatesForSectionsIny(posx, posy, j);

          correlPot = correlPot + ( Ujy_container(i,j) )*(kineticCor + J*advectiveCor);
        }
      }
      U1y.coeffRef(i,i)= 1.0+J*dt*(hbar*hbar/(my*dy*dy)+ W(posx,posy) + b_y*correlPot)/((cdouble)2.0*hbar);
      U2y.coeffRef(i,i)= 1.0-J*dt*(hbar*hbar/(my*dy*dy)+ W(posx,posy) + b_y*correlPot)/((cdouble)2.0*hbar);
      //G_J_y(i,0) = correlPot.real();//nG_J_y(i,1) = correlPot.imag();

    }
    U1x.makeCompressed();
    U1y.makeCompressed();
    U2x.makeCompressed();
    U2y.makeCompressed();
    //LU decomposition done
    LUsolverx.compute(U1x);
    if(LUsolverx.info()!=Success) {
      cout << "LUx decomposition FAILED!" << endl;
      return 1;
    }
    U2psix= U2x*psiX;
    psiX = LUsolverx.solve(U2psix); //the wavefunction of the next time iteration is generated
    LUsolvery.compute(U1y);
    if(LUsolvery.info()!=Success) {
      cout << "LUy decomposition FAILED!" << endl;
      return 1;
    }
    U2psiy= U2y*psiY;
    psiY = LUsolvery.solve(U2psiy);
    if( it%outputDataEvery == 0){ //then we output the data
      probabDataFile <<"KA-Norm_x=" << Nx<<";|psi(traj)|^2="<<prob_in_traj_x<< endl<<probDensityx << endl << endl<<endl;
      probabDataFile<<"KA-Norm_y=" << Ny<<";|psi(traj)|^2="<<prob_in_traj_y <<endl << probDensityy << endl << endl<<endl;
      for(int j=0; j<=lastjUsedInItx; ++j){
        DATA_sumChiInfo<<j<<" "<<sumaChisx(j)<<endl;
      }
      DATA_sumChiInfo<<endl<<endl;
      DATA_chiInfo <<abs(Chijx_container.leftCols(lastjUsedInItx+1))<<endl<<endl<<endl;
      DATA_G_J_x<<G_J_x<<endl<<endl<<endl;
      //DATA_G_J_y<<G_J_y<<endl<<endl<<endl;
      DATA_KinAdv_x << KinAdv_x<<endl<<endl<<endl;
      //DATA_KinAdv_y << KinAdv_y<<endl<<endl<<endl;
      //DATA_XO_Re_Uj_x << Ujx_container.real() << endl<<endl<<endl;
      //DATA_XO_Im_Uj_x << Ujx_container.imag() << endl << endl << endl;
    }
  } //end time iteration loop
  for(int it=0; it<=timeIts; ++it){
    if( it%outputDataEvery == 0){ //then we output the data
      trajDataFile << traj[it][0] << " " << traj[it][1] << " 0 "<<traj[it][2]<< " "<< traj[it][3]<< endl;
      if(traj[it][0]>=xBound){trajProportionCrossed(it)+=1;}
    }
  }
} //end TrajectoryNumber loop
probabDataFile.close();
trajDataFile.close();
DATA_chiInfo.close();
DATA_sumChiInfo.close();
//DATA_XO_Re_Uj_x.close();
//DATA_XO_Im_Uj_x.close();

//We output the shape of the potential in order to be able to plot it
ofstream potentialToPlot, trajProps;
//we define some output finnes parameters in case it is not necessary to plot the potential to full accuracy (make the algorithm faster)
double potentialPlotFinness=$14$;
int enoughStepx=xDivs*potentialPlotFinness;
int enoughStepy=yDivs*potentialPlotFinness;
double* posArx=new double[xDivs+1];
double* posAry=new double[yDivs+1];
for(int i=0; i<=xDivs; i+=enoughStepx){
  posArx[i]=xgrid(i);
}
for(int j=0; j<=yDivs; j+=enoughStepy){
  posAry[j]=ygrid(j);
}
potentialToPlot.open("DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHuang_tINDEP.txt");
trajProps.open("DATA_XO_KA_trajProps_k=$28$.txt");
trajProportionCrossed = (1/(double)numTrajs)*trajProportionCrossed;
// we also output the trajectory proportion that crossed x=0 at each time - in order to compute the transmission -
for(int it=0; it<=timeIts; ++it){
    if(it%outputDataEvery ==0){trajProps<<trajProportionCrossed(it)<<endl;
  }
}
for(int i=0; i<=xDivs; i+=enoughStepx){
  for(int j=0; j<=yDivs; j+=enoughStepy){
      potentialToPlot << posArx[i] << " " << posAry[j] << " " << W(posArx[i], posAry[j])<< endl;
  }potentialToPlot<<endl;
}potentialToPlot<<endl;

potentialToPlot.close();
trajProps.close();

for(int i=0; i<=(timeIts+1); i++){ delete[] traj[i]; }
delete[] traj;

return 0;
}
