$$13$1$CODE_Simulator_CN_Conditional_WFs_with_XO_Trajectories.cpp$$

/* Script to generate the closest possible conditonal wave-functions to the positions
of the trajectories evolved using the XO algorithm.

The idea is to extract from the traejctory output file the computed trajectories using
the XO alg., extract the full wavefunction at each time computed with the 2D CN method
and to slice for each trajectory the full WF in those x and y to get the corresponding
"exact" SPCWF. This also allows the computation of the "exact" U_j(x,t) expected for
that trajectory position y_a(t) for an exact WF.
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <complex>
#include "LIB_dcomplex.h" // Macro for using dcomplex as std::complex<double> and J as the complex 0.0+i*1.0
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>

using namespace std::complex_literals;
using namespace Eigen;
using namespace std;
#define PI 3.141592653589793238463
#define INF 1000000.0


// Declare necessary variables

int xDivs=$3$, yDivs=$4$, jmax=$5$, timeIts=$10$, trajNum=$11$, lastjUsedInItx=0;
double xmin=$6$, ymin=$7$, xmax=$8$, ymax=$9$, chiSumTolerance=$13$, mx=$14$, my=$15$, hbar=$16$, dx, dy, fractionalx, fractionaly;

//The "analytical" expression for the adiabatic states for each x as a function of y and j are expressed
double eigenstatesForSectionsInx(double y, double x, int j){ //the so called psi_x^j(y)
$12$
}

int main(int argNum, char **argVec){

dx = (xmax-xmin)/xDivs;
dy = (ymax-ymin)/yDivs;

// Declare INPUT FILES and necessary variables and arrays ###################################

ifstream DATA_traj_XO_File, DATA_wf_CN_File, DATA_read_Derivatives;

DATA_traj_XO_File.open("$1$"); // de argumento!
DATA_wf_CN_File.open("$2$");

ofstream DATA_conditWF_CN, DATA_CN_Re_Uj, DATA_CN_Im_Uj, DATA_chiInfo, DATA_sumChiInfo, DATA_G_J_Conditional, DATA_Full_WF_G_J, DATA_derivatives;

DATA_Full_WF_G_J.open("DATA_Full_WF_G_J_CN.txt");
DATA_conditWF_CN.open("DATA_probDensity_conditional_WF_CN.txt");
DATA_CN_Re_Uj.open("DATA_CN_Re_Uj.txt");
DATA_CN_Im_Uj.open("DATA_CN_Im_Uj.txt");

DATA_chiInfo.open("DATA_chiInfo_CN_w_KA_trajs.txt");
DATA_sumChiInfo.open("DATA_sumChiInfo_CN_with_KA_traj.txt");
DATA_G_J_Conditional.open("DATA_G_J_with_CN_and_KA_trajs.txt");
DATA_derivatives.open("DATA_derivatives_S_R.txt");

// Declare the arrays to store the trajectories

double** traj = new double*[timeIts+1];
for (int i=0; i<=timeIts; ++i){ traj[i]= new double[4];} //the trajectory will include both position x,y and velocity x,y

double posx;

int gridPoints=(xDivs+1)*(yDivs+1), wholex, wholey;
ArrayXcd WF(gridPoints), psiX(xDivs+1), psiY(yDivs+1), auxComplexVector(gridPoints), conjPsi(gridPoints);
cdouble arrayEl, Uj;

ArrayXd probDensityx(xDivs+1), probDensityy(yDivs+1), abs2WF(gridPoints), absWF(gridPoints), dS_dy(gridPoints), dR_dy(gridPoints), d2S_dy2(gridPoints), d2R_dy2(gridPoints), C(gridPoints), dC_dy(gridPoints), d2C_dy2(gridPoints);
double N_x, N_y, sumaParaChisx, wholexf, wholeyf;

MatrixXd x_y_Gx_Jx(gridPoints, 4), G_J_x_cond(xDivs+1, 2);

ArrayXd xgrid(xDivs+1), ygrid(yDivs+1), sumaChisx(jmax+1);
for(int k=0; k<=xDivs; ++k){xgrid(k)=xmin+k*dx;}
for(int k=0; k<=yDivs;++k){ygrid(k)=ymin+dy*k;}

ArrayXXcd Ujx_container(xDivs+1, jmax+1), Chijx_container(xDivs+1, jmax+1);
string line;

// Start trajectory and time iterations #########################################
// TODO: Due to the order in which traejctories are recorded, there are lots of repeated computations (like the extraction of WF for every j-th iteration of each trajectory or the obtention of dS/dy,...d2R/dy2 for each iteration) - All of them are repeated as many times as trajectories we have
// One could save them just once in enormous arrays of numIt columns of gridPoints rows. However these might be way too big and cause RAM overflow. Thus, a better sterategy might be to save it into output files and the reading them only for doing the conditioning. The other solution would be to compute the stuff for the same iteration of every particle at once (reverse the order of the time iteration and trajectory iteration loops).
// This however would make things complicated when trying to plot the potentials together with the XOKA results.


// For each computed trajectory do
for(int trajIdx=0; trajIdx<trajNum; trajIdx++){

  // Extract the whole trajectory
  for(int tIt=0; tIt<=timeIts; tIt++){
    getline( DATA_traj_XO_File, line);
    istringstream iss(line);
    iss >> traj[tIt][0] >> traj[tIt][1] >> wholexf >> traj[tIt][2] >> traj[tIt][3]; //In the file each line is: x y 0 vx vy
  }

  DATA_wf_CN_File.clear();
  DATA_wf_CN_File.seekg(0);
  getline( DATA_wf_CN_File, line); // First line has simulation parameters

  // For each computed time iteration do
  for(int tIt=0; tIt<=timeIts; tIt++){
      // Extract the full complex WF for this iteration
      getline( DATA_wf_CN_File, line); // the iteration number is written here
      for(int xIt=0; xIt<gridPoints; xIt++){
        getline(DATA_wf_CN_File,line);
        istringstream ss(line);
        ss>>arrayEl;
        WF(xIt)=arrayEl;
      }
      // Condition the WF to the current tIt trajectory position in x and in y
      fractionalx = std::modf((traj[tIt][0]-xmin)/dx, &wholexf);
      fractionaly = std::modf((traj[tIt][1]-ymin)/dy, &wholeyf);

      wholex = wholexf;
      wholey = wholeyf;

      if(wholex>=xDivs){wholex=xDivs-2;}else if(wholex<0){wholex=0;} //make the correction to avoid indexing outside the matrix
      if(wholey>=yDivs){wholey=yDivs-2;}else if(wholey<0){wholey=0;}

      // psiX will be WF(x,t; y=y_a(t)), while pisY is WF(y,t; x= x_a(t))
      for (int yIt=0; yIt<yDivs; yIt++){
         psiY(yIt) = (1-fractionalx)*WF(wholex*(yDivs+1)+ yIt) + fractionalx*WF((wholex+1)*(yDivs+1)+ yIt);}
      for (int xIt=0; xIt<xDivs; xIt++){
         psiX(xIt) = (1-fractionaly)*WF(xIt*(yDivs+1)+ wholey) + fractionaly*WF(xIt*(yDivs+1)+ wholey+1);}

      // Compute the norms of the conditional WF-s
      probDensityx = abs2(psiX);
      probDensityy = abs2(psiY);

      N_x=0.5*(probDensityx(0)+probDensityx(xDivs));
      for(int i=1; i<xDivs; ++i){N_x+=probDensityx(i);}
      N_x*=dx;
      N_y=0.5*(probDensityy(0)+probDensityy(yDivs));
      for(int i=1; i<yDivs; ++i){N_y+=probDensityy(i);}
      N_y*=dy;

      // Compute the U^j(x,t) for this trajectory y_a(t) and this time for each j in (0,...,jmax)

      lastjUsedInItx=-1.0;
      sumaParaChisx=0.0;

      for (int j=0; j<=jmax; ++j){
          for(int i=0; i<=xDivs; ++i){
            posx = xgrid(i);
            //we get the Uj for this x and this j

            Uj=0.5*(eigenstatesForSectionsInx(ymin,posx,j)*psiY(0) + eigenstatesForSectionsInx(ymax,posx,j)*psiY(yDivs));
            for(int k=1; k<yDivs; ++k){ Uj=Uj+eigenstatesForSectionsInx(ygrid(k),posx,j)*psiY(k);}
            Ujx_container(i, j)=Uj*dy/((cdouble) sqrt(N_x*N_y));
            Chijx_container(i,j)=Ujx_container(i, j)*psiX(i);
          }
          lastjUsedInItx=j;
          sumaParaChisx+=abs2(Chijx_container.col(j)).sum();
          if((sumaChisx(j)=sumaParaChisx*dx)>=chiSumTolerance){break;}
      }

      // Output the results
      // the norm and the conditional wf probability densities
      DATA_conditWF_CN <<"CN-Norm_x=" << N_x<<endl<<probDensityx << endl << endl<<endl;
      DATA_conditWF_CN<<"CN-Norm_y=" << N_y<<endl << probDensityy << endl << endl<<endl;

      // output the U_jx
      DATA_CN_Re_Uj << Ujx_container.real() << endl << endl << endl;
      DATA_CN_Im_Uj << Ujx_container.imag() << endl << endl << endl;


      // output chi_j_x partial sums
      for(int j=0; j<=lastjUsedInItx; ++j){
        DATA_sumChiInfo<<j<<" "<<sumaChisx(j)<<endl;
      }
      DATA_sumChiInfo<<endl<<endl;

      // output chi_j_x coefficient modulous
      DATA_chiInfo <<abs(Chijx_container.leftCols(lastjUsedInItx+1))<<endl<<endl<<endl;

      //COMPUTE "EXACT" G,J correlation potentials for this trajectory position and velocity
      /* The strategy is the following:
       - Only compute teh derivatives of S and R once, in the first trajectory iteration and
        output them to a file.
       - For the rest of trajectories simply extract the derivatives form the file
      */
      abs2WF = abs2(WF); //there is a coefficient wise inverse! maybe it works for complex
      absWF = abs(WF);

      if(trajIdx==0){
          //conjPsi = inverse(WF);
          conjPsi = conj(WF);
          C = log(absWF);
          // We first compute the numerical derivatives we will need
          // dS(x,y,t)/dy = hbar/R^2*Im(Psi_conj dPsi(x,y,t)/dy) = hbar*Im(Psi^-1 dPsi(x,y,t)/dy)
          // dR(x,y,t)/dR

      for(int i=0; i<=xDivs; i++){
          for(int j=0; j<=yDivs;j++){
              if(j==0){
                  auxComplexVector(i*(yDivs+1))=conjPsi(i*(yDivs+1))*(WF(i*(yDivs+1)+1)-WF(i*(yDivs+1)))/dy;

                  dR_dy(i*(yDivs+1)) = (absWF(i*(yDivs+1)+1)-absWF(i*(yDivs+1)))/dy;

                  //dC_dy(i*(yDivs+1)) = (C(i*(yDivs+1)+1)-C(i*(yDivs+1)))/dy;
              } else if(j==yDivs){
                  auxComplexVector(i*(yDivs+1)+yDivs)=conjPsi(i*(yDivs+1)+yDivs)*(WF(i*(yDivs+1)+yDivs)-WF(i*(yDivs+1)+yDivs-1))/dy;

                  dR_dy(i*(yDivs+1)+yDivs) = (absWF(i*(yDivs+1)+yDivs)-absWF(i*(yDivs+1)+yDivs-1))/dy;

                  //dC_dy(i*(yDivs+1)+yDivs) = (C(i*(yDivs+1)+yDivs)-C(i*(yDivs+1)+yDivs-1))/dy;
              } else if(j==1 || j==(yDivs-1)){
                  auxComplexVector(i*(yDivs+1)+j)=conjPsi(i*(yDivs+1)+j)*(WF(i*(yDivs+1)+j+1)-WF(i*(yDivs+1)+j-1))/(2.0*dy);

                  dR_dy(i*(yDivs+1)+j) = (absWF(i*(yDivs+1)+j+1)-absWF(i*(yDivs+1)+j-1))/(2.0*dy);

                  //dC_dy(i*(yDivs+1)+j) = (C(i*(yDivs+1)+j+1)-C(i*(yDivs+1)+j-1))/(2.0*dy);

              }else{
                  auxComplexVector(i*(yDivs+1)+j)=conjPsi(i*(yDivs+1)+j)*(-WF(i*(yDivs+1)+j+2)+8.0*WF(i*(yDivs+1)+j+1)-8.0*WF(i*(yDivs+1)+j-1)+WF(i*(yDivs+1)+j-2))/(12.0*dy);

                  dR_dy(i*(yDivs+1)+j) = (-absWF(i*(yDivs+1)+j+2)+8.0*absWF(i*(yDivs+1)+j+1)-8.0*absWF(i*(yDivs+1)+j-1)+absWF(i*(yDivs+1)+j-2))/(12.0*dy);

                  //dC_dy(i*(yDivs+1)+j) = (-C(i*(yDivs+1)+j+2)+8.0*C(i*(yDivs+1)+j+1)-8.0*C(i*(yDivs+1)+j-1)+C(i*(yDivs+1)+j-2))/(12.0*dy);
              }
          }
      }
      // imaginary part is extracted and dS/dy = J/|psi|^2 obtained
      dS_dy = hbar*imag(auxComplexVector);

      //d^2S(x,y,t)/dy^2 = d/dy(dS/dy)
      //d^2R(x,y,t)/dy^2 = d/dy(dR/dy)
      // We will use a fourth order second derivative formula for d^2R(x,y,t)/dy^2 as
      // we have an explicit R(x,y,t). However, for S we do not, so we will compute a
      // first derivative of the fourth order to dS/dy to obtain d2S/dy2
      for(int i=0; i<=xDivs; i++){
          for(int j=0; j<=yDivs;j++){
              if(j==0){
                  d2S_dy2(i*(yDivs+1))=(dS_dy(i*(yDivs+1)+1)-dS_dy(i*(yDivs+1)))/dy;

                  d2R_dy2(i*(yDivs+1)) = (dR_dy(i*(yDivs+1)+1)-dR_dy(i*(yDivs+1)))/dy;

                  //d2C_dy2(i*(yDivs+1)) = (dC_dy(i*(yDivs+1)+1)-dC_dy(i*(yDivs+1)))/dy;
              } else if(j==yDivs){
                  d2S_dy2(i*(yDivs+1)+yDivs)=(dS_dy(i*(yDivs+1)+yDivs)-dS_dy(i*(yDivs+1)+yDivs-1))/dy;

                  d2R_dy2(i*(yDivs+1)+yDivs) = (dR_dy(i*(yDivs+1)+yDivs)-dR_dy(i*(yDivs+1)+yDivs-1))/dy;

                  //d2C_dy2(i*(yDivs+1)+yDivs) = (dC_dy(i*(yDivs+1)+yDivs)-dC_dy(i*(yDivs+1)+yDivs-1))/dy;
              } else if(j==1 || j==(yDivs-1)){
                  d2S_dy2(i*(yDivs+1)+j)=(dS_dy(i*(yDivs+1)+j+1)-dS_dy(i*(yDivs+1)+j-1))/(2.0*dy);

                  d2R_dy2(i*(yDivs+1)+j) = (absWF(i*(yDivs+1)+j+1)-2*absWF(i*(yDivs+1)+j)+absWF(i*(yDivs+1)+j-1))/(dy*dy);

                  //d2C_dy2(i*(yDivs+1)+j) = (C(i*(yDivs+1)+j+1)-2*C(i*(yDivs+1)+j)+C(i*(yDivs+1)+j-1))/(dy*dy);
              }else{
                  d2S_dy2(i*(yDivs+1)+j)=(-dS_dy(i*(yDivs+1)+j+2)+8.0*dS_dy(i*(yDivs+1)+j+1)-8.0*dS_dy(i*(yDivs+1)+j-1)+dS_dy(i*(yDivs+1)+j-2))/(12.0*dy);

                  d2R_dy2(i*(yDivs+1)+j) = (-absWF(i*(yDivs+1)+j+2)+16.0*absWF(i*(yDivs+1)+j+1)-30.0*absWF(i*(yDivs+1)+j)+16.0*absWF(i*(yDivs+1)+j-1)-absWF(i*(yDivs+1)+j-2))/(12.0*dy*dy);

                  //d2C_dy2(i*(yDivs+1)+j) = (-C(i*(yDivs+1)+j+2)+16.0*C(i*(yDivs+1)+j+1)-30.0*C(i*(yDivs+1)+j)+16.0*C(i*(yDivs+1)+j-1)-C(i*(yDivs+1)+j-2))/(12.0*dy*dy);
              }
          }
      }

      // We output the computations to a file
      x_y_Gx_Jx.col(0) = dS_dy;
      x_y_Gx_Jx.col(1) = d2S_dy2;
      x_y_Gx_Jx.col(2) = dR_dy;
      x_y_Gx_Jx.col(3) = d2R_dy2;

      DATA_derivatives << x_y_Gx_Jx<<endl;

      // get it ready for the true output - will only be repeated in the first trajectory
      for(int ix=0; ix<=xDivs; ix++){
          for (int iy=0; iy<=yDivs;iy++){
              x_y_Gx_Jx(ix*(yDivs+1)+iy, 0) = xgrid(ix);
              x_y_Gx_Jx(ix*(yDivs+1)+iy, 1) = ygrid(iy);
          }
      }

  } else{
      for(int k=0; k<gridPoints; ++k){
          getline( DATA_read_Derivatives, line);
          istringstream iss(line);
          iss >> dS_dy(k) >> d2S_dy2(k) >> dR_dy(k) >> d2R_dy2(k);
      }
  }

    // Compute G_x(x,y,t) and J_x(x,y,t) for this trajectory
    x_y_Gx_Jx.col(2) = dS_dy*(dS_dy*(1.0/(2.0*my)) - abs2WF*traj[tIt][3])-hbar*hbar/(2.0*my)*d2R_dy2*absWF*abs2WF;
    x_y_Gx_Jx.col(3) = hbar*dR_dy*absWF*( abs2WF*traj[tIt][3]-2*dS_dy/my )-(hbar/2.0*my)*abs2WF*d2S_dy2;

    /*
    x_y_Gx_Jx.col(2)=dS_dy;
    x_y_Gx_Jx.col(3)=dR_dy;

    x_y_Gx_Jx.col(2) = dS_dy*(dS_dy*(1.0/(2.0*my)) - traj[tIt][3])-hbar*hbar/(2.0*my)*d2R_dy2/absWF;
    x_y_Gx_Jx.col(3) = hbar*dR_dy/absWF*( traj[tIt][3]-dS_dy/my )-(hbar/2.0*my)*d2S_dy2;
    */
      // We get the index j corresponding to the current point in the trajectory y(t). If i and j were saved the other way around, slicing the matrices in a certain j would be way simpler -they would be seqNuential pieces-. As they are not, one would need a previous step to slice it trivially.
      // Fortunately, one can slice an eigen Array using a list of indices: dS_dy(seqN(j_slice, xDivs+1, yDivs+1)) meaning take xDivs+1 elements of dS_dy starting from j_slice and jumping yDivs+1 away per index
      if (fractionaly>0.5){
          wholey++; //wholey already contains the correct index j_slice
      }

      // Conditioning the expressions and using the trajectory velocity we can compute G_x(x,t):=G(x,t;y(t)) and J_x(x,t):=J(x,t;y(t))
      G_J_x_cond.col(0) = x_y_Gx_Jx(seqN(wholey, xDivs+1, yDivs+1), 2);

      G_J_x_cond.col(1) = x_y_Gx_Jx(seqN(wholey, xDivs+1, yDivs+1), 3);

      DATA_G_J_Conditional << G_J_x_cond << endl <<endl <<endl;
      //DATA_Full_WF_G_J << x_y_Gx_Jx << endl << endl << endl;
      for(int ix=0; ix<=xDivs; ix++){
          DATA_Full_WF_G_J <<  x_y_Gx_Jx.block(ix*yDivs+ix, 0, yDivs+1, 4)<< endl<<endl;
      }
      DATA_Full_WF_G_J<< endl;

  } //end timeIts
  if(trajIdx==0){
     // this only in the first trajectory iteration
    DATA_derivatives.close();
    DATA_read_Derivatives.open("DATA_derivatives_S_R.txt");
    }else{
        DATA_read_Derivatives.clear();
        DATA_read_Derivatives.seekg(0);
    }
} //end trajectories

for(int i=0; i<=timeIts; i++){ delete[] traj[i]; }
delete[] traj;

DATA_traj_XO_File.close();
DATA_wf_CN_File.close();

DATA_conditWF_CN.close();
DATA_CN_Re_Uj.close();
DATA_CN_Im_Uj.close();

DATA_sumChiInfo.close();
DATA_chiInfo.close();

DATA_G_J_Conditional.close();
DATA_read_Derivatives.close();
DATA_Full_WF_G_J.close();

return 0;
}
