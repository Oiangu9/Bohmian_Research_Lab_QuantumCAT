#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(int argNum, char **argVec){
    if (argNum<10){
        cout << "Error while reading the arguments. Too few arguments introduced? \n";
        return -1;
    }
    ofstream writtenFile;
        
    writtenFile.open("CODE_2D_WF_EnergyCalculator.cpp");
    if(!writtenFile.good()){
        cout << "Error while creating the code files!\n";
        return 1;
    }
    
    writtenFile << "// 2D WF Energy Caluculator by using the Sandwich method- obv TIME INDEPENDENT POTENTIALS\n#include <iostream>\n#include <iomanip>\n#include <fstream>\n#include <string>\n#include <complex>\n#include \"LIB_dcomplex.h\" // Macro for using dcomplex as std::complex<double> and J as the complex 0.0+i*1.0\n#include <eigen3/Eigen/Core>\n#include <eigen3/Eigen/Sparse>\n#include <eigen3/Eigen/SparseLU>\n#include <cmath>\nusing namespace std::complex_literals;\nusing namespace Eigen;\nusing namespace std;\nconst double PI = 3.141592653589793238463; const double EULER=2.718281828459045;\n//defining the particle parameters\ndouble m1=" << argVec[9] << ", m2=" << argVec[10] << ", hbar=" << argVec[11] << ";\n//define the SPACE GRID\ndouble x1min=" << argVec[5] << ", x1max=" << argVec[6] << ", x2min=" << argVec[7] << ", x2max=" << argVec[8] << ", dt=" <<argVec[12]<<";\ncdouble Energy=0.0;\nint nx1=" << argVec[3] << ", nx2=" << argVec[4] << ";\n//We define the initial WF\ncdouble psi0(double x, double y){\n " << argVec[1] << ";\n}\n//We define the Potential energy field TIME INDEPENDENT\ndouble V(double x, double y){\n " << argVec[2] << ";\n}\nint main(){\n\ndouble dx1=(x1max-x1min)/nx1, dx2=(x2max-x2min)/nx2;\nint gridPoints=(nx1+1)*(nx2+1);\n\n//Build the planed WavePacket (which should be seen as a 2d array\nVectorXcd psi_ini(gridPoints); //nx, ny intervalos, then nx+1, ny+1 ptos por cada dim\nfor(int ix=0; ix<=nx1; ix++){\nfor (int iy=0; iy<=nx2;iy++){\npsi_ini(ix*(nx2+1)+iy)=psi0(x1min+ix*dx1, x2min+iy*dx2);\n}\n}\n//Build the propagators U1 and U2\nSparseMatrix<cdouble> UL(gridPoints, gridPoints);\nUL.reserve(VectorXi::Constant(gridPoints,5));\ncdouble a=J*dt*hbar/(4*m1*dx1*dx1), d=J*dt*hbar/(4*m2*dx2*dx2);\nint j;\nfor(int ix=0;ix<=nx1;ix++){\nfor(int iy=0; iy<=nx2;iy++){\nj=ix*(nx2+1)+iy;\nUL.insert(j,j)=1.0+J*dt*(hbar*hbar*(1/(m1*dx1*dx1)+1/(m2*dx2*dx2))+V(x1min+ix*dx1,x2min+iy*dx2))/((cdouble)2*hbar);\nif(iy!=nx2) UL.insert(j,j+1)=-d;\nif(iy!=0) UL.insert(j,j-1)=-d;\nif(ix>0) UL.insert(j,j-nx2-1)=-a;\nif(ix<nx1) UL.insert(j,j+nx2+1)=-a;\n}\n}\nUL.makeCompressed();\n// UR is exactly de complex conjugate of UL\nSparseMatrix<cdouble> UR(gridPoints, gridPoints);\nUR.reserve(VectorXi::Constant(gridPoints,5));\nUR = UL.conjugate();\nUR.makeCompressed();\n// Make the LU decomposition of UL^transposed\nSparseLU<SparseMatrix<cdouble>> LUsolver;\nLUsolver.compute(UL.transpose());\nif(LUsolver.info()!=Success) {\ncout << \"LU decomposition FAILED!\" << endl;\nreturn 1;\n}\n//We are seeking for E = psi^t UL^-1 UR psi\n//we get first conjpsi^t * UL^-1=v, by v^t=UL^t^-1 * conjpsi <-> UL^t v^t = conjpsi <-> if lu= LU(UL^t) then v^t=lu.solve(conjpsi)\n \nVectorXcd vt(gridPoints),URpsi(gridPoints) ; \nvt=LUsolver.solve(psi_ini.conjugate());\n//cout<<psi_ini<<endl<<endl<<psi_ini.conjugate();\nURpsi=UR*psi_ini;\nEnergy = vt.dot(URpsi);\n//cout<<vt<<endl<<endl<<URpsi;\nofstream outsideFile;\noutsideFile.open(\"DATA_ObtainedEnergy.txt\");\noutsideFile << std::setprecision(17);\noutsideFile << Energy<<endl; cout<<endl<<endl <<\"Calculated Energy = \"<<Energy<<endl<<endl;\noutsideFile.close();\nreturn 0;\n}\n";
    
    
    
    writtenFile.close();
    return 0;
}
