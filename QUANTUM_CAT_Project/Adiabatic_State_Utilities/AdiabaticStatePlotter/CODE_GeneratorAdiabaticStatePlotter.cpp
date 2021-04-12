#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(int argNum, char **argVec){
    if (argNum<9){
        cout << "Error while reading the arguments. Too few arguments introduced? \n";
        return -1;
    }
    //int option;
    //sscanf(argVec[10], "%d", &option);

    ofstream writtenFile;
    writtenFile.open("CODE_AdiabaticStatePlotter.cpp");
    if(!writtenFile.good()){
        cout << "Error while creating the code files!\n";
        return 1;
    }


    writtenFile << "#include <iostream>\n#include <iomanip>\n#include <fstream>\n#include <string>\n#include <complex>\n#include \"LIB_dcomplex.h\" // Macro for using dcomplex as std::complex<double> and J as the complex 0.0+i*1.0\n#include <eigen3/Eigen/Core>\n#include <eigen3/Eigen/Sparse>\n#include <eigen3/Eigen/SparseLU>\n#include <cmath>\n#include <random>\nusing namespace std::complex_literals;\nusing namespace Eigen;\nusing namespace std;\n#define PI 3.141592653589793238463\n#define INF 1000000.0\nint nx="<<argVec[4]<<", ny="<<argVec[5]<<", jmax="<<argVec[3]<<";\ndouble xmin="<<argVec[6]<<", xmax="<<argVec[7]<<", ymin="<<argVec[8]<<", ymax="<<argVec[9]<<", dx, dy;\ncdouble EigenstatesForSectionsInx(double y, double x, int j){ "<<argVec[1]<<"}\ndouble potential(double x, double y){"<<argVec[2]<<"}\nArrayXd posx(nx+1), posy(ny+1);\nint main(){\nofstream plotFile, potentialFile;\nplotFile.open(\"DATA_adiabaticStatePlot.txt\");\n//Prepare the grid in x y\ndx=(xmax-xmin)/nx;\ndy=(ymax-ymin)/ny;\nfor(int ix=0; ix<=nx; ix++){\nposx(ix) = xmin+ix*dx;\n}\nfor(int iy=0; iy<=ny; iy++){\nposy(iy) = ymin+iy*dy;\n}\nfor(int ix=0; ix<=nx; ++ix){\nfor(int iy=0; iy<=ny; ++iy){\nplotFile << posy(iy)<<\" \";\nfor(int j=0; j<=jmax; ++j){\nplotFile << real(EigenstatesForSectionsInx(posy(iy), posx(ix), j))<<\" \";\n}\nplotFile<<endl;\n}\nplotFile<<endl<<endl; }\nplotFile.close();\npotentialFile.open(\"DATA_potentialToPlot.txt\");\nfor(int i=0; i<=nx; i+=1){\nfor(int j=0; j<=ny; j+=1){\npotentialFile << posx(i) << \" \" << posy(j) << \" \" << potential(posx(i), posy(j))<< endl;\n}potentialFile<<endl;}\npotentialFile.close();\nreturn 0;\n}\n";
    
    return 0;
    }
