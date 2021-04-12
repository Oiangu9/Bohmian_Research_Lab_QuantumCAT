#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

int main(int argNum, char **argVec){
    char c;
    int maxArgNum;
    ifstream readFile;
    ofstream writtenFile;
    std::ostringstream s;


    s <<  "./CODE_FILE_GENERATORS/" << argVec[2]; // The name of the codeFileGenerator file
    writtenFile.open(s.str());
    s.clear();
    s.str("");


    s << "./RAW_SCRIPTS/" << argVec[1]; // The name of the file to read
    readFile.open(s.str());
    s.clear();
    s.str("");

    /* it is expected that the file beggins with
          $$18$1$CODE_Simulator_blahblah.cpp$$
      or any number variant and name variant
      where this means:
         $$lastArgVecNum$optionNumber$nameForGeneratedScript$$
      for the second and third etc options lastArgVecNum must not be introduced:
         $$optionNum$nameForGeneratedScript$$
    */
    readFile >> c>>c;

    if (c!='$'){ cout << "The file to format is not correctly styled. It should begin with something like $$1$CODE_Simulator_blahblah.cpp$$\n\n"; return 1;}
    do{ readFile>>c; s<<c; }while(readFile.peek()!='$');
    std::istringstream iss(s.str());// Maximum number of arguments that the codeFileGenerator should have (for slicing the argVec)
    iss >>maxArgNum;
    s.clear();
    s.str("");
    readFile >> c;

    /*writtenFile << "#include <iostream>\n#include <fstream>\n#include <string>\n#include <sstream>\n\nusing namespace std;\n\nint main(int argNum, char **argVec){\nif (argNum<"<<maxArgNum+1<<"){\n\tcout << \"Error while reading the arguments. Too few arguments introduced? \\n\";\n\treturn -1;\n}\nint option=1;\nsscanf(argVec["<<maxArgNum+1<<"], \"%d\", &option);\nofstream writtenFile;\n";

     In order to allow multioption and multi argnum in a sam codefile genrator, we have to sacrifice
    this strict argnum check, still we will check that at least two arguments are introduced
    namely, an argument s.t. the codeFileGenerator makes sense and the option number
    */

    writtenFile << "#include <iostream>\n#include <fstream>\n#include <string>\n#include <sstream>\n\nusing namespace std;\n\nint main(int argNum, char **argVec){\nif (argNum<2){\n\tcout << \"Error while reading the arguments. Too few arguments introduced? \\n\";\n\treturn -1;\n}\nint option=1;\nsscanf(argVec[argNum-1], \"%d\", &option);\nofstream writtenFile;\n";


    writtenFile << "\nif(option==";
    do{ readFile>>c; writtenFile<<c; }while(readFile.peek()!='$');
    readFile >> c;
    writtenFile << "){\n\twrittenFile.open(\"";

    do{ readFile>>c; writtenFile<<c; }while(readFile.peek()!='$');
    writtenFile<< "\");\n\twrittenFile << \"";

    readFile >> c >> c; // the last $$

    do{
        if (readFile.peek()=='\n') writtenFile << "\\n";
        if(readFile.peek()==' ') writtenFile <<" ";
        readFile >> c;
        if( readFile.eof() ) break;
        if(c=='"') writtenFile <<"\\\"";
        else if(c=='$'){
          if(readFile.peek()!='$') { // A number between single $ signs indicates an argVec
            writtenFile<< '"'<<"<< argVec[";
            do{ readFile >> c;  writtenFile << c; }while(readFile.peek()!='$');
            writtenFile<<"] <<"<<'"';
            readFile >> c;
          }else{ // A number between double $$ signs indicates a different option for the code file generator
            readFile >> c; //the other $
            writtenFile << "\";\n\n}else if(option==";
            do{ readFile>>c; writtenFile<<c; }while(readFile.peek()!='$');
            readFile >> c;
            writtenFile << "){\n\twrittenFile.open(\"";

            do{ readFile>>c; writtenFile<<c; }while(readFile.peek()!='$');
            writtenFile<< "\");\n\twrittenFile << \"";

            readFile >> c >> c; // the last $$
          }
        }else{ writtenFile << c;}
    }while(1);

    writtenFile << "\";\n}\nwrittenFile.close();\n\nreturn 0;\n}";

    writtenFile.close();
    readFile.close();
    return 0;
}
