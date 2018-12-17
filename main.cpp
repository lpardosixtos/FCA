#include <iostream>
#include "FCA.hpp"
using namespace std;

int main(int argsc, char** argsv){
    //Lee argumentos
    string::size_type sz;
    int imNum=stod((string) argsv[1], &sz);

    //Instacía la clase FCA
    FCA images1(argsv, imNum, 3);

    //Separa las matrices
    images1.unmix();

    //Imprime los resultados en dos archivos pgm
    vector<string> names(imNum, string(argsv[2]));
    for(int i=1; i<=imNum; i++){
        names[i-1]+=to_string(i);
        names[i-1]+=".pgm";
    }
    images1.printToPgm(names);
    return 0;
}
