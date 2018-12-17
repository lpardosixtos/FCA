#include "LAVector.hpp"
#include <iostream>
#include <fstream>

//Constructores y destructor

LAVector::LAVector(){
    vect=nullptr;
    sze=0;
}

LAVector::LAVector(int sz1){
    vect=new double[sz1];
    sze=sz1;
}

LAVector::LAVector(int sz1, double val){
    vect=new double[sz1];
    sze=sz1;
    for(int i=0; i<sze; i++) vect[i]=val;
}

LAVector::LAVector(std::string fileName){
    std::ifstream fin;
    fin.open(fileName);
    fin >> sze;
    vect=new double[sze];
    for(int i=0; i<sze; i++) fin >> vect[i];
    fin.close();
}

LAVector::~LAVector(){
    delete[] vect;
}

//-----------------------------------------------------

//Setters
void LAVector::resize(int n){
    if(sze==n) return;
    if(vect!=nullptr){
        delete[] vect;
        sze=n;
        vect=new double[sze];
    }
    else{
        sze=n;
        vect=new double[sze];
    }
}

void LAVector::fill(double val){
    for(int i=0; i<sze; i++) vect[i]=val;
}

void LAVector::readFromPgm(std::string fileName){
    std::ifstream fin;
    fin.open(fileName);
    std::string trash;
    fin >> trash;
    int aux1, aux2;
    fin >> aux1 >> aux2;
    int sz=aux1*aux2;
    this->resize(sz);
    fin >> trash;
    for(int i=0; i<sze; i++){
        fin >> vect[i];
    }
    fin.close();
}

//-------------------------------------------

//Métodos de copiado e impresión

void LAVector::copy(LAVector& v){
    v.resize(sze);
    for(int i=0; i<sze; i++) v[i]=vect[i];
}


void LAVector::print(){
    std::cout << vect[0];
    for(int i=1; i<sze; i++) std::cout << " " << vect[i];
    std::cout << "\n";
}

void LAVector::print(std::string fileName){
    std::ofstream fout;
    fout.open(fileName, std::fstream::trunc);
    fout << sze << "\n";
    for(int i=0; i<sze; i++) fout << vect[i] << "\n";
    fout.close();
}

void LAVector::printPgm(std::string fileName, int rows, int cols){
    std::ofstream fout;
    fout.open(fileName, std::fstream::trunc);
    fout << "P2\n";
    fout << cols << " "<< rows << "\n255\n";
    for(int i=0; i<sze; i++){
        fout << vect[i] << "\n";
    }
    fout.close();
}

//--------------------------------------------------

//Getters
int LAVector::size(){
    return sze;
}

double& LAVector::operator[](int index){
    return vect[index];
}


//------------------------------------------------

//Operaciones de álgebra lineal

double LAVector::dot(LAVector& v1){
    double prod=0;
    for(int i=0; i<sze; i++){
        prod+=vect[i]*v1[i];
    }
    return prod;
}

void LAVector::operator+=(LAVector& v){
    for(int i=0; i<sze; i++) vect[i]+=v[i];
}

void LAVector::operator-=(LAVector& v){
    for(int i=0; i<sze; i++) vect[i]-=v[i];
}

void LAVector::operator*=(double val){
    for(int i=0; i<sze; i++) vect[i]*=val;
}

void LAVector::escalarProd(double val, LAVector& v){
    v.resize(sze);
    for(int i=0; i<sze; i++) v[i]=val*vect[i];
}

void LAVector::sumTo(double val, LAVector& v){
    for(int i=0; i<sze; i++) v[i]+=val*vect[i];
}

double LAVector::norm2(){
    double norm=0;
    for(int i=0; i<sze; i++) norm+=vect[i]*vect[i];
    return sqrt(norm);
}

void LAVector::normalize(){
    double aux=this->norm2();
    for(int i=0; i<sze; i++) vect[i]/=aux;
}


//-----------------------------------------------------
//Static functions

void LAVector::EtoEProd(LAVector& v1, LAVector &v2, LAVector& prod){
    prod.resize(v1.size());
    for(int i=0; i<v1.size(); i++){
        prod[i]=v1[i]*v2[i];
    }
}

void LAVector::sum(LAVector& v1, LAVector& v2, LAVector &suma){
    suma.resize(v1.size());
    for(int i=0; i<v1.size(); i++){
        suma[i]=v1[i]+v2[i];
    }
}

void LAVector::linComb(double a1, LAVector& v1, double a2, LAVector& v2, LAVector &suma){
    suma.resize(v1.size());
    for(int i=0; i<v1.size(); i++){
        suma[i]=a1*v1[i]+a2*v2[i];
    }
}
