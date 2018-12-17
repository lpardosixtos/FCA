#include "tensor.hpp"
#include <iostream>
#include <cmath>

tensor::tensor(int s, int r, int c){
    sz=s;
    rows=r;
    cols=c;
    T=new LAMatrix[sz];
    for(int i=0; i<sz; i++) T[i].setShape(r, c);
}

tensor::tensor(int s){
    sz=s;
    T=new LAMatrix[sz];
}

tensor::~tensor(){
    if(T) delete[] T;
}


//-----------------------------------------------------------------
//Getters

LAMatrix& tensor::operator[](int index){
    return T[index];
}

int tensor::size(){
    return sz;
}

int tensor::getRows(){
    return rows;
}

int tensor::getCols(){
    return cols;
}

//----------------------------------------------------------------
//Setters

void tensor::resize(int s){
    if(T) delete T;
    sz=s;
    T=new LAMatrix[s];
}

void tensor::resize(int s, int r, int c){
    if(T) delete T;
    sz=s;
    rows=r;
    cols=c;
    T=new LAMatrix[s];
    for(int i=0; i<s; i++) T[i].setShape(r,c);
}

void tensor::readFromPgms(char** names, int s, int start){
    this->resize(s);
    for(int i=0; i<sz; i++) T[i].readFromPgm((std::string) names[i+start]);
    rows=T[0].getRows();
    cols=T[0].getCols();
}

//----------------------------------------------------------------
//Método de copia e impresion

void tensor::copy(tensor& C){
    C.resize(sz);
    for(int i=0; i<sz; i++) (*this)[i].copy(C[i]);
}

void tensor::print(){
    for(int i=0; i<sz; i++){
        T[i].print();
        std::cout << "\n";
    }
}

//-----------------------------------------------------------------
//Operaciones de álgebra

void tensor::prod(LAMatrix& Mat, tensor& prod){
    prod.resize(sz);
    LAMatrix auxM;
    LAMatrix sumM(rows, cols);
    for(int k=0; k<Mat.getRows(); k++){
        sumM.fill(0.0);
        for(int i=0; i<size(); i++){
            (*this)[i].prod(Mat[k][i], auxM);
            sumM+=auxM;
        }
        sumM.copy(prod[k]);
    }
}

void tensor::operator+=(tensor& T1){
    for(int i=0; i<sz; i++) (*this)[i]+=T1[i];
}

//Resultado de hacer el producto de Kronecker de Mat por la identidad de tamaño rows y
//luego multiplicar por el tensor
void tensor::outerProdProd(LAMatrix& Mat){
    tensor auxT(sz);
    LAMatrix auxM;
    for(int i=0; i<Mat.getRows(); i++){
        auxT[i].setShape(rows, cols, 0.0);
        for(int j=0; j<sz; j++){
            (*this)[j].prod(Mat[i][j], auxM);
            auxT[i]+=auxM;
        }
    }
    for(int i=0; i<sz; i++) auxT[i].copy((*this)[i]);
}
