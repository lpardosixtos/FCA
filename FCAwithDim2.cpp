#include "FCA.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

//Constructor y destructor
FCA::FCA(char** argsv, int imNum, int start){
    Data=new tensor;
    Data->readFromPgms(argsv, imNum, start);
    sz=Data->size();
}

FCA::~FCA(){
    if(Data) delete Data;
    if(Mean) delete Mean;
    if(T) delete T;
    if(TInv) delete TInv;
}


//-------------------------------------------------------------
//Métodos para FCA
void FCA::center(){
    LAMatrix ones(Data->getRows(), Data->getCols(), 1.0);
    if(Mean) delete Mean;
    Mean=new tensor(sz, Data->getRows(), Data->getCols());
    for(int i=0; i<sz; i++){
        (*Data)[i].prod(ones, (*Mean)[i]);
        (*Mean)[i]*=1.0/double(Data->getCols());
        (*Data)[i]-=(*Mean)[i];
    }
}

void FCA::covariance(LAMatrix& covM){
    covM.setShape(sz, sz);
    LAMatrix auxM;
    LAMatrix auxTrans;
    for(int i=0; i<sz; i++){
        for(int j=i; j<sz; j++){
            (*Data)[j].transpose(auxTrans);
            (*Data)[i].prod(auxTrans, auxM);
            covM[i][j]=1.0/double(Data->getRows())*auxM.trace();
            covM[j][i]=covM[i][j];
        }
    }
}

void FCA::buildT(){
    if(T) delete T;
//Inicializa Matrices y vectores auxiliares
    T=new LAMatrix;
    TInv=new LAMatrix;
    LAMatrix Ut;
    LAMatrix U;
    LAVector sigma;

//Calcula y diagonaliza la matriz de covarianza usando Jacobi
    LAMatrix covM;
    covariance(covM);
    covM.Jacobi(sigma, Ut);

//Guarda la matriz U de la factorización e incializa más matrices auxiliares
    Ut.transpose(U);
    LAMatrix auxM;
    LAMatrix auxM1;
    Ut.copy(auxM);
    U.copy(auxM1);

//Calcula T y TInv
    for(int i=0; i<sigma.size(); i++) {
        double aux;
        aux=1.0/sqrt((sigma)[i]);
        auxM[i]*=aux;
        auxM1[i]*=sqrt(sigma[i]);
    }
    Ut.prod(auxM1, *TInv);
    U.prod(auxM, *T);
}

//Realiza el blanqueado de los datos
void FCA::whitening(){
    buildT();
    Data->outerProdProd(*T);
    Mean->outerProdProd(*T);

}

//Calcula la gamma para descenso de gradiente
double FCA::gamma(LAVector& last, LAVector& now, LAVector& gLast, LAVector& gNow){
    LAVector auxV;
    now.copy(auxV);
    auxV-=last;
    LAVector auxG;
    gLast.copy(auxG);
    auxG-=gNow;
    double aux=auxV.dot(auxG);
    aux/=auxG.norm2();
    return aux;
}

//Calcula el gradiente de dimensión 2
void FCA::gradientDim2(LAVector& Q, LAVector& grad, double t[]){
    double aux1=Q[0]*Q[0]*t[0]+2.0*Q[0]*Q[1]*t[1]+Q[1]*Q[1]*t[2];
    double aux2=Q[1]*Q[1]*t[0]-2.0*Q[1]*Q[0]*t[1]+Q[0]*Q[0]*t[2];

    double der1=2.0*Q[0]*t[0]+2.0*Q[1]*t[1];
    double der2=-2.0*Q[1]*t[1]+2.0*Q[0]*t[2];
    grad[0]=2.0*aux1*der1+2.0*aux2*der2;

    der1=2.0*Q[0]*t[1]+2.0*Q[1]*t[2];
    der2=2.0*Q[1]*t[0]-2.0*Q[0]*t[1];
    grad[1]=2.0*aux1*der1+2.0*aux2*der2;
}

void FCA::gradient(LAVector& q, LAVector& grad, double** t){
    double val=0.0;
    for(int i=0; i<sz; i++){
        grad[i]=0.0;
        for(int j=i+1; j<sz; j++){
            val+=2.0*q[i]*q[j]*t[i][j];
            grad[i]+=2.0*q[j]*t[i][j];
        }
        val+=q[i]*q[i]*t[i][i];
        grad[i]+=2.0*q[i]*t[i][i];
    }
    for(int i=0; i<sz; i++) grad[i]*=2.0*val;
}


//Optimización de dimensión 2
void FCA::optimizeDim2(LAMatrix& maximum){
    maximum.setShape(2,2);
//Calculo de las trazas de los tres productos posibles
    double t[3];
    LAMatrix auxM;
    (*Data)[0].prod((*Data)[0], auxM);
    t[0]=auxM.trace();
    (*Data)[0].prod((*Data)[1], auxM);
    t[1]=auxM.trace();
    (*Data)[1].prod((*Data)[1], auxM);
    t[2]=auxM.trace();

    std::cout << t[0] << " "  << t[1] << " " << t[2] << "\n";

//Calcula las matrices Q y grads que guardarán las matrices actuales y anteriores
    LAVector grads[2];
    LAVector Qs[2];
    for(int i=0; i<2; i++){
        grads[i].resize(sz);
        Qs[i].resize(sz);
        Qs[i].fill(0.0);
    }
    Qs[0][0]=1;
    Qs[1][1]=1;
    gradientDim2(Qs[0], grads[0], t);
    gradientDim2(Qs[1], grads[1], t);
    grads[0].print();
    grads[1].print();
    LAVector auxV;
    for(int i=2; i<100000; i++){
        int index=i&1;
        double gam=gamma(Qs[index], Qs[1-index], grads[index], grads[1-index]);
        grads[1-index].escalarProd(gam, auxV);
        Qs[1-index].copy(Qs[index]);
        Qs[index]+=auxV;
        Qs[index].normalize();
    }

//Guarda el máximo obtenido
    maximum[0][0]=Qs[1][0];
    maximum[0][1]=-Qs[1][1];
    maximum[1][0]=Qs[1][1];
    maximum[1][1]=Qs[1][0];
}




//Calcula la matriz original, recibe la matriz Q de la factorización de A
void FCA::toNormal(LAMatrix& Q){
    LAMatrix Qt;
    Q.transpose(Qt);
    Data->outerProdProd(Qt);
    Mean->outerProdProd(Qt);

    *Data+=*Mean;
}

//Normaliza los datos para valores entre 0 y 255
void FCA::Normalize255(){
    double mini, maxi, l;
    for(int i=0; i<sz; i++){
        mini=1e+10;
        maxi=-1e+10;
        for(int k=0; k<Data->getRows(); k++){
            for(int j=0; j<Data->getCols(); j++){
                mini=std::min((*Data)[i][k][j], mini);
                maxi=std::max((*Data)[i][k][j], maxi);
            }
        }
        l=maxi-mini;
        for(int k=0; k<Data->getRows(); k++){
            for(int j=0; j<Data->getCols(); j++){
                (*Data)[i][k][j]=floor(((*Data)[i][k][j]-mini)*255.0/l);
            }
        }
    }
}

//Imprime en archivos Pgm
void FCA::printToPgm(std::vector<std::string> names){
    for(int i=0; i<sz; i++) (*Data)[i].printToPgm(names[i]);
}

//Separa 2 matrices
void FCA::unmixDim2(){
    center();
    whitening();
    LAMatrix maxim;
    optimizeDim2(maxim);
    maxim.print();
    toNormal(maxim);
    Normalize255();
    for(int i=0; i<Data->getRows(); i++){
         for(int j=0; j<Data->getCols(); j++){
             (*Data)[0][i][j]=255-(*Data)[0][i][j];
             (*Data)[1][i][j]=255-(*Data)[1][i][j];
         }
    }
}

void FCA::removeProy(LAMatrix& Qt, int index, LAVector& q){
    LAVector sum(sz, 0.0);
    LAVector auxV;
    for(int i=0; i<index; i++){
        Qt[i].escalarProd(q.dot(Qt[i]), auxV);
        sum+=auxV;
    }
    q-=sum;
}

void FCA::optimize(LAMatrix& maximum){
    LAMatrix Qt(sz, sz, 0.0);
    LAMatrix auxM;
    double** t=new double*[sz];
    for(int i=0; i<sz; i++) t[i]=new double[sz];
    for(int i=0; i<sz; i++){
        for(int j=i; j<sz; j++){
            (*Data)[i].prod((*Data)[j], auxM);
            t[i][j]=auxM.trace();
            t[j][i]=t[i][j];
        }
    }
    LAVector grads[2];
    LAVector Qs[2];
    for(int i=0; i<2; i++){
        grads[i].resize(sz);
        Qs[i].resize(sz);
    }
    for(int k=0; k<sz; k++){
        for(int i=0; i<2; i++) Qs[i].fill(0.0);
        Qs[0][0]=1;
        Qs[1][1]=1;
        gradient(Qs[0], grads[0], t);
        gradient(Qs[1], grads[1], t);
        grads[0].print();
        grads[1].print();
        std::cout << "grads\n";

        LAVector auxV;
        for(int i=2; i<10000; i++){
            int index=i&1;
            removeProy(Qt, k, Qs[1-index]);
            double gam=gamma(Qs[index], Qs[1-index], grads[index], grads[1-index]);
            grads[1-index].escalarProd(gam, auxV);
            Qs[1-index].copy(Qs[index]);
            Qs[index]+=auxV;
            removeProy(Qt, k, Qs[index]);
            Qs[index].normalize();
        }
        Qs[1].copy(Qt[k]);
    }
    Qt.transpose(maximum);
    for(int i=0; i<sz; i++) delete[] t[i];
    delete[] t;
}

void FCA::unmix(){
    center();
    whitening();
    LAMatrix maxim;
    optimize(maxim);
    maxim.print();
    toNormal(maxim);
    Normalize255();
    for(int i=0; i<Data->getRows(); i++){
         for(int j=0; j<Data->getCols(); j++){
             (*Data)[0][i][j]=255-(*Data)[0][i][j];
             (*Data)[1][i][j]=255-(*Data)[1][i][j];
         }
    }
}
