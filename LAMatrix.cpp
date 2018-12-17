#include "LAMatrix.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

//El método agregados para el examen se encuentra al final.

//Constructors and Destructor

LAMatrix::LAMatrix(){
    M=nullptr;
    rows=0;
    cols=0;
}

LAMatrix::LAMatrix (int r, int c){
    M=new LAVector*[r];
    rows=r;
    cols=c;
    for(int i=0; i<rows; i++) M[i]=new LAVector(cols);
}

LAMatrix::LAMatrix(int r, int c, char f){
    M=new LAVector*[r];
    rows=r;
    cols=c;
    for(int i=0; i<rows; i++) M[i]=new LAVector(cols);
    form=f;
}

LAMatrix::LAMatrix (int r, int c, double val){
    M=new LAVector*[r];
    rows=r;
    cols=c;
    for(int i=0; i<rows; i++){
        M[i]=new LAVector(cols);
        for(int j=0; j<cols; j++) (*(M[i]))[j]=val;
    }
}

LAMatrix::LAMatrix(std::string fileName){
    std::ifstream fin;
    fin.open(fileName);
    fin>> rows >> cols;
    M=new LAVector*[rows];
    for(int i=0; i<rows; i++){
        M[i]=new LAVector(cols);
        for(int j=0; j<cols; j++) fin>>(*(M[i]))[j];
    }
    fin.close();
}



LAMatrix::~LAMatrix(){
    if(M){
        for(int i=0; i<rows; i++) delete M[i];
        delete[] M;
    }
}

//------------------------------------------------------------------------------------

//Calcula la transpuesta de la inversa

void LAMatrix::inverse(LAMatrix& Inv){
    //Factorización PLU
    LAMatrix L, U;
    Inv.setShape(cols, rows);
    int* p;
    this->factPLU(p, L, U);

    //Resuelve los n sistemas
    for(int j=0; j<rows; j++){
        LAVector B(rows, 0.0);
        B[j]=1.0;
        LAVector b(rows);
        for(int i=0; i<rows; i++){
            b[i]=B[p[i]];
        }
        LAVector x(rows);
        L.solve(b, x);
        LAVector aux(rows);
        U.solve(x, Inv[j]);
    }
}

//----------------------------------
//Getters

int LAMatrix::getRows(){ return rows; }
int LAMatrix::getCols(){ return cols; }

void LAMatrix::copy(LAMatrix& Copy){
    Copy.setShape(rows,cols);
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++) Copy[i][j]=(*(M[i]))[j];
    }
}

void LAMatrix::copyT(LAMatrix& Copy){
    Copy.setShape(cols, rows);
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++) Copy[j][i]=(*(M[i]))[j];
    }
}

LAVector& LAMatrix::operator[](int i){
    return *(M[i]);
}

void LAMatrix::print(){
    for(int i=0; i<rows; i++){
        std::cout << (*M[i])[0];
        for(int j=1; j<cols; j++) std::cout << " "  << (*M[i])[j];
        std::cout << "\n";
    }
}

void LAMatrix::print(std::string fileName){
    std::ofstream fout;
    fout.open(fileName, std::fstream::trunc);
    fout << rows << " "  << cols << "\n";
    for(int i=0; i<rows; i++){
        fout << (*M[i])[0];
        for(int j=1; j<cols; j++) fout << " "  << (*M[i])[j];
        fout << "\n";
    }
    fout.close();
}

void LAMatrix::printToPgm(std::string fileName){
    std::ofstream fout;
    fout.open(fileName, std::fstream::trunc);
    fout << "P2\n";
    fout << cols << " "  << rows << "\n";
    fout << 255 << "\n";
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++) fout << (*M[i])[j] << "\n";;
    }
    fout.close();
}

//---------------------------
//Setters
void LAMatrix::setShape(int r, int c){
    if(M!=nullptr){
        if(rows!=r){
            for(int i=0; i<rows; i++) delete M[i];
            delete[] M;
            M=new LAVector*[r];
            for(int i=0; i<r; i++) M[i]=new LAVector(c);
            rows=r;
            cols=c;
            return;
            if(cols!=c){
                for(int i=0; i<rows; i++) delete M[i];
                for(int i=0; i<rows; i++) M[i]=new LAVector(c);
                cols=c;
            }
            return;
        }
    }
    M=new LAVector*[r];
    rows=r;
    cols=c;
    for(int i=0; i<rows; i++) M[i]=new LAVector(cols);
    return;
}

void LAMatrix::setShape(int r){
    if(M!=nullptr){
        if(rows!=r){
            for(int i=0; i<rows; i++) delete M[i];
            delete[] M;
            M=new LAVector*[r];
            rows=r;
            return;
        }
    }
    M=new LAVector*[r];
    rows=r;
    return;
}

void LAMatrix::fill(double val){
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++) (*(M[i]))[j]=val;
    }
}

void LAMatrix::setShape(int r, int c, double val){
    if(M!=nullptr){
        if(rows!=r){
            for(int i=0; i<rows; i++) delete M[i];
            delete[] M;
            M=new LAVector*[r];
            for(int i=0; i<r; i++) M[i]=new LAVector(c);
            rows=r;
            cols=c;
            this->fill(val);
            return;
            if(cols!=c){
                for(int i=0; i<rows; i++) delete M[i];
                for(int i=0; i<rows; i++) M[i]=new LAVector(c);
                cols=c;
            }
            this->fill(val);
            return;
        }
    }
    M=new LAVector*[r];
    rows=r;
    cols=c;
    for(int i=0; i<rows; i++) M[i]=new LAVector(cols);
    this->fill(val);
    return;
}

void LAMatrix::setForm(char c){form=c;}

void LAMatrix::transpose(LAMatrix& T){
    T.setShape(cols, rows);
    for(int i=0; i<cols; i++){
        for(int j=0; j<rows; j++){
            T[i][j]=(*M[j])[i];
        }
    }
}

void LAMatrix::readFromPgms(int imageNum, char** argsv, int start, int r, int c){
    this->setShape(imageNum, r*c);
    for(int i=0; i<imageNum; i++){
        std::ifstream fin;
        fin.open((std::string) argsv[start+i]);
        std::string trash;
        fin >> trash >> trash >> trash >>trash;
        for(int j=0; j<cols; j++) fin >> (*M[i])[j];
        fin.close();
    }
}

void LAMatrix::readFromPgm(std::string fileName){
    std::ifstream fin;
    fin.open(fileName);
    std::string trash;
    fin >>trash;
    int r, c;
    fin >> c >> r;
    fin >> trash;
    this->setShape(r,c);
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++) fin >> (*M[i])[j];
    }
    fin.close();
}

void LAMatrix::toInt(){
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++) (*M[i])[j]=floor((*M[i])[j]);
    }
}

//--------------------------------------------------
//Basic methods

void LAMatrix::operator*=(double val){
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++) (*M[i])[j]*=val;
    }
}

void LAMatrix::operator+=(LAMatrix& M1){
    for(int i=0; i<rows; i++) *M[i]+=M1[i];
}
void LAMatrix::operator-=(LAMatrix& M1){
    for(int i=0; i<rows; i++) *M[i]-=M1[i];
}

void LAMatrix::prod(LAVector& v, LAVector& ans){
    ans.resize(rows);
    for(int i=0; i<rows; i++){
        ans[i]=M[i]->dot(v);
    }
}

void LAMatrix::prod(LAMatrix &B, LAMatrix &C){
    int r=rows, c=B.getCols();
    C.setShape(r,c);
    for(int i=0; i<r; i++){
        //#pragma omp parallel for
        for(int j=0; j<c; j++){
            double sum=0;
            for(int k=0; k<B.getRows(); k++){
                sum+=(*M[i])[k]*B[k][j];
            }
            C[i][j]=sum;
        }
    }
}

void LAMatrix::prodLeft(LAMatrix&B, LAMatrix &C){
    int r=B.getRows(), c=cols;
    C.setShape(r,c);
    for(int i=0; i<r; i++){
        //#pragma omp parallel for
        for(int j=0; j<c; j++){
            double sum=0;
            for(int k=0; k<B.getCols(); k++){
                sum+=(*M[k])[j]*B[i][k];
            }
            C[i][j]=sum;
        }
    }
}

void LAMatrix::prodT(LAMatrix &B, LAMatrix&C){
    int r=rows, c=B.getRows();
    C.setShape(r,c);
    for(int i=0; i<r; i++){
        //#pragma omp parallel for
        for(int j=0; j<c; j++){
            double sum=0;
            for(int k=0; k<B.getCols(); k++){
                sum+=(*M[i])[k]*B[j][k];
            }
            C[i][j]=sum;
        }
    }
}

void LAMatrix::prodLeftT(LAMatrix&B, LAMatrix &C){
    int r=B.getCols(), c=cols;
    C.setShape(r,c);
    for(int i=0; i<r; i++){
        //#pragma omp parallel for
        for(int j=0; j<c; j++){
            double sum=0;
            for(int k=0; k<B.getRows(); k++){
                sum+=(*M[k])[j]*B[k][i];
            }
            C[i][j]=sum;
        }
    }
}

void LAMatrix::prod(double val, LAMatrix& prod){
    prod.setShape(rows, cols);
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++) prod[i][j]=val*(*M[i])[j];
    }
}
//----------------------------------------

//Linear Algebra methods

//Elemental operations

void LAMatrix::addOwnRows(int i, int j, double fact){
	for(int k=0; k<cols; k++) {
		double aux=fact*(*M[i])[k];
		(*M[j])[k]+=aux;
	}
}
void LAMatrix::addOwnCols(int i, int j, double fact){
	for(int k=0; k<rows; k++){
		double aux=fact*(*M[k])[i];
		(*M[k])[j]+=aux;
	}
}

void LAMatrix::addOwnRows(int i, int j, LAVector& B, double fact){
	B[j]+=B[i]*fact;//También se altera el valor del vector objetivo
	for(int k=0; k<cols; k++) {
		double aux=fact*(*M[i])[k];
		(*M[j])[k]+=aux;
	}
}

void LAMatrix::swapRows(int n, int m){
	LAVector* aux=M[n];
	M[n]=M[m];
	M[m]=aux;
}

void LAMatrix::swapCols(int n, int m){
	double aux;
	for(int i=0; i<rows; i++){
		aux=(*M[i])[n];
		(*M[i])[n]=(*M[i])[m];
		(*M[i])[m]=aux;
	}
}

void LAMatrix::swapRows(int n, int m, LAVector& B){
	double aux1=B[n];//También se alteran las posiciones del vector objetivo
	B[n]=B[m];
	B[m]=aux1;
	LAVector* aux=M[n];
	M[n]=M[m];
	M[m]=aux;
}

void LAMatrix::swapCols(int n, int m, int* pos){
	int auxint=pos[n];//Se tiene que llevar el estado actual de las posiciones del vector solución
	pos[n]=pos[m];
	pos[m]=auxint;
	double aux;
	for(int i=0; i<rows; i++){
		aux=(*M[i])[n];
		(*M[i])[n]=(*M[i])[m];
		(*M[i])[m]=aux;
	}
}

void LAMatrix::multRow(int n, double fact){
	for(int i=0; i<cols; i++){
		(*M[n])[i]*=fact;
	}
}
void LAMatrix::multCol(int n, double fact){
	for(int i=0; i<rows; i++){
		(*M[i])[n]*=fact;
	}
}

double LAMatrix::trace(){
    double sum=0.0;
    for(int i=0; i<rows; i++) sum+=(*M[i])[i];
    return sum;
}

//Factorizations

bool LAMatrix::factLU(LAMatrix&L, LAMatrix& U){
	L.setShape(rows, cols);
	U.setShape(rows, cols);
	L.setForm('L');
	U.setForm('U');
	for(int l=0; l<rows; l++){
		double sum=0;
		for(int i=0; i<l; i++){
			for(int k=0; k<i; k++) sum+=L[i][k]*U[k][l];
			U[i][l]=(*M[i])[l]-sum;
			L[i][l]=0.0;
		}
		sum=0;
		for(int j=0; j<l; j++){
			for(int k=0; k<j; k++) sum+=L[l][k]*U[k][j];
			L[l][j]=((*M[l])[j]-sum)/U[j][j];
			U[l][j]=0.0;
		}
		sum=0;
		for(int k=0; k<l; k++) sum+=L[l][k]*U[k][l];
		L[l][l]=1;
		U[l][l]=(*M[l])[l]-sum;
		if(U[l][l]<1e-8) return false;
	}
	return true;
}

bool LAMatrix::factPLU(int*& p, LAMatrix &L, LAMatrix& U){
	L.setShape(rows, cols);
	L.setForm('L');
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++) L[i][j]=0;
	}
	this->copy(U);
	U.setForm('U');
	p=new int[rows];
	for(int i=0; i<rows; i++) p[i]=i;
	for(int l=0; l<rows-1; l++){
		double maxi=fabs(U[l][l]);
		int maxIndex=l;
		for(int i=l+1; i<rows; i++){
			if(fabs(U[i][l])>maxi){
				maxi=fabs(U[i][l]);
				maxIndex=i;
			}
		}
		if(maxIndex!=l){
			U.swapRows(maxIndex,l);
			L.swapRows(maxIndex,l);
			int aux=p[maxIndex];
			p[maxIndex]=p[l];
			p[l]=aux;
		}
		for(int i=l+1; i<rows; i++){
			double factor=-U[i][l]/U[l][l];
			U.addOwnRows(l,i, factor);
			L[i][l]=-factor;
		}
	}
	for(int i=0; i<L.getRows(); i++) L[i][i]++;
	return true;
}

bool LAMatrix::factQR(LAMatrix& Q, LAMatrix& R){
    Q.setShape(cols, rows);
    R.setShape(rows, cols);
    R.setForm('U');
    LAMatrix AT;
    this->copyT(AT);
    for(int j=0; j<cols; j++){
        LAVector sum(rows, 0.0);
        for(int i=0; i<j; i++){
            R[i][j]=AT[j].dot(Q[i]);
            for(int k=0; k<rows; k++) sum[k]+=R[i][j]*Q[i][k];
        }
        sum-=AT[j];
        R[j][j]=sum.norm2();
        if(fabs(R[j][j])<1e-8) return false;
        for(int k=0; k<rows; k++) Q[j][k]=-sum[k]/R[j][j];
    }
    return true;
}


//Solvers

bool LAMatrix::solveD(LAVector& B, LAVector& X){
	bool aux=false;
	for(int i=0; i<rows; i++){
        if(fabs((*M[i])[i])<1e-8) return false;
        X[i]=B[i]/(*M[i])[i];
    }
    return true;
}

bool LAMatrix::solveL(LAVector& B, LAVector& X){
	if(fabs((*M[0])[0])<1e-8) return false;
    X[0]=B[0];
    X[0]/=(*M[0])[0];
    for(int i=1; i<rows; i++){
        if(fabs((*M[i])[i])<1e-8) return false;
        X[i]=B[i];
        double aux=0;
        for(int j=0; j<i; j++){
            aux+=(*M[i])[j]*X[j];
        }
        X[i]-=aux;
        X[i]/=(*M[i])[i];
    }
    return true;
}

bool LAMatrix::solveU(LAVector& B, LAVector& X){
	if(fabs((*M[rows-1])[cols-1])<1e-8) return false;
    X[cols-1]=B[rows-1]/(*M[rows-1])[cols-1];
    for(int i=cols-2; i>=0; i--){
        if(fabs((*M[i])[i])<1e-8) return false;
        X[i]=B[i];
        double aux=0;
        for(int j=cols-1; j>i; j--){
            aux+=(*M[i])[j]*X[j];
        }
        X[i]-=aux;
        X[i]/=(*M[i])[i];
    }
    return true;
}

bool LAMatrix::solve(LAVector& B, LAVector&  X){//Manda a llamar al método correspondiente, no se pueden llamar individualmente
	if(form=='D') return solveD(B, X);
	if(form=='L') return solveL(B, X);
	if(form=='U') return solveU(B, X);
	LAMatrix L, U;
	int* p;
	this->factPLU(p,L,U);
	LAVector x0(L.getCols());
	LAVector b(rows);
	for(int i=0; i<rows; i++){
		b[i]=B[p[i]];
	}
	bool flag=L.solve(b, x0);
	if(!flag) return false;
	flag=U.solve(x0, X);
	return flag;
}


//Methods for eigenvalues and eigenvectors

double LAMatrix::powerMax(int iters){
	LAVector x[2];
    x[0].resize(cols);
    x[1].resize(cols);
    x[0].fill(0.0);
	x[0][0]=1;
	double now=0;
	double last;
    int cont=0;
	do{
        int index=cont&1;
        x[index].normalize();
        this->prod(x[index], x[1-index]);
        last=now;
        now=x[1-index].dot(x[1-index]);
        double aux=x[index].dot(x[1-index]);
        now/=aux;
        cont++;
    }while(fabs(now-last)>1e-10 && cont<iters);
	return now;
}

double LAMatrix::powerMin(int iters){
    LAVector x[2];
    x[0].resize(cols);
    x[1].resize(cols);
    x[0].fill(0.0);
	x[0][0]=1;
    LAVector auxX[2];
    auxX[0].resize(cols);
    auxX[1].resize(cols);
	double now=0;
	double last;
	LAMatrix L, U;
	int* p;
	this->factPLU(p,L,U);
    int cont=0;
	do{
        int index=cont&1;
        x[index].normalize();
		for(int i=0; i<cols; i++) auxX[0][i]=x[index][p[i]];
		L.solve(auxX[0], auxX[1]);
		U.solve(auxX[1], x[1-index]);
        last=now;
        now=x[1-index].dot(x[1-index]);
        double aux=x[index].dot(x[1-index]);
        now/=aux;
        cont++;
    }while(fabs(now-last)>1e-10);
	return 1.0/now;
}

void LAMatrix::Jacobi(LAVector& eigenvals, LAMatrix &Id){
    LAMatrix C;
    this->copy(C);
    eigenvals.resize(rows);
    Id.setShape(rows, cols, 0.0);
    for(int i=0; i<rows; i++) Id[i][i]=1.0;
    while(true){
        double maxi=0;
        int iMax, jMax;
        for(int i=0; i<rows; i++){
            for(int j=i+1; j<cols; j++){
                if(maxi<fabs(C[i][j])){
                    maxi=fabs(C[i][j]);
                    iMax=i;
                    jMax=j;
                }
            }
        }
        if(maxi<1e-5) break;
        double a=C[iMax][iMax], b=C[iMax][jMax], d=C[jMax][jMax];
        double theta=atan2(2*b, a-d)/2;
        double c=cos(theta), s=sin(theta);
        C[iMax][iMax]=a*c*c+2*b*c*s+d*s*s;
        C[jMax][jMax]=a*s*s-2*b*c*s+d*c*c;
        C[iMax][jMax]=b*(c*c-s*s)+(d-a)*s*c;
        C[jMax][iMax]=b*(c*c-s*s)+(d-a)*s*c;
        double a1=Id[iMax][iMax], b1=Id[iMax][jMax], d1=Id[jMax][jMax], b2=Id[jMax][iMax];
        Id[iMax][iMax]=a1*c+b2*s;
        Id[iMax][jMax]=d1*s+b1*c;
        Id[jMax][iMax]=b2*c-a1*s;
        Id[jMax][jMax]=-b1*s+d1*c;

        for(int k=0; k<iMax; k++){
            double p=C[k][iMax], q=C[k][jMax];
            C[k][iMax]=p*c+q*s;
            C[iMax][k]=p*c+q*s;
            C[k][jMax]=q*c-p*s;
            C[jMax][k]=q*c-p*s;
            double r=Id[iMax][k], t=Id[jMax][k];
            Id[iMax][k]=r*c+t*s;
            Id[jMax][k]=t*c-r*s;
        }
        for(int k=iMax+1; k<jMax; k++){
            double p=C[k][iMax], q=C[k][jMax];
            C[k][iMax]=p*c+q*s;
            C[iMax][k]=p*c+q*s;
            C[k][jMax]=q*c-p*s;
            C[jMax][k]=q*c-p*s;
            double r=Id[iMax][k], t=Id[jMax][k];
            Id[iMax][k]=r*c+t*s;
            Id[jMax][k]=t*c-r*s;
        }
        for(int k=jMax+1; k<rows; k++){
            double p=C[k][iMax], q=C[k][jMax];
            C[k][iMax]=p*c+q*s;
            C[iMax][k]=p*c+q*s;
            C[k][jMax]=q*c-p*s;
            C[jMax][k]=q*c-p*s;
            double r=Id[iMax][k], t=Id[jMax][k];
            Id[iMax][k]=r*c+t*s;
            Id[jMax][k]=t*c-r*s;
        }
    }
    for(int i=0; i<rows; i++) eigenvals[i]=C[i][i];
}

void LAMatrix::subspaceMax(LAVector &eigVals, LAMatrix& Phi){
    int cont=0;
    while(true){
        LAMatrix B;
        LAMatrix P;
        this->prodT(Phi, P);
        bool flag=false;
        for(int i=0; i<P.getCols(); i++){
            double aux=P[0][i]/Phi[i][0];
            for(int j=1; j<P.getRows(); j++){
                if(fabs(aux-P[j][i]/Phi[i][j])>1e-3){
                    flag=true;
                    break;
                }
            }
        }
        if(!flag) break;
        Phi.prod(P, B);
        cont++;
        B.Jacobi(eigVals, P);
        LAMatrix Aux;
        Phi.copy(Aux);
        Aux.prodLeft(P, Phi);

        for(int i=0; i<Phi.getRows(); i++){
            LAVector temp;
            for(int j=0; j<i; j++){
                Phi[j].escalarProd((Phi[i].dot(Phi[j])), temp);
                Phi[i]-=temp;
            }
            Phi[i].normalize();
            Phi[i].copy(temp);
            this->prod(temp, Phi[i]);
            for(int j=0; j<i; j++){
                Phi[j].escalarProd((Phi[i].dot(Phi[j])), temp);
                Phi[i]-=temp;
            }
            Phi[i].normalize();
        }
    }
    LAVector auxVec;
    for(int i=0; i<Phi.getRows(); i++){
        this->prod(Phi[i], auxVec);
        eigVals[i]=auxVec[0]/Phi[i][0];
    }
}

void LAMatrix::subspaceMin(LAVector &eigVals, LAMatrix& Phi){
    int cont=0;
    LAMatrix L, U;
    int* p=nullptr;
    this->factPLU(p,L,U);
    LAVector auxX[2];
    auxX[0].resize(cols);
    auxX[1].resize(cols);
    while(true){
        LAMatrix B;
        LAMatrix P;
        this->prodT(Phi, P);
        bool flag=false;
        for(int i=0; i<P.getCols(); i++){
            double aux=P[0][i]/Phi[i][0];
            for(int j=1; j<P.getRows(); j++){
                if(fabs(aux-P[j][i]/Phi[i][j])>1e-5){
                    flag=true;
                    break;
                }
            }
        }
        if(!flag) break;
        Phi.prod(P, B);

        cont++;
        B.Jacobi(eigVals, P);
        LAMatrix Aux;
        Phi.copy(Aux);
        Aux.prodLeft(P, Phi);
        for(int i=0; i<Phi.getRows(); i++){
            LAVector temp;
            for(int j=0; j<i; j++){
                Phi[j].escalarProd((Phi[i].dot(Phi[j])), temp);
                Phi[i]-=temp;
            }
            Phi[i].normalize();
            Phi[i].copy(temp);
            for(int j=0; j<cols; j++) auxX[0][j]=Phi[i][p[j]];
    		L.solve(auxX[0], auxX[1]);
    		U.solve(auxX[1], Phi[i]);
            for(int j=0; j<i; j++){
                Phi[j].escalarProd((Phi[i].dot(Phi[j])), temp);
                Phi[i]-=temp;
            }
            Phi[i].normalize();
        }
    }
    LAVector auxVec;
    for(int i=0; i<Phi.getRows(); i++){
        this->prod(Phi[i], auxVec);
        eigVals[i]=auxVec[0]/Phi[i][0];
    }
}

void LAMatrix::conjugGrad(LAVector& b, LAVector& x){
    LAVector r;
    this->prod(x, r);
    r-=b;
    r.escalarProd(-1.0, r);
    LAVector p;
    r.copy(p);
    double aux;
    while(true){
        //Cálculo de alpha
        double alpha=p.dot(r);
        LAVector auxVec;
        this->prod(p, auxVec);
        alpha/=p.dot(auxVec);
        //Calculo de la siguiente x
        p.escalarProd(alpha, auxVec);
        x+=auxVec;
        //Calculo de la siguiente r
        this->prod(p, auxVec);
        auxVec.escalarProd(alpha, auxVec);
        r-=auxVec;
        if(r.norm2()<1e-7) break;
        //Calculo de la siguiete p
        aux=p.dot(r)/p.dot(p);
        p.escalarProd(aux,p);
        p+=r;
    }
}

//----------------------------------------------
//Static methods


void LAMatrix::prod(LAMatrix& A, LAVector& v, LAVector& ans){
    ans.resize(A.getRows());
    for(int i=0; i<A.getRows(); i++){
        ans[i]=A[i].dot(v);
    }
}

void LAMatrix::prod(LAMatrix& A, LAMatrix &B, LAMatrix &C){
    int r=A.getRows(), c=B.getCols();
    C.setShape(r,c);
    for(int i=0; i<r; i++){
        //#pragma omp parallel for
        for(int j=0; j<c; j++){
            double sum=0;
            for(int k=0; k<B.getRows(); k++){
                sum+=A[i][k]*B[k][j];
            }
            C[i][j]=sum;
        }
    }
}
