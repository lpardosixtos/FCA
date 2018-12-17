#ifndef TENSOR_HPP_INCLUDED
#define TENSOR_HPP_INCLUDED

#include "LAMatrix.hpp"
#include "LAVector.hpp"

class tensor{
    LAMatrix* T=nullptr;
    int sz;
    int rows, cols;
public:
//Constructores y destructor
    tensor(){};
    tensor(int, int, int);
    tensor(int);
    ~tensor();

//Setters y métodos para redimensionar el tensor
    void resize(int, int, int);
    void resize(int);

//Getters
    int size();
    int getRows();
    int getCols();
    LAMatrix& operator[](int);

//Métodos de copia, impresión y lectura
    void readFromPgms(char**, int, int);
    void copy(tensor&);
    void print();

//Métodos y operaciones de algebra
    void operator+=(tensor&);
    void prod(LAMatrix&, tensor&);
    void outerProdProd(LAMatrix&);
};



#endif
