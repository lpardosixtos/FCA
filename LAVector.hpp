#ifndef LAVECTOR_HPP_INCLUDED
#define LAVECTOR_HPP_INCLUDED

#include <omp.h>
#include <cmath>
#include <string>

class LAVector{
    double* vect=nullptr;
    int sze=0;

public:

//Constructores y destructores

    LAVector();//Constructor vacío
    LAVector(int);//Inicializa el tamaño
    LAVector(int, double);//Inicializa con tamaño y llena de un valor
    LAVector(std::string);//Lee vector de archivo
    ~LAVector();

//Métodos de copia e impresión

    void copy(LAVector&); //Crea copia en el argumento
    void print();//Imprime en consola
    void print(std::string);//Imprime en archivo
    void printPgm(std::string, int, int);//Imprime en archivo pgm

//Setters

    void fill(double);//Llena de un valor
    void resize(int);
    void readFromPgm(std::string);//Lee el vector de archivo pgm

//Getters
    double& operator[](int);
    int size();

//operaciones entre vectores

    double dot(LAVector&);
    void operator+=(LAVector&);
    void operator-=(LAVector&);
    void operator*=(double);
    void escalarProd(double, LAVector&);//Guarda el producto por escalar en el segundo argumento
    void sumTo(double, LAVector&);//Suma el vector actual multiplicado por el escalar al segundo alrgumento

//Normas
    double norm2();
    void normalize();//Normaliza el vector

//Métodos estáticos
    static void EtoEProd(LAVector&, LAVector&, LAVector&);
    static void sum(LAVector&, LAVector&, LAVector&);
    static void linComb(double, LAVector&, double, LAVector&, LAVector&);
};
#endif
