#ifndef LAMATRIX_HPP_INCLUDED
#define LAMATRIX_HPP_INCLUDED

#include <string>
#include "LAVector.hpp"


class LAMatrix{
protected:
    LAVector** M; //Matriz
    int rows;
    int cols;
    char form='0'; //Guarda si la matriz es diagonal, TSUP, TINF
	bool solveD(LAVector&, LAVector&);//Resuelve el sistema si la matriz es diagonal
	bool solveL(LAVector&, LAVector&);//Resuelve el sistema si la matriz es triangular inferior
	bool solveU(LAVector&, LAVector&);//Resuelve el sistema si la matriz es triangular superior

public:

//Constructores y destructores
    LAMatrix(); //Constructor vacío
    LAMatrix(int, int); //Constructor con dimensiones
    LAMatrix(int, int, char); //Constructor con dimensiones y forma
    LAMatrix(int, int, double); //Constructor con dimensiones y llena de un valor
    LAMatrix(std::string); //Lee de archivo en formato tradicional
    ~LAMatrix();

//Getters
    int getRows();
    int getCols();
    LAVector& operator[](int);

//Setters
    void setShape(int);//Inicializa únicamente el número de filas
    void setShape(int, int);//Da dimensiones
    void setShape(int, int, double);//Da dimensiones y llena de un valor
    void setForm(char);
    void fill(double);//Llena de un valor
    void readFromPgms(int, char**, int, int, int);
    void readFromPgm(std::string);
    void toInt();

//Métodos para imprimir y copiar
    void print();//Imprime en consola
    void print(std::string);//Imprime en archivo
    void copy(LAMatrix&);//Copia a la matriz argumento
    void copyT(LAMatrix&);//Copia la traspuesta a la matriz argumento
    void transpose(LAMatrix&);//Guarda la transpuesta en la matriz argumento
    void printToPgm(std::string);

//Inversa de la Matriz
    void inverse(LAMatrix&); //Guarda en el argumento la traspuesta de la inversa

//Operaciones elementales

    void addOwnRows(int, int, double=1);//Suma un renglón multiplicado por un factor a otro renglón
	void addOwnCols(int, int, double=1);//Suma una columna multiplicada por un factor a otra columna
	void addOwnRows(int, int, LAVector&, double=1);//Suma un renglón multiplicado por un factor a otro renglón, en la matriz extendida
	void swapRows(int, int);//Intercambia dos renglones
	void swapCols(int, int);//Intercambia dos columnas
	void swapRows(int, int, LAVector&);//Intercambia dos renglones de la matriz extendida
	void swapCols(int, int, int*);//Intercambia dos columnas de la matriz extendida, guarda los cambios del vector solución
	void multRow(int, double);//Multiplica un renglón por un factor
	void multCol(int, double);//Multiplica una columna por un factor
    double trace();

//Factorizaciones

    bool factLU(LAMatrix&, LAMatrix&);
    bool factPLU(int*&, LAMatrix&, LAMatrix&);
    bool factQR(LAMatrix&, LAMatrix&);

//Solver

    bool solve(LAVector&, LAVector&);//Resuelve el sistema, internamente depende del tipo
    void conjugGrad(LAVector&, LAVector&);//Método de gradiente conjugado, guarda resultado en el segundo argumento

//Operaciones entre matrices y vectores
//Todos los metodos guardan el resultado en el segundo argumento
    void operator*=(double);//Producto por escalar
    void operator+=(LAMatrix&);//Suma de dos matrices
    void operator-=(LAMatrix&);//Resta de dos matrices
    void prod(double, LAMatrix&);//Producto por escalar
    void prod(LAVector&, LAVector&);//Producto por vector
    void prod(LAMatrix&, LAMatrix&);//Producto por matriz
    void prodLeft(LAMatrix&, LAMatrix&);//Producto por matriz por la izquierda
    void prodT(LAMatrix&, LAMatrix&);//Producto por la traspuesta del primer argumento
    void prodLeftT(LAMatrix&, LAMatrix&);//Producto por la izquierda por la traspuesta del primer argumento

//Methods for eigenvalues and eigenvectors

    double powerMax(int=1<<30);//Método de potencia
    double powerMin(int=1<<30);//Método de potencia inversa
    void Jacobi(LAVector&, LAMatrix&);//Método de Jacobi
    void subspaceMax(LAVector&, LAMatrix&);//Método del subespacio
    void subspaceMin(LAVector&, LAMatrix&);//Método del subespacio para los menores valores


// Métodos estáticos

    static void prod(LAMatrix&, LAVector&, LAVector&); //Producto por vector
    static void prod(LAMatrix&, LAMatrix&, LAMatrix&); //Priducto por matriz

};


#endif
