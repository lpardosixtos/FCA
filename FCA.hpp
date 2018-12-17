#ifndef FCA_HPP
#define FCA_HPP
#include <string>
#include <vector>
#include "tensor.hpp"

class FCA{
    int sz;//Número de datos
    tensor* Data=nullptr;//Tensor que guarda la imágenes
    tensor* Mean=nullptr;//Tensor que guarda la media de cada imagen
    LAMatrix* T=nullptr;//Factor T de la matriz de mezclado
    LAMatrix* TInv=nullptr;//Inversa de T

    void buildT();//Método para construir la matriz T y TInv
    //Método para calcular el parámetro de descenso de gradiente
    double gamma(LAVector&, LAVector&, LAVector&, LAVector&);
    void gradient(LAVector&, LAVector&, double**);
    void removeProy(LAMatrix&, int, LAVector&);

public:
    FCA(char**, int, int);//Constructor que inicializa el tensor de los archivos pgm
    ~FCA();
    void center();//Centra los datos
    void covariance(LAMatrix&);//Calcula la matriz de covarianza de los datos
    void whitening();//Aplica el proceso de whitening a los datos
    void toNormal(LAMatrix&);//Cálculo de las imágenes originales para el caso de dimensión 2
    void printToPgm(std::vector<std::string>);//Imprime las imágenes en archivos pgm
    void Normalize255();//Normaliza cada matriz para que tome valores entre 0 y 255

    void optimize(LAMatrix&);
    void unmix();
};






#endif
