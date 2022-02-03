/* GAUSS-LEGENDRE IN C++

Created using the F.F.CAMPOS' Algotithm in "Algoritmos Numéricos", LTC 3rd

Author: Davi Ferreira Santiago
    School of Engineering, Universidade Federal de Minas Gerais
Date: February 2nd, 2022 */

#include <iostream>
#include <cmath>
#include <iomanip>

#include "GaussLegendre.hpp"

using namespace std;

int main(){

    double a;
    double b;
    int n;

    cin >> a >> b >> n;

    gauss_legendre(a, b, n);

    return 0;

}