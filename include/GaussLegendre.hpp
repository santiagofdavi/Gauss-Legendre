/* GAUSS-LEGENDRE IN C++

Created using the F.F.CAMPOS' Algotithm in "Algoritmos Numéricos", LTC 3rd

Author: Davi Ferreira Santiago
    School of Engineering, Universidade Federal de Minas Gerais
Date: February 2nd, 2022 */

#ifndef GAUSSLEGENDRE_H
#define GAUSSLEGENDRE_H

double avaliacao_funcao(double &x);

void abs_pes(int &n, double W[], double T[], int *Info);

void print(int &i, double &x, double &y, double &w, double &t);

void gauss_legendre(double &a, double &b, int &n);

#endif