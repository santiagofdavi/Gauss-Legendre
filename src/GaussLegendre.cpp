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

double avaliacao_funcao(double &x){

    return exp(x) + 3 * x * x + x +1;

}

void abs_pes(int &n, double W[], double T[], int *Info){

    double Toler = pow(10, -15);
    unsigned IterMax = 30;
    int m = floor(n / 2);

    double delta;
    unsigned Iter;

    double z;

    double p1;
    double Pz;
    double p0;
    double DPz;

    for(unsigned i = 1; i <= m; i++){

        z = cos((i - 0.25) / (n + 0.5) * M_PI);

        delta = 1 + Toler;
        Iter = 0;

        while(1){

            p1 = 1;
            Pz = z;

            for(unsigned k = 2; k <= n; k++){

                p0 = p1;
                p1 = Pz;
                Pz = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;

            }

            DPz = n * (p1 - z * Pz) / (1 - pow(z, 2));

            if(fabs(delta) <= Toler || Iter == IterMax){

                break;

            }

            delta = Pz/DPz;
            z -= delta;
            Iter++;

        }

        if(fabs(delta) <= Toler){

            T[i - 1] = -1 * z;
            T[n - i] = z;
            W[i - 1] = 2 / ((1 - pow(z, 2)) * pow(DPz, 2));
            W[n - i] = W[i - 1];

        }

        else{

            T[i - 1] = 0;
            T[n - i] = 0;
            W[i - 1] = 0;
            W[n - i] = 0;
            *Info++;

        }

    }

    if(n % 2 != 0){

        double n_alt = n;

        T[m] = 0;
        W[m] = M_PI / 2.0 * pow(tgamma((n_alt + 1) / 2.0) / tgamma(n_alt / 2.0 + 1), 2);

    }

    cout << "- Abscissas:  ";

    for(unsigned i = 0; i < m + 1; i++){

        cout << setprecision(10) << fixed << T[i] << "  ";

        if((i + 1) % 5 == 0){

            cout << endl << "              ";

        }

    }

    cout << endl << "- Pesos    :  ";

    for(unsigned i = 0; i < m + 1; i++){

        cout << setprecision(10) << fixed << W[i] << "  ";

        if((i + 1) % 5 == 0){

            cout << endl << "              ";

        }

    }
    
    cout << endl << endl;

}

void print(int &i, double &x, double &y, double &w, double &t){

    cout << setprecision(5) << fixed << "\t" << i << "     " << t << "     " 
         << x << "     " << y << "     " << w << endl;

}

void gauss_legendre(double &a, double &b, int &n){

    cout << "\nIntegracao numerica via Gauss-Legendre com " << n << " pontos\n" << endl;

    double Integral = 0;

    int Info = 0;

    int m = n + 1;

    double T[m];
    double W[m];

    abs_pes(n, W, T, &Info);

    if(Info != 0){

        return;

    }

    cout << "\ti\tt(i)\t    x(i)       f(x(i))       W(i)" << endl;

    double ba2 = (b - a) / 2;

    double x;
    double y;

    for(int i = 1; i <= n; i++){

        x = a + ba2 * (T[i - 1] + 1);
        y = avaliacao_funcao(x);
        Integral += y * W[i - 1];

        print(i, x, y, W[i - 1], T[i - 1]);

    }

    Integral *= ba2;

    cout << "\nIntegral = " << Integral << endl << endl;

}