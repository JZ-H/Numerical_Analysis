#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.1415926535898

int main() {
    int L = 10000, T = 3, m = 100000, n = 300;
    double h = 0.1, k = 0.01, a = 1;
    double lambda = a * a * k / (h * h);
    double* l = (double*)malloc(sizeof(double) * m);
    double* u = (double*)malloc(sizeof(double) * m);
    double* z = (double*)malloc(sizeof(double) * m);
    double* w = (double*)malloc(sizeof(double) * m);
    
    w[0] = 10000;
    for (int i = 1; i < m - 1; i++) {
        w[i] = 0;
    }

    l[0] = 1 + 2 * lambda;
    u[0] = -lambda / l[0];
    for (int i = 1; i < m - 2; i++) {
        l[i] = 1 + 2 * lambda + lambda * u[i - 1];
        u[i] = -lambda / l[i];
    }

    l[m - 2] = 1 + 2 * lambda + lambda * u[m - 3];
    for (int i = 0; i < n; i++) {
        double t = (i + 1) * k;
        z[0] = w[0] / l[0];
        for(int j = 1;j<m-1;j++)
            z[j] = (w[j]+lambda*z[j-1])/l[j];
        w[m-2] = z[m-2];
        for(int j = m-3;j>=0;j--)
            w[j] = z[j]-u[j]*w[j+1];
        //printf("t = %lf\n", t);
        for(int j = 0;j<m-1;j++){
            double x = (j+1)*h;
            if(t==3)
            printf("x = %lf, w[%d] = %lf\n",x,j+1,w[j]+1);
            if(t==3&&j>40)
            break;
        }
    }

    return 0;

}