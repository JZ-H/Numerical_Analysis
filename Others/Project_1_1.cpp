#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.1415926535898

int main() {
    int L = 10, T = 3, m = 100, n = 30;
    double h = 0.1, k = 0.1, a = 1;
    double lambda = a * a * k / (h * h);
    double* l = (double*)malloc(sizeof(double) * (m-1));
    double* u = (double*)malloc(sizeof(double) * (m-1));
    double* z = (double*)malloc(sizeof(double) * (m-1));
    double* w = (double*)malloc(sizeof(double) * (m-1));
    for (int i = 0; i < m - 1; i++) {
        w[i] = -1;
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
            if(t==1.6)
            printf("u(x = %.3lf, t = %.3lf) = %lf\n",x,t,w[j]+1);

        }
    }
    return 0;

}