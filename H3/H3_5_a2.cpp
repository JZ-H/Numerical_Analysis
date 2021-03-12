#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Max_Iterations 1000
#define TOL 0.001

double Norm(double* x,int r){
    double n = fabs(x[0]);
    for(int i = 1;i<r;i++)
        if(fabs(x[i])>n)
            n = fabs(x[i]);
    return n;
}

void Print(double* x, int r){
    for(int i = 0;i<r;i++)
        printf("x%d = %.5lf\n", i+1, x[i]);
}

int main(){
    int k = 1,r = 3, c = 4;
    double xo[3] = {0,0,0};
    double x[3] = {0,0,0};
    double temp[3] = {0,0,0};
    double a[3][3] = {{3,-1,1},{3,6,2},{3,3,7}};
    double b[3] = {1,0,4};
    double sum  = 0;
    while(k<Max_Iterations){
        for(int i = 0;i<r;i++){
            sum = 0;
            for(int j = 0;j<r;j++){
                if(j==i)
                    continue;
                else
                    sum += a[i][j]*xo[j];
            }
            x[i] = 1.0/a[i][i]*(b[i]-sum);
        }
        printf("--- After %d Iterations ---\n",k);
        Print(x,r);
        for(int i = 0;i<r;i++)
            temp[i] = x[i] - xo[i];
        if(Norm(temp, r)<TOL){
            printf("--- Finished After %d Iterations ---\n",k);
            Print(x,r);
            printf("The procedure was successful.");
            break;
        }
        k++;
        for(int i = 0;i<r;i++)
            xo[i] = x[i];
    }
    if(k>=Max_Iterations)
        printf("Maximum number of iterations exceeded");
    return 0;
}