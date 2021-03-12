#include <stdio.h>
#include <math.h>
#define PI 3.1415926535898

double g(double x,double c){
    double y = 3*x*x / c - exp(x) / c + x;
    return y;
}

void FP(double p0,double c){
    double p;
    int i = 1;
    while(i<10000){
        p  = g(p0,c);
        printf("%d: %.2lf\n",i,p);
        if(fabs(p-p0)<0.01){
            printf("x = %.2lf\n", p);
            break;
        }
        i++;
        p0 = p;
    }
    if(i>=10000)
        printf("CAN NOT CONVERGE!!");
}

int main(){
    printf("##########\n[-1, 0]\n");
    FP(0, 3);
    printf("##########\n[0.5, 1.5]\n");
    FP(0.5,-3);
    printf("##########\n[3, 4]\n");
    FP(3,15);
    return 0;
}