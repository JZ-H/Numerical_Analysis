#include <stdio.h>
#include <math.h>
#define PI 3.1415926535898

double g(double x){
    double y = sin(PI*x)/2+x/4+x;
    return y;
}

int main(){
    double p, p0 = 1;
    int i = 1;
    while(i<10000){
        p  = g(p0);
        printf("%d: %.2lf\n",i,p);
        if(fabs(p-p0)<0.01){
            printf("x = %.2lf\n", p);
            return 0;
        }
        i++;
        p0 = p;
    }
    printf("CAN NOT CONVEGER!!");
    return 0;
}