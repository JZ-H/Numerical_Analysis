#include <stdio.h>
#include <math.h>

double g(double x){
    double y = 3*x*x / 10 - exp(x) / 10 + x;
    return y;
}

int main(){
    double p, p0 = 3.5;
    int i = 1;
    while(i<100){
        p  = g(p0);
        printf("%d: %.2lf\n",i,p);
        if(fabs(p-p0)<0.01){
            printf("ANS: %.2lf\n", p);
            return 0;
        }
        i++;
        p0 = p;
    }
    printf("CAN NOT CONVEGER!!");
    return 0;
}