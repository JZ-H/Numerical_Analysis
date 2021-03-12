#include <stdio.h>
#include <math.h>

double Func(double x){
    double y = exp(x)-x*(x-3)-2;
    return y;
}
/*

double Func(double x){
    double y = x*cos(x)-x*(2*x-3)-1;
    return y;
}

*/
int main(){
    double e, p, a, b,FP,FA;
    int i = 1;
    a = 0;
    b = 1;
    FA = Func(a);
    while(1){
        p = a + (b-a) / 2;
        printf("%d: %.5lf\n",i,p);
        FP = Func(p);
        if(FP == 0 || (b-a)/2  < 0.00001){
            printf("x = %.5lf\n",p);
            return 0;
        }
        if(FA*FP>0){
            a = p ;
            FA = FP;
        }
        else
            b = p;
        i++;
    }
    return 0;
}