#include <stdio.h>
#include <math.h>

double Func(double x){
    double y = x*cos(x)-x*(2*x-3)-1;
    return y;
}


void Bisection(double a,double b){
    double e, p,FP,FA;
    int i = 1;
    FA = Func(a);
    while(i<10000){
        p = a + (b-a) / 2;
        printf("%d: %.5lf\n",i,p);
        FP = Func(p);
        if(FP == 0 || (b-a)/2 < 0.00001){
            printf("x = %.5lf\n",p);
            break;
        }
        if(FA*FP>0){
            a = p ;
            FA = FP;
        }
        else
            b = p;
        i++;
    }
    if(i>=10000)
        printf("CAN NOT CONVERGE!!");
}

int main(){
    printf("##########\n[0.2, 0.3]\n");
    Bisection(0.2,0.3);
    printf("##########\n[1.2, 1.3]\n");
    Bisection(1.2,1.3);
    return 0;
}