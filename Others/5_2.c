/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <math.h>
#define Max_Iterations 10000

//设定待求解的函数
double Func(double x){
    double y = x*cos(x)-x*(2*x-3)-1;
    return y;
}

//二分法求解
void Bisection(double a,double b){
    double e, p,FP,FA;
    int i = 1;
    FA = Func(a);
    while(i<Max_Iterations){
        p = a + (b-a) / 2;
        //打印第i次迭代结果
        printf("p%d: %.6lf\n",i,p);
        FP = Func(p);
        if(FP == 0 || (b-a)/2 < 0.00001){
            //打印最终结果
            printf("x = p%d = %.6lf\n",i,p);
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
    //判断迭代次数超过上限，若超过，则迭代失败
    if(i>=Max_Iterations)
        printf("Failed after %d iterations",i);
}

int main(){
    //用二分法求解[0.2, 0.3]区间上的解
    printf("##########\n[0.2, 0.3]\n");
    Bisection(0.2,0.3);
    
    //用二分法求解[1.2, 1.3]区间上的解
    printf("##########\n[1.2, 1.3]\n");
    Bisection(1.2,1.3);
    return 0;
}

