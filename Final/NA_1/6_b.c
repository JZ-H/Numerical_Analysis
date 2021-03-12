/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <math.h>
#define Max_Iterations 10000

//设定待求解的函数
double g(double x,double c){
    //为使不动点迭代方法下该过程可以收敛，将原方程f(x)=0变换为f(x)/c+x=x
    //对于不同范围的根求解需要设定相应的c值
    double y = 3*x*x / c - exp(x) / c + x;
    return y;
}

//不动点方法求解
void FP(double p0,double c){
    double p;
    int i = 1;
    while(i<Max_Iterations){
        p = g(p0,c);//计算得到pi
        //打印第i次迭代结果
        printf("p%d: %.3lf\n",i,p);
        if(fabs(p-p0)<0.01){
            //打印最终结果
            printf("x = p%d = %.3lf\n",i, p);
            break;
        }
        i++;
        p0 = p;//更新p0
    }
    //判断迭代次数超过上限，若超过，则迭代失败
    if(i>=Max_Iterations)
        printf("Failed after %d iterations",i);
}

int main(){
    //计算[-1, 0]内的解
    printf("##########\n[-1, 0]\n");
    FP(0, 3);//设定c为3，使其能收敛
    
    //计算[0.5, 1.5]内的解
    printf("##########\n[0.5, 1.5]\n");
    FP(0.5,-3);//设定c为-3，使其能收敛
    
    //计算[3, 4]内的解
    printf("##########\n[3, 4]\n");
    FP(3,15);////设定c为15，使其能收敛
    return 0;
}


