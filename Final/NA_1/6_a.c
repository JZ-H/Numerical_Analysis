/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <math.h>
#define PI 3.1415926535898
#define Max_Iterations 10000

//设定待求解的函数
double g(double x){
    //为使不动点迭代方法下该过程可以收敛，将原函数左右两端除以4后再加x
    double y = sin(PI*x)/2+x/4+x;
    return y;
}

int main(){
    //不动点方法求解
    double p, p0 = 1;
    int i = 1;
    while(i<Max_Iterations){
        p  = g(p0);//计算得到pi
        //打印第i次迭代结果
        printf("p%d: %.3lf\n",i,p);
        if(fabs(p-p0)<0.01){
            //打印最终结果
            printf("x = p%d = %.3lf\n",i, p);
            return 0;
        }
        i++;
        p0 = p;//更新p0
    }
    //判断迭代次数超过上限，若超过，则迭代失败
    if(i>=Max_Iterations)
        printf("Failed after %d iterations",i);
    return 0;
}
