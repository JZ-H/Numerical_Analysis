/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Max_Iterations 1000 //设定最大迭代次数
#define TOL 0.001 //设定迭代终止误差限

//计算无穷范数
double Norm(double* x,int r){
    double n = fabs(x[0]);
    for(int i = 1;i<r;i++)
        if(fabs(x[i])>n)
            n = fabs(x[i]);
    return n;
}

//打印结果
void Print(double* x, int r){
    for(int i = 0;i<r;i++)
        printf("x%d = %.5lf\n", i+1, x[i]);
}

int main(){
    //初始化，写入相关数据
    int k = 1,r = 3, c = 4;
    double xo[3] = {0,0,0};
    double x[3] = {0,0,0};
    double temp[3] = {0,0,0};
    double a[3][3] = {{-2,1,0.5},{1,-2,-0.5},{0,1,2}};
    double b[3] = {4,-4,0};
    double sum  = 0;

    //the Jacobi method
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
        //打印前三次迭代结果
        if(k<=3){
            printf("--- After %d Iterations ---\n",k);
            Print(x,r);
        }
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
    //如果迭代次数超过上限，则输出迭代失败
    if(k>=Max_Iterations)
        printf("Maximum number of iterations exceeded.");
    return 0;
}