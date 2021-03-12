/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Max_Iterations 1000 //设定最大可迭代次数
#define PI 3.1415926535898 //定义pi值
#define UNK 3 //设定未知量的数量
#define TOL 0.05 //设定误差

//计算向量二范数
double Norm(double* x){
    double n = 0;
    for(int i = 1;i<UNK;i++)
        n += x[i]*x[i];
    return sqrt(n);
}

//打印结果
void Print(double* x, double g1){
    for(int i = 0;i<UNK;i++)
        printf("x%d = %.5lf\n",i+1,x[i]);
    printf("g = %.5lf\n",g1);
}

//设定g函数
double g(double* x){
    double *y = (double*)malloc(sizeof(double)*UNK),z = 0;
    y[0] = 10*x[0]-2*x[1]*x[1]+x[1]-2*x[2]-5;
    y[1] = 8*x[1]*x[1]+4*x[2]*x[2]-9;
    y[2] = 8*x[1]*x[2]+4;
    for(int i = 0;i<UNK;i++)
        z += y[i]*y[i];
    return z;
}

//设定g函数的梯度
double* nabla_g(double* x,double* y){
    y[0] = 20*(10*x[0]-2*x[1]*x[1]+x[1]-2*x[2]-5);
    y[1] = -8*x[1]*(10*x[0]-2*x[1]*x[1]+x[1]-2*x[2]-5)+32*x[1]*(8*x[1]*x[1]+4*x[2]*x[2]-9)+16*x[2]*(8*x[1]*x[2]+4);
    y[2] = -4*(10*x[0]-2*x[1]*x[1]+x[1]-2*x[2]-5)+16*x[2]*(8*x[1]*x[1]+4*x[2]*x[2]-9)+16*x[1]*(8*x[1]*x[2]+4);  
    return y;
}

int main(){
    //相关变量初始化和动态空间开辟
    double *x , *z, *temp, gn,g0,g1,g2,g3,z0,an,a0,a1,a2,a3,h1,h2,h3;
    x = (double*)malloc(sizeof(double)*UNK);
    z = (double*)malloc(sizeof(double)*UNK);
    temp = (double*)malloc(sizeof(double)*UNK);
    int k = 1, f = 0;
    
    //设定初始x = {0,0,0}
    for(int i = 0;i<UNK;i++)
        x[i] = 0;

    while(k<Max_Iterations){
        g1 = g(x);
        nabla_g(x,z);
        z0=Norm(z);
        if(z0==0){
            printf("Zero Gradient\n");
            Print(x,g1);
            printf("The procedure completed, might have a minimum.");
            break;
        }

        //将z转化为单位向量
        for(int i = 0;i<UNK;i++)
            z[i] = z[i] / z0;
        a1 = 0;
        a3 = 1;
        for(int i = 0;i<UNK;i++)
            temp[i] = x[i]-a3*z[i];
        g3 = g(temp);

        while(g3>=g1){
            a3 = a3/2;
            for(int i = 0;i<UNK;i++)
                temp[i] = x[i]-a3*z[i];
            g3 = g(temp); 
            if(a3<TOL/2.0){
                printf("No Likely Improvement\n");
                Print(x,g1);
                printf("The procedure completed, might have a minimum.");
                f = 1;
                break;
            }
        }

        //如果在上面的循环中已输出结果，则退出迭代
        if(f == 1)
            break;
            
        a2 = a3/2;
        for(int i = 0;i<UNK;i++)
            temp[i] = x[i]-a2*z[i];
        g2 = g(temp);
        
        h1 = (g2-g1)/a2;
        h2 = (g3-g2)/(a3-a2);
        h3 = (h2-h1)/a3;

        a0 = 0.5*(a2-h1/h3);
        for(int i = 0;i<UNK;i++)
            temp[i] = x[i]-a0*z[i];
        g0 = g(temp);

        gn = g0<g3?g0:g3;
        an = gn==g0?a0:a3;
        for(int i = 0;i<UNK;i++)
            x[i] = x[i]-an*z[i];
        if(fabs(gn-g1)<TOL){
            printf("Finished\n");
            Print(x,gn);
            printf("The procedure was successful.");
            break;
        }
        k++;
    }
    //超出迭代最高次数，迭代失败
    if(k>=Max_Iterations){
        printf("Maximum iterations exceeded\n");
        printf("The procedure was unsuccessful.");
    }
    free(temp);
    free(x);
    free(z);
    return 0;
}
