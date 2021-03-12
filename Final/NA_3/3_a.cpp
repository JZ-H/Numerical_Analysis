/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//找到矩阵第i+1行中的绝对值最大值
double Find_Max(double** A, int i, int r){
    double max = fabs(A[i][0]);
    for(int j = 1;j < r;j++)
        if(fabs(A[i][j])>max)
            max = fabs(A[i][j]); 
    return max;
}

//第j+1行减去乘上m后的i+1行，目的在于使第j+1行且第i+1列的元素为0
void Subtract_Row(double *Aj,double m,double *Ai,int c){
    for(int i = 0;i<c;i++)
        Aj[i] -= m*Ai[i];
}

//交换第p+1行和第ip行
void Swap_Row(double**A,int p,int ip,int c){
    double t;
    for(int j =0;j<c;j++){
        t = A[p][j];
        A[p][j] = A[ip][j];
        A[ip][j] = t;
    }
}

//
int Find_P(double **A,double* s,int i,int r){
    int p = i;
    double max = fabs(A[i][i]) / s[i];
    for(int k = i;k<r;k++)
        if(fabs(A[k][i]) / s[k] > max){
            p = k;
            max = fabs(A[k][i]) / s[k];
        }
    return p;
}

//通过反向替换求解x
double* Backward_Substitution(int r,double** A,double* ans){
    double sum_temp;
    ans[r-1] = A[r-1][r]/A[r-1][r-1];
    for(int i = r-2;i>=0;i--){
        sum_temp = 0;
        for(int j = i+1;j<r;j++)
            sum_temp+=A[i][j]*ans[j];
        ans[i]=(A[i][r] - sum_temp)/A[i][i];
        }
    return ans;
}

//打印结果
void Print(double *x, int r){
    for(int i = 0;i<r;i++)
        printf("x%d = %lf\n",i+1,x[i]);
}


int main(){
    //初始化，读入矩阵数据
    int r = 2,c = 3,p;
    double m,*x,*temp;
    x = (double*)malloc(sizeof(double)*r);

    double A1[3] = {0.03, 58.9, 59.2};
    double A2[3] = {5.31, -6.10, 47.0};
    double* A[2] = {A1,A2};
    double s[2];
    for(int i = 0;i<r;i++){
        s[i] = Find_Max(A, i, r);
        if(s[i] == 0 && printf("No Unique Solution Exists!"))
            return 0;
    }

    //Gaussian elimination
    for(int i =0;i<r-1;i++){
        p = Find_P(A, s, i, r);
        //若找不到满足条件的p，方程无法求解返回NULL
        if(A[p][i] == 0 && printf("No Unique Solution Exists!"))
            return 0;
        //若p+1不是第i+1行，则需要将其与第i+1行互换
        if(p!=i){
            temp = A[p];
            A[p] = A[i];
            A[i] = temp;
        }
        //将矩阵主对角线以下部分消为0
        for(int j = i+1;j<r;j++){
            m = A[j][i]/A[i][i];
            Subtract_Row(A[j],m,A[i],c);      
        }
    }

    //开始反向替换，通过变换后的矩阵求解x
    if(A[r-1][r-1]==0&& printf("No Unique Solution Exists!"))
        return 0;
    Backward_Substitution(r,A,x);
    Print(x, r);
    return 0;    
}

