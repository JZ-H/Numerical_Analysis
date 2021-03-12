#include <stdio.h>
#include <stdlib.h>

void Read(int r, int c,double **p){
    for(int i = 0;i < r;i++)
        for(int j = 0;j < c;j++)
            scanf("%lf",&p[i][j]);
}

void Subtract_Row(double *Aj,double m,double *Ai,int c){
    for(int i = 0;i<c;i++)
        Aj[i] -= m*Ai[i];
}

void Swap_Row(double**A,int p,int ip,int c){
    double t;
    for(int j =0;j<c;j++){
        t = A[p][j];
        A[p][j] = A[ip][j];
        A[ip][j] = t;
    }
}

int Find_Min(int r,int i,double **A){
    int min = -1;
    for(int j = i;j<r;j++)
        if(A[j][i]!=0){
            min = j;
            break;
        }
    return min;
}

double Sum(int j,int n,double** A,double *x){
    double s=0;
    for(int i = j;i<=n;i++)
        s+=A[j-1][i]*x[i];
    return s;
}

void Print(double *x, int r){
    for(int i = 0;i<r;i++)
        printf("x%d = %.4lf\n",i+1,x[i]);
}

void Print_M(double ** A,int r,int c){
    for(int i = 0;i<r;i++){
        for(int j =0;j<c;j++)
            printf("%.4lf ",A[i][j]);
        printf("\n");
    }
    for(int i = 0;i<r*2;i++)
        printf("#");
    printf("\n");
}

int main(){
    // 初始化，读入矩阵数据
    int r,c,p;
    double m,*x;
    scanf("%d %d",&r,&c);
    if(r!=c-1 && printf("Input Error!"))
        return 0;
    double **A=(double**)malloc(sizeof(double*)*r);
    x = (double*)malloc(sizeof(double)*r);
    for(int i = 0;i<r;i++)
        A[i] = (double*)malloc(sizeof(double)*c);
    Read(r,c,A); 
    
    // Gaussian elimination
    for(int i =0;i<r-1;i++){
        p = Find_Min(r,i,A);
        if(p==-1 && printf("No Unique Solution Exists!"))
            return 0;
        if(p!=i)
            Swap_Row(A,p,i,c);
        for(int j = i+1;j<r;j++){
            m = A[j][i]/A[i][i];
            Subtract_Row(A[j],m,A[i],c);      
        }
        //Print_M(A,r,c);
    }

    // 通过变换后的矩阵求解x
    if(A[r-1][r-1]==0&& printf("No Unique Solution Exists!"))
        return 0;
    x[r-1] = A[r-1][r]/A[r-1][r-1];
    for(int i = r-2;i>=0;i--)
        x[i]=(A[i][r] - Sum(i+1, r-1,A,x))/A[i][i];
    Print(x,r);

    // 释放动态内存
    free(x);
    for(int i = 0;i<r;i++)
        free(A[i]);
    return 0;    
}