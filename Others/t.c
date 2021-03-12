#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Max_Iterations 10 //设定最大可迭代次数
#define PI 3.141592653589793 //定义pi值
#define UNK 4 //设定未知量的数量
#define TOL 0.1 //设定误差

void Subtract_Row(double *Aj,double m,double *Ai,int c);
void Swap_Row(double**A,int p,int ip,int c);
int Find_Min(int r,int i,double **A);
double* Backward_Substitution(int r,double** A,double* ans);
double* Gauss(double** x,double* y,double* ans);
void Print(double* x);
double* F(double* x,double* y);
double** J(double* x,double** y);

int main(){
    double  *y, *neg_f,**j;
    j = (double**)malloc(sizeof(double*)*UNK);
    for(int i = 0;i<UNK;i++)
        j[i] = (double*)malloc(sizeof(double)*UNK);
    y = (double*)malloc(sizeof(double)*UNK);
    neg_f = (double*)malloc(sizeof(double)*UNK);
    for(int i = 0;i<4;i++)
        for(int x = 0;x<5;x++)
            if(x!=4)
                scanf("%lf",&j[i][x]);
            else
                scanf("%lf",&neg_f[i]);
    Gauss(j,neg_f,y);
    return 0;
}



//设定F函数
double* F(double* x,double* y){
    y[0] = 3*x[0]-cos(x[1]*x[2])-1.0/2.0;
    y[1] = 4*x[0]*x[0]-625*x[1]*x[1]+2*x[1]-1;
    y[2] = exp(-x[0]*x[1])+20*x[2]+10.0*PI/3.0-1;
    return y;
}

//设定Jacobian矩阵
double** J(double* x,double** y){
    y[0][0] = 3;
    y[0][1] = x[2]*sin(x[1]*x[2]);
    y[0][2] = x[1]*sin(x[1]*x[2]);
    y[1][0] = 8*x[0];
    y[1][1] = -1250*x[1]+2;
    y[1][2] = 0;
    y[2][0] = -x[1]*exp(-x[0]*x[1]);
    y[2][1] = -x[0]*exp(-x[0]*x[1]);
    y[2][2] = 20;
    return y;
}

//打印x的结果
void Print(double* x){
    for(int i = 0;i<UNK;i++)
        printf("x%d = %.6lf\n",i+1,x[i]);
}


/*--- 以下为 Gaussian elimination 相关函数 ---*/

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

//在第i+1列找到最小的j(i<=j<r)且满足该列第j+1行的元素非0
//如果没有返回-1
int Find_Min(int r,int i,double **A){
    int min = -1;
    for(int j = i;j<r;j++)
        if(A[j][i]!=0){
            min = j;
            break;
        }
    return min;
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

//Gaussian elimination主体函数
double* Gauss(double** x,double* y,double* ans){
    //初始化，读入矩阵数据，r为行数，c为列数
    int r = UNK,c = UNK + 1,p;
    double m;

    //构建增广矩阵
    double **A=(double**)malloc(sizeof(double*)*r);
    for(int i = 0;i<r;i++)
        A[i] = (double*)malloc(sizeof(double)*c);
    for(int i = 0;i<r;i++)
        for(int j = 0;j<c;j++){
            if(j != c-1)
                A[i][j] = x[i][j];
            else
                A[i][j] = y[i];
        }
    //Gaussian elimination
    for(int i =0;i<r-1;i++){
        p = Find_Min(r,i,A);
        //若找不到满足条件的p，方程无法求解返回NULL
        if(p==-1 && printf("No Unique Solution Exists!"))
            return NULL;
        //若p不在第i+1行，则需要将其换到i+1行
        if(p!=i)
            Swap_Row(A,p,i,c);
        //将矩阵主对角线以下部分消为0
        for(int j = i+1;j<r;j++){
            m = A[j][i]/A[i][i];
            Subtract_Row(A[j],m,A[i],c);      
        }
    }

    //开始反向替换，通过变换后的矩阵求解x
    if(A[r-1][r-1]==0&& printf("No Unique Solution Exists!"))
        return NULL;
    for(int i = 0;i<4;i++){
    for(int j = 0;j<5;j++)
        printf("%.5lf ",A[i][j]);printf("\n");}
    Backward_Substitution(r,A,ans);
    for(int i = 0;i<r;i++)
        free(A[i]);
    for(int i = 0;i<4;i++)
        printf("%.5lf ",ans[i]);
    return ans;    
}