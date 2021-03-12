#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Max_Iterations 10 //设定最大可迭代次数
#define PI 3.141592653589793 //定义pi值
#define UNK 3 //设定未知量的数量
#define TOL 0.1 //设定误差
void Subtract_Row(double *Aj,double m,double *Ai,int c);
void Swap_Row(double**A,int p,int ip,int c);
int Find_Min(int r,int i,double **A);
double Sum(int j,int n,double** A,double *x);
double* Gauss(double** x,double* y,double* ans);

//计算向量无穷范数
double Norm(double* x){
    double max = x[0];
    for(int i = 1;i<UNK;i++)
        if(fabs(x[i])>max)
            max = fabs(x[i]);
    return max;
}

void Print(double* x){
    for(int i = 0;i<UNK;i++)
        printf("x%d = %.5lf\n",i+1,x[i]);
}

double* F(double* x,double* y){
    y[0] = 3*x[0]-cos(x[1]*x[2])-1.0/2.0;
    y[1] = 4*x[0]*x[0]-625*x[1]*x[1]+2*x[1]-1;
    y[2] = exp(-x[0]*x[1])+20*x[2]+10.0*PI/3.0-1;
    return y;
}

double** J(double* x){
    double** y = (double**)malloc(sizeof(double*)*UNK);
    for(int i = 0;i<UNK;i++)
        y[i] = (double*)malloc(sizeof(double)*UNK);
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

int main(){
    double *x , *y, *neg_f;
    x = (double*)malloc(sizeof(double)*UNK);
    y = (double*)malloc(sizeof(double)*UNK);
    neg_f = (double*)malloc(sizeof(double)*UNK);
    for(int i = 0;i<UNK;i++)
        x[i] = 0;
    int i = 1;
    while(i<Max_Iterations){
        F(x, neg_f);
        for(int i = 0;i<UNK;i++)
            neg_f[i] = -neg_f[i];
        Gauss(J(x),neg_f,y);
        for(int i = 0;i<UNK;i++)
            x[i] = x[i] + y[i];
        printf("---After %d Iterations---\n",i);
        Print(x);
        if(i==2){
            break;
        }
        i++;
    }
    return 0;
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

//在第i+1列找到最小的j(i<=j<r)满足该列第j+1行的元素非0
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

//求解公式中的一部分
//计算第j行从j+1列到最后一列的每个元素上相应x_{j+1}的和
double Sum(int j,int n,double** A,double *x){
    double s=0;
    for(int i = j;i<=n;i++)
        s+=A[j-1][i]*x[i];
    return s;
}

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
    
    // Gaussian elimination
    for(int i =0;i<r-1;i++){
        p = Find_Min(r,i,A);
        //若找不到满足条件的p，方程无法求解返回NULL
        if(p==-1 && printf("No Unique Solution Exists!"))
            return NULL;
        //若p不在第i+1行，则需要将其换到i+1行
        if(p!=i)
            Swap_Row(A,p,i,c);
        //消元
        for(int j = i+1;j<r;j++){
            m = A[j][i]/A[i][i];
            Subtract_Row(A[j],m,A[i],c);      
        }
    }

    // 通过变换后的矩阵求解x
    if(A[r-1][r-1]==0&& printf("No Unique Solution Exists!"))
        return NULL;
    ans[r-1] = A[r-1][r]/A[r-1][r-1];
    for(int i = r-2;i>=0;i--)
        ans[i]=(A[i][r] - Sum(i+1, r-1,A,ans))/A[i][i];
    for(int i = 0;i<r;i++)
        free(A[i]);
    return ans;    
}