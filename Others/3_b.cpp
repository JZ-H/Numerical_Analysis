/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592653589793 //����piֵ
#define UNK 3 //�趨δ֪��������
#define Iterations 2 //�趨��������

//��������
void Subtract_Row(double *Aj,double m,double *Ai,int c);
void Swap_Row(double**A,int p,int ip,int c);
int Find_Min(int r,int i,double **A);
double* Backward_Substitution(int r,double** A,double* ans);
double* Gauss(double** x,double* y,double* ans);
void Print(double* x);
double* F(double* x,double* y);
double** J(double* x,double** y);

int main(){
    //��ر�����ʼ���Ͷ�̬�ռ俪��
    double *x , *y, *neg_f,**j_matrix;
    j_matrix = (double**)malloc(sizeof(double*)*UNK);
    for(int i = 0;i<UNK;i++)
        j_matrix[i] = (double*)malloc(sizeof(double)*UNK);
    x = (double*)malloc(sizeof(double)*UNK);
    y = (double*)malloc(sizeof(double)*UNK);
    neg_f = (double*)malloc(sizeof(double)*UNK);
    int k = 1;

    //�趨��ʼx = {0,0,0}
    for(int i = 0;i<UNK;i++)
        x[i] = 0;
    
    while(1){
        //����F(x)��Jacobian����
        F(x, neg_f);
        J(x,j_matrix);

        //��ⷽ�� Jy = -F �е�y
        for(int i = 0;i<UNK;i++)
            neg_f[i] = -neg_f[i];
        Gauss(j_matrix,neg_f,y);

        for(int i = 0;i<UNK;i++)
            x[i] = x[i] + y[i];

        //��ӡ�������
        printf("---After %d Iterations---\n",k);
        Print(x);
        if(k==Iterations){
            break;
        }
        k++;
    }
    free(j_matrix);
    free(y);
    free(x);
    free(neg_f);
    return 0;
}

//�趨F����
double* F(double* x, double* y){
    y[0] = x[0]*x[0]+x[1]-37;
    y[1] = x[0]-x[1]*x[1]-5;
    y[2] = x[0]+x[1]+x[2]-3;
    return y;
}

//�趨Jacobian����
double** J(double* x, double** y){
    y[0][0] = 2*x[0];
    y[0][1] = 1;
    y[0][2] = 0;
    y[1][0] = 1;
    y[1][1] = -2*x[1];
    y[1][2] = 0;
    y[2][0] = 1;
    y[2][1] = 1;
    y[2][2] = 1;
    return y;
}

//��ӡx�Ľ��
void Print(double* x){
    for(int i = 0;i<UNK;i++)
        printf("x%d = %.6lf\n",i+1,x[i]);
}


/*--- ����Ϊ Gaussian elimination ��غ��� ---*/

//��j+1�м�ȥ����m���i+1�У�Ŀ������ʹ��j+1���ҵ�i+1�е�Ԫ��Ϊ0
void Subtract_Row(double *Aj,double m,double *Ai,int c){
    for(int i = 0;i<c;i++)
        Aj[i] -= m*Ai[i];
}

//������p+1�к͵�ip��
void Swap_Row(double**A,int p,int ip,int c){
    double t;
    for(int j =0;j<c;j++){
        t = A[p][j];
        A[p][j] = A[ip][j];
        A[ip][j] = t;
    }
}

//�ڵ�i+1���ҵ���С��j(i<=j<r)��������е�j+1�е�Ԫ�ط�0
//���û�з���-1
int Find_Min(int r,int i,double **A){
    int min = -1;
    for(int j = i;j<r;j++)
        if(A[j][i]!=0){
            min = j;
            break;
        }
    return min;
}

//ͨ�������滻���x
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

//Gaussian elimination���庯��
double* Gauss(double** x,double* y,double* ans){
    //��ʼ��������������ݣ�rΪ������cΪ����
    int r = UNK,c = UNK + 1,p;
    double m;

    //�����������
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
        //���Ҳ�������������p�������޷���ⷵ��NULL
        if(p==-1 && printf("No Unique Solution Exists!"))
            return NULL;
        //��p���ڵ�i+1�У�����Ҫ���任��i+1��
        if(p!=i)
            Swap_Row(A,p,i,c);
        //���������Խ������²�����Ϊ0
        for(int j = i+1;j<r;j++){
            m = A[j][i]/A[i][i];
            Subtract_Row(A[j],m,A[i],c);      
        }
    }

    //��ʼ�����滻��ͨ���任��ľ������x
    if(A[r-1][r-1]==0&& printf("No Unique Solution Exists!"))
        return NULL;
    Backward_Substitution(r,A,ans);
    for(int i = 0;i<r;i++)
        free(A[i]);
    return ans;    
}