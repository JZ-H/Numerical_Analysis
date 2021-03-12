/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//�ҵ������i+1���еľ���ֵ���ֵ
double Find_Max(double** A, int i, int r){
    double max = fabs(A[i][0]);
    for(int j = 1;j < r;j++)
        if(fabs(A[i][j])>max)
            max = fabs(A[i][j]); 
    return max;
}

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

//��ӡ���
void Print(double *x, int r){
    for(int i = 0;i<r;i++)
        printf("x%d = %lf\n",i+1,x[i]);
}


int main(){
    //��ʼ���������������
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
        //���Ҳ�������������p�������޷���ⷵ��NULL
        if(A[p][i] == 0 && printf("No Unique Solution Exists!"))
            return 0;
        //��p+1���ǵ�i+1�У�����Ҫ�������i+1�л���
        if(p!=i){
            temp = A[p];
            A[p] = A[i];
            A[i] = temp;
        }
        //���������Խ������²�����Ϊ0
        for(int j = i+1;j<r;j++){
            m = A[j][i]/A[i][i];
            Subtract_Row(A[j],m,A[i],c);      
        }
    }

    //��ʼ�����滻��ͨ���任��ľ������x
    if(A[r-1][r-1]==0&& printf("No Unique Solution Exists!"))
        return 0;
    Backward_Substitution(r,A,x);
    Print(x, r);
    return 0;    
}

