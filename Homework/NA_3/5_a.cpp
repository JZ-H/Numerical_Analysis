/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Max_Iterations 1000 //�趨����������
#define TOL 0.001 //�趨����������

//���������
double Norm(double* x,int r){
    double n = fabs(x[0]);
    for(int i = 1;i<r;i++)
        if(fabs(x[i])>n)
            n = fabs(x[i]);
    return n;
}

//��ӡ���
void Print(double* x, int r){
    for(int i = 0;i<r;i++)
        printf("x%d = %.5lf\n", i+1, x[i]);
}

int main(){
    //��ʼ����д���������
    int k = 1,r = 3, c = 4;
    double xo[3] = {0,0,0};
    double x[3] = {0,0,0};
    double temp[3] = {0,0,0};
    double a[3][3] = {{3,-1,1},{3,6,2},{3,3,7}};
    double b[3] = {1,0,4};
    double sum1  = 0, sum2 = 0, sum = 0;

    //�ֱ������ַ������е���
    for(int method = 0;method<2;method++){
        if(method == 0)
            printf("--- Using the Jacobi method ---\n");
        else
            printf("--- Using the Gauss-Seidel method ---\n");
        while(k<Max_Iterations){
            if(method == 0)
            //Jacobi method�µĵ�������
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
            else
            //Gauss-Seidel method�µĵ�������
                for(int i = 0;i<r;i++){
                    sum1 = 0;
                    sum2 = 0;
                    for(int j = 0;j<i;j++)
                        sum1 += a[i][j]*x[j];
                    for(int j = i+1;j<r;j++)
                        sum2 += a[i][j]*xo[j];
                    x[i] = 1.0/a[i][i]*(b[i]-sum1-sum2);
                }
            //printf("--- After %d Iterations ---\n",k);
            //Print(x,r);
            for(int i = 0;i<r;i++)
                temp[i] = x[i] - xo[i];
            if(Norm(temp, r)<TOL){
                printf("--- Finished Successfully After %d Iterations ---\n",k);
                Print(x,r);
                break;
            }
            k++;
            for(int i = 0;i<r;i++)
                xo[i] = x[i];
        }
        //������������������ޣ����������ʧ��
        if(k>=Max_Iterations)
            printf("Maximum number of iterations exceeded");
    }
    return 0;
}