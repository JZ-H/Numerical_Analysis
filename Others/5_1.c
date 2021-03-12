/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <math.h>
#define Max_Iterations 10000

//�趨�����ĺ���
double Func(double x){
    double y = exp(x)-x*(x-3)-2;
    return y;
}

int main(){
    double e, p, a, b,FP,FA;
    int i = 1;
    //���ַ�����
    a = 0;
    b = 1;
    FA = Func(a);
    while(i<Max_Iterations){
        p = a + (b-a) / 2;
        //��ӡ��i�ε������
        printf("p%d: %.6lf\n",i,p);
        FP = Func(p);
        if(FP == 0 || (b-a) < 0.00001){
            //��ӡ���ս��
            printf("x = p%d = %.6lf\n",i,p);
            return 0;
        }
        if(FA*FP>0){
            a = p ;
            FA = FP;
        }
        else
            b = p;
        i++;
    }
    //�жϵ��������������ޣ��������������ʧ��
    if(i>=Max_Iterations)
        printf("Failed after %d iterations",i);
    return 0;
}