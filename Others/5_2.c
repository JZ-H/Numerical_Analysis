/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <math.h>
#define Max_Iterations 10000

//�趨�����ĺ���
double Func(double x){
    double y = x*cos(x)-x*(2*x-3)-1;
    return y;
}

//���ַ����
void Bisection(double a,double b){
    double e, p,FP,FA;
    int i = 1;
    FA = Func(a);
    while(i<Max_Iterations){
        p = a + (b-a) / 2;
        //��ӡ��i�ε������
        printf("p%d: %.6lf\n",i,p);
        FP = Func(p);
        if(FP == 0 || (b-a)/2 < 0.00001){
            //��ӡ���ս��
            printf("x = p%d = %.6lf\n",i,p);
            break;
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
}

int main(){
    //�ö��ַ����[0.2, 0.3]�����ϵĽ�
    printf("##########\n[0.2, 0.3]\n");
    Bisection(0.2,0.3);
    
    //�ö��ַ����[1.2, 1.3]�����ϵĽ�
    printf("##########\n[1.2, 1.3]\n");
    Bisection(1.2,1.3);
    return 0;
}

