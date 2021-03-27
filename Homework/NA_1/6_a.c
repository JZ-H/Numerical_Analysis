/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <math.h>
#define PI 3.1415926535898
#define Max_Iterations 10000

//�趨�����ĺ���
double g(double x){
    //Ϊʹ��������������¸ù��̿�����������ԭ�����������˳���4���ټ�x
    double y = sin(PI*x)/2+x/4+x;
    return y;
}

int main(){
    //�����㷽�����
    double p, p0 = 1;
    int i = 1;
    while(i<Max_Iterations){
        p  = g(p0);//����õ�pi
        //��ӡ��i�ε������
        printf("p%d: %.3lf\n",i,p);
        if(fabs(p-p0)<0.01){
            //��ӡ���ս��
            printf("x = p%d = %.3lf\n",i, p);
            return 0;
        }
        i++;
        p0 = p;//����p0
    }
    //�жϵ��������������ޣ��������������ʧ��
    if(i>=Max_Iterations)
        printf("Failed after %d iterations",i);
    return 0;
}
