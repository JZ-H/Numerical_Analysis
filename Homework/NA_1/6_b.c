/* - - - - coding: GB 2312 - - - - */

#include <stdio.h>
#include <math.h>
#define Max_Iterations 10000

//�趨�����ĺ���
double g(double x,double c){
    //Ϊʹ��������������¸ù��̿�����������ԭ����f(x)=0�任Ϊf(x)/c+x=x
    //���ڲ�ͬ��Χ�ĸ������Ҫ�趨��Ӧ��cֵ
    double y = 3*x*x / c - exp(x) / c + x;
    return y;
}

//�����㷽�����
void FP(double p0,double c){
    double p;
    int i = 1;
    while(i<Max_Iterations){
        p = g(p0,c);//����õ�pi
        //��ӡ��i�ε������
        printf("p%d: %.3lf\n",i,p);
        if(fabs(p-p0)<0.01){
            //��ӡ���ս��
            printf("x = p%d = %.3lf\n",i, p);
            break;
        }
        i++;
        p0 = p;//����p0
    }
    //�жϵ��������������ޣ��������������ʧ��
    if(i>=Max_Iterations)
        printf("Failed after %d iterations",i);
}

int main(){
    //����[-1, 0]�ڵĽ�
    printf("##########\n[-1, 0]\n");
    FP(0, 3);//�趨cΪ3��ʹ��������
    
    //����[0.5, 1.5]�ڵĽ�
    printf("##########\n[0.5, 1.5]\n");
    FP(0.5,-3);//�趨cΪ-3��ʹ��������
    
    //����[3, 4]�ڵĽ�
    printf("##########\n[3, 4]\n");
    FP(3,15);////�趨cΪ15��ʹ��������
    return 0;
}


