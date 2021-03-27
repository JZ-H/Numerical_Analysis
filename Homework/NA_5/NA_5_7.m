clear; clc;
x = [0,2,4,5];
y = [6,8,14,20];

% �������a0��a1�л��õ���ֵ
sum_x2 = sum(x.^2);
sum_x = sum(x);
sum_y = sum(y);
sum_xy = sum(x.*y);
len = length(x);

% �����a0��a1
a0 = (sum_x2*sum_y-sum_xy*sum_x)/(len*sum_x2-sum_x^2);
a1 = (sum_xy*len-sum_x*sum_y)/(len*sum_x2-sum_x^2);
fprintf("y = %.6fx + %.6f\n",a1,a0)

% �������ϳ���yֵ
y_hat = a1*x+a0;
E = sum((y-y_hat).^2);
fprintf("E = %.6f\n",E)