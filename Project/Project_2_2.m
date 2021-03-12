% coding: utf-8
% description: 应用Crank-Nicolson方法求解限定源扩散问题
clear;clc;
% 相关参数的初始化
L=10000; T=3; m=100000; n=30;
h=L/m; k=T/n; a=1;
lambda = a * a * k / (h * h);
l=zeros(m-1,1);
u=zeros(m-1,1);
z=zeros(m-1,1);
w=zeros(m,1);
for i = 1:m
    % 将L=5000处设为扩散开始处，并把t=0.1时各点浓度数值作为初始值
    w(i) = 1/sqrt(pi*0.1)*exp(-((50000-i)*0.1)^2/(4*0.1));
end

% 初始化用于存储结果数据的向量
x_ = zeros(3 / h * 3 / h,1);	%记录x值数据
y_ = zeros(3 / h * 3 / h,1);	%记录t值数据
u_ = zeros(3 / h * 3 / h,1);	%记录u值数据

% 记录t=0.1时的数据
for j = 1:3 / h
    x_(j) = j*h;
    y_(j) = 0.1;
    u_(j) = w(50000+j);
end

% 求解线性方程组
l(1) = 1 + lambda;
u(1) = -lambda/(2*l(1));
for i = 2:m-2
    l(i) = 1 + lambda + lambda*u(i-1)/2;
    u(i) = -lambda/(2*l(i));
end
l(m-1) = 1 + lambda + lambda*u(m-2)/2;
for i = 1:n - 1
    t = (i+1)*k;
    z(1) = ((1-lambda)*w(1)+lambda/2*w(2))/l(1);
    for j = 2:m-1
        z(j)=((1-lambda)*w(j)+lambda/2*(w(j+1)+w(j-1)+z(j-1)))/l(j);
    end
    w(m-1) = z(m-1);
    for j = m-2:-1:1
        w(j) = z(j) - u(j)*w(j+1);
    end

    % 记录计算结果
    for j = 1:3 / h
        x_(j+i*3 / h) = j*h;
        y_(j+i*3 / h) = t;
        u_(j+i*3 / h) = w(50000+j);
    end
end

% 输出结果
r_name=["x";"t";"u"];
result_table = table(x_,y_,u_,'VariableNames',r_name);
writetable(result_table,"Project_2_2.2.csv");

% 作图，应用散乱点插值方法
[X,T]=meshgrid(linspace(min(x_),max(x_)),linspace(min(y_),max(y_)));
U=griddata(x_,y_,u_,X,T,'v4');
mesh(X,T,U);
xlabel('x');
ylabel('t');
zlabel('u');