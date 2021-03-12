% coding: utf-8
% description: 计算Crank-Nicolson法求解限定源扩散问题的误差
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
x_ = zeros(3 / h, 3 / h);	%记录x值数据
y_ = zeros(3 / h, 3 / h);	%记录t值数据
u_ = zeros(3 / h, 3 / h);	%记录u值数据

% 记录t=0.1时的数据
for j = 1:3 / h
    x_(1,j) = j*h;
    y_(1,j) = k;
    u_(1,j) = w(50000+j);
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
        x_(i+1,j) = j*h;
        y_(i+1,j) = t;
        u_(i+1,j) = w(50000+j);
    end
end

% 计算解析解数值的矩阵
x=(0.1:0.1:3);
y=(0.1:0.1:3);
[X,Y] = meshgrid(x,y);
Z=1./sqrt(pi.*Y).*exp(-X.^2./(4.*Y));

% 计算绝对误差
Z1 = abs(Z-u_);

% 计算相对误差
Z2 = abs(Z-u_)./Z;

% 作图
subplot(2,1,1),mesh(X,Y,Z1);
xlabel('x');
ylabel('t');
zlabel('Absolute Error');
subplot(2,1,2),mesh(X,Y,Z2);
xlabel('x');
ylabel('t');
zlabel('Relative Error');