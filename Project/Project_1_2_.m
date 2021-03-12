clear;clc;
% 相关参数的初始化
L=30000; T=3; m=1000000; n=100;
h=L/m; k=T/n; a=1;
lambda = a * a * k / (h * h);
l=zeros(m-1,1);
u=zeros(m-1,1);
z=zeros(m-1,1);
w=zeros(m,1);
w(m) = -1;
for i = 1:m - 1 
    w(i) = -1;
end

% 初始化用于存储结果数据的向量
x_ = zeros(3 / h, 3 / h);	%记录x值数据
y_ = zeros(3 / h, 3 / h);	%记录t值数据
u_ = zeros(3 / h, 3 / h);	%记录u值数据

% 通过线性系统迭代求解
l(1) = 1 + lambda;
u(1) = -lambda/(2*l(1));
for i = 2:m-2
    l(i) = 1 + lambda + lambda*u(i-1)/2;
    u(i) = -lambda/(2*l(i));
end
l(m-1) = 1 + lambda + lambda*u(m-2)/2;
for i = 1:n
    t = i*k;
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
        x_(i,j) = j*h;
        y_(i,j) = t;
        u_(i,j) = w(j)+1;
    end
end

% 作图;
mesh(x_,y_,u_);
xlabel('x');
ylabel('t');
zlabel('u');