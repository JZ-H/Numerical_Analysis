% coding: utf-8
% description: 计算限定源扩散问题的解析解数值
clear;clc;
x=(0.1:0.05:3);
y=(0.1:0.05:3);
[X,Y] = meshgrid(x,y);
% 计算u值的矩阵
Z=1./sqrt(pi.*Y).*exp(-X.^2./(4.*Y));
% 作图
mesh(X,Y,Z);
xlabel('x');
ylabel('t');
zlabel('u');