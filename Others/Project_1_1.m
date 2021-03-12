clear;clc;
L=10000; T=3; m=100000; n=30;
h=L/m; k=T/n; a=1;
lambda = a * a * k / (h * h);
l=zeros(m-1,1);
u=zeros(m-1,1);
z=zeros(m-1,1);
w=zeros(m-1,1);
t=zeros(m-1,1);
x_ = zeros(3 / h * 3 / h,1);
y_ = zeros(3 / h * 3 / h,1);
z_ = zeros(3 / h * 3 / h,1);
for i = 1:m - 1
    w(i) = -1;
end

l(1) = 1 + 2 * lambda;
u(1) = -lambda / l(1);

for i = 2:m - 2
    l(i) = 1 + 2 * lambda + lambda * u(i - 1);
    u(i) = -lambda / l(i);
end

l(m - 1) = 1 + 2 * lambda + lambda * u(m - 2);
for i = 1:n 
    t = i * k;
    z(1) = w(1) / l(1);
    for j = 2:m-1
        z(j) = (w(j)+lambda*z(j-1))/l(j);
    end
    w(m-1) = z(m-1);
    for j = m-2:-1:1
        w(j) = z(j)-u(j)*w(j+1);
    end
    %printf("t = %lf\n", t);
    %for j = 1:m-1
        %x(j) = j*h;
        %if t==3
        %    break;
        %    fprintf("u(x = %.3f, t = %.3f) = %.6f\n",x,t,w(j)+1);
        %    if t==3 && j>300
        %        break;
        %    end
        %end
    %end
    for j = 1:3 / h
        x_(j+(i-1)*3 / h) = j*h;
        y_(j+(i-1)*3 / h) = t;
        z_(j+(i-1)*3 / h) = w(j)+1;
    end

end

[X,T]=meshgrid(linspace(min(x_),max(x_)),linspace(min(y_),max(y_)));
U=griddata(x_,y_,z_,X,T,'v4');
mesh(X,T,U);
xlabel('X');
ylabel('T');
zlabel('U');


