clear; clc;
m = 6;
n = 6;
N = m + n;
a = zeros(m+1,1);
q = zeros(m+1,1);
p = zeros(n+1,1);
b = zeros(N,N+1);

% 将ai初始化为sin(x)麦克劳林展开式的系数
for i = 1:N+1
    a(i) = ((-1)^(ceil(i/2)+1))*(mod(i,2)*sin(0)+mod(i-1,2)*cos(0))/factorial(i-1);
end
q(1) = 1;
p(1) = a(1);

% 构建方程组
for i = 1:N
    for j = 1:i-1
        if j<=n
            b(i,j) = 0;
        end
    end
    if i<= n
        b(i,i) = 1;
    end
    for j = i+1:N
        b(i,j) = 0;
    end
    for j = 1:i
        if j<=m
            b(i,n+j)=-a(i-j+1);
        end
    end
    for j = n+i+1:N
        b(i,j) = 0;
    end
    b(i,N+1) = a(i+1);
end

% 求解方程组
for i = n+1:N-1
    % 确定主元
    k = i;
    for j = i+1:N
        if abs(b(j,i)) > abs(b(k,i))
            k = j;
        end
    end
    if b(k,i) == 0
        fprintf("The system is singular")
        return;
    end
    if k ~= i	%k不为i时进行交换
        for j = i:N+1
            temp = b(i,j);
            b(i,j) = b(k,j);
            b(k,j) = temp;
        end
    end
    % 消元
    for j = i+1:N
        xm = b(j,i)/b(i,i);
        for k = i+1:N+1
            b(j,k) = b(j,k) - xm*b(i,k);
        end
        b(j,i) = 0;
    end
end
if b(N,N) == 0
    fprintf("The system is singular")
    return;
end

% 求解出qi、pi
if m>0
    q(m+1) = b(N,N+1)/b(N,N);
end
for i = N-1:-1:n+1
    sum = 0;
    for j = i+1:N
        sum = sum + b(i,j)*q(j-n+1);
    end
    q(i-n+1) = (b(i,N+1)-sum)/b(i,i);
end
for i = n:-1:1
    sum = 0;
    for j = n+1:N
        sum = sum + b(i,j)*q(j-n+1);
    end
    p(i+1) = b(i,N+1)-sum;
end

% 输出结果
for i = 1:n+1
    fprintf("p%d = %.6f  ",i-1,p(i))
end
fprintf("\n")
for i = 1:n+1
    fprintf("q%d = %.6f  ",i-1,q(i))
end
fprintf("\n")