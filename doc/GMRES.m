% main
clear; clc; close all 
A = pascal(4);
b = [0 0 0 0]';
x0 = [1 0 0 0]';
%
[V,R,H,res] = bGMRES(A,b,x0)
%
r0 = b-A*x0;beta = norm(r0);
be = zeros(4,1);be(1) = beta;
be
[T,bk] = givens( H,be )
%
newT = zeros(4,4);
for i = 1:4
    for j = 1:i
        newT(j,i) = T(j,i);
    end
end
newT
x4 = inv(newT)*bk
x5 = backward( newT,bk )
%
realSolution = inv(A)*b
sol1 = x0+V(:,1:4)*x4
sol2 = x0+V(:,1:4)*x5

clear; clc; close all 
A = sprandsym(10,0.7);
A
b = [0 0 0 0 0 0 0 0 0 0]';
x0 = [1 0 0 0 0 0 0 0 0 0]';
%
[V,R,H,res] = bGMRES(A,b,x0)
%
r0 = b-A*x0;beta = norm(r0);
be = zeros(10,1);be(1) = beta;
be
[T,bk] = givens( H,be )
%
newT = zeros(10,10);
for i = 1:10
    for j = 1:i
        newT(j,i) = T(j,i);
    end
end
newT
x4 = inv(newT)*bk
x5 = backward( newT,bk )
%
realSolution = inv(A)*b
sol1 = x0+V(:,1:10)*x4
sol2 = x0+V(:,1:10)*x5

%Ax = b
function [V,R,H,res] = bGMRES(A,b,x0)
    %bGMRES:basic GMRES method
    %Input: x0:初值;A为mxm矩阵,b为解
    %Output: res为残差
    [m, ~] = size(A);
    R = Inf(m,m);%R为剩余向量
    H = zeros(m+1,m);V = zeros(m,m+1);%A*V=V*H
    %设定初值
    r0 = b-A*x0;
    V(:,1) = r0./norm(r0);
    for j = 1:m
        R(:,j) = A*V(:,j);
        for i = 1:j
            H(i,j) = R(:,j)'*V(:,i);
            R(:,j) = R(:,j) - H(i,j)*V(:,i);
        end
        H(j+1,j) = norm(R(:,j));
        res = H(j+1,j);
        if abs(H(j+1,j)) < 1e-10
            sprintf('done without residual')
            break;
        else
            V(:,j+1) = R(:,j)./H(j+1,j);
        end
    end
end

function [T,bk] = givens( H,b )
%givens: 通过givens变换化上Hessenborg阵为上三角矩阵
%  化 Hx=b 为 Tx=c
    [~,n] = size(H);
    %提取
    Ht = H(n+1,:);
    H = H(1:n,:);
    b = b(1:n);
    %Rotate Matrix need to recrate every iteration
    %R = eye(n,n);%Rotate Matrix
    for k = 1:n-1
        R = eye(n,n);
        down = (H(k,k)^2+H(k+1,k)^2)^(1/2);
        s = H(k+1,k)/down;
        c = H(k,k)/down;
        R(k:k+1,k:k+1) = [c,s;-s,c];
        R
        H = R*H
        H(k+1,k) = 0;
        b = R*b
    end
    T = [H;Ht];
    bk = b;
end

function x = backward( A,b )
%backward:   上三角矩阵回带过程
% 求解Ax=b
    [~,c] = size(A);
    A = A(1:c,:);
    x = zeros(c,1);
    x(end) = b(end)/A(c,c);
    for k = 2:c
        V = x(c-k+2:c);
        x(c-k+1) = (b(c-k+1)-A(c-k+1,c-k+2:end)*V)/A(c-k+1,c-k+1);
    end
end
