function U=matlab_summer_3_pde(N)
%输入学号并找出最大和次大值
%H=input('请输入学号，每个数字后加一个空格')
M = [3 1 7 0 1 0 4 2 5 7];
M=sort(M);
C=M(8);
A=M(10);
a=max(2,A);
b=max(1,C);

%解方程
%细度%x对区间进行剖分
N=input('N=')
% N=5;
n=N-1;
h=1/n;
%xi=i*h;
%yj=j*h

%AU=F+H
%,以N=51为例最后做运算的是49^2*49^2的方阵，即（(n-1)^2*(n-1)^2）*（(n-1)^2*1）=（(n-1)^2*1）+（(n-1)^2*1）规模

% 求A三对角块矩阵
I=eye(n-1);
C=(diag(4*ones(1,n-1))-diag(ones(1,n-2),1)-diag(ones(1,n-2),-1));
A=kron(diag(ones(1,(n-2)),1)+diag(ones(1,(n-2)),-1),I)+kron(diag(ones(1,n-1)),C);
%
%录入边界条件
U=zeros(n+1,n+1);
U(1,:)=0;
U(:,1)=0;%带入边界条件（最后一个元素其实在条件给定的时候是矛盾的（不唯一））
for j=0:n
    U(n+1,j+1)=sin(j*h*b);%带入边界条件
    U(j+1,n+1)=sin(j*h*a);%带入边界条件
end

%构造F
F=zeros((n-1)^2,1);
for i=1:n-1%xi=i*h
    for j=1:n-1%yj=j*h*j
        F(i*j)=sin(a*h*i)*sin(b*h*j);
    end
end

%构造H
H=zeros((n-1)^2,1)
for j=1:n-2
    H((n-1)*j)=sin(b*i*h);
end
for j=1:n-2
    H((n-1)*(n-2)+j)=sin(j*h*a);
end
H((n-1)^2)=sin((n-1)*h*a)+sin((n-1)*h*b);
H=H*(1/h^2);

%得到微分方程组AX=F+K,用左除求解
X=A\(F+H);
%构造F

%还原
for i=1:n-1%把X分为N-1个/组
    for j=1:n-1%对每一组按列从前到后塞到矩阵U中
        U(j+1,i+1)=X(i*j);
    end
end

% 输出结果
disp('方程的解为');

[X,Y]=meshgrid(0:h:1,0:h:1);
Z=U;
mesh(X,Y,Z);%线框图
xlabel('x')%x轴标记
ylabel('y')%y轴标记
zlabel('U')%z轴标记
title('三维有限差分PDE求解')%标题
figure();
surf(X,Y,Z);%表面图
xlabel('x')%x轴标记
ylabel('y')%y轴标记
zlabel('U')%z轴标记
title('三维有限差分PDE求解')%标题
% shading interp;
