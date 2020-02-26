function U=matlab_summer_3_pde(N)
%����ѧ�Ų��ҳ����ʹδ�ֵ
%H=input('������ѧ�ţ�ÿ�����ֺ��һ���ո�')
M = [3 1 7 0 1 0 4 2 5 7];
M=sort(M);
C=M(8);
A=M(10);
a=max(2,A);
b=max(1,C);

%�ⷽ��
%ϸ��%x����������ʷ�
N=input('N=')
% N=5;
n=N-1;
h=1/n;
%xi=i*h;
%yj=j*h

%AU=F+H
%,��N=51Ϊ��������������49^2*49^2�ķ��󣬼���(n-1)^2*(n-1)^2��*��(n-1)^2*1��=��(n-1)^2*1��+��(n-1)^2*1����ģ

% ��A���Խǿ����
I=eye(n-1);
C=(diag(4*ones(1,n-1))-diag(ones(1,n-2),1)-diag(ones(1,n-2),-1));
A=kron(diag(ones(1,(n-2)),1)+diag(ones(1,(n-2)),-1),I)+kron(diag(ones(1,n-1)),C);
%
%¼��߽�����
U=zeros(n+1,n+1);
U(1,:)=0;
U(:,1)=0;%����߽����������һ��Ԫ����ʵ������������ʱ����ì�ܵģ���Ψһ����
for j=0:n
    U(n+1,j+1)=sin(j*h*b);%����߽�����
    U(j+1,n+1)=sin(j*h*a);%����߽�����
end

%����F
F=zeros((n-1)^2,1);
for i=1:n-1%xi=i*h
    for j=1:n-1%yj=j*h*j
        F(i*j)=sin(a*h*i)*sin(b*h*j);
    end
end

%����H
H=zeros((n-1)^2,1)
for j=1:n-2
    H((n-1)*j)=sin(b*i*h);
end
for j=1:n-2
    H((n-1)*(n-2)+j)=sin(j*h*a);
end
H((n-1)^2)=sin((n-1)*h*a)+sin((n-1)*h*b);
H=H*(1/h^2);

%�õ�΢�ַ�����AX=F+K,��������
X=A\(F+H);
%����F

%��ԭ
for i=1:n-1%��X��ΪN-1��/��
    for j=1:n-1%��ÿһ�鰴�д�ǰ������������U��
        U(j+1,i+1)=X(i*j);
    end
end

% ������
disp('���̵Ľ�Ϊ');

[X,Y]=meshgrid(0:h:1,0:h:1);
Z=U;
mesh(X,Y,Z);%�߿�ͼ
xlabel('x')%x����
ylabel('y')%y����
zlabel('U')%z����
title('��ά���޲��PDE���')%����
figure();
surf(X,Y,Z);%����ͼ
xlabel('x')%x����
ylabel('y')%y����
zlabel('U')%z����
title('��ά���޲��PDE���')%����
% shading interp;
