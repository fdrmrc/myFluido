clear all
close all

 m=200;
 b=pi/2;
for n=1:m
         x(n,1)=b*((n-1)/(m-1))^2; %genera distribuzione di nodi di tipo parabolico
end
h=diff(x);

% Costruzione matrice discretizzante D2
v=-1./(h(1:m-2).*h(2:m-1)); %diag princ 
w=1./((h(1:m-2)).*(h(1:m-2) + h(2:m-1))); %sotto diag
s=1./((h(2:m-1)).*(h(1:m-2) + h(2:m-1))); %sopra diag
D2=2*spdiags([[w;0;0],[0;v;0],[0;0;s]],[-1,0,1],m,m);
%Costruzione D1
d=1./(h(1:m-2)+h(2:m-1));
D1=spdiags([[-d;0;0],[0;0;d]],[-1,1],m,m);
%Impostazione delle BD
D2(1,1:2)=[2,0];
D1(1,1:2)=[0,0];
D2(m,m-1:m)=[0,2];
D1(m,m-1:m)=[0,0];
%Funzionale da azzerare:
b= [sin(0.5*x(1:m-1)).*cos(1.5*x(1:m-1));1]; %termine noto
F=@(u) 0.5*(D2*u) + diag(u)*(D1*u) - b;
J=@(u) 0.5*D2 + diag(D1*u)+diag(u)*D1; %jacobiano calcolato analiticamente


u=ones(m,1);
tol=max(h)^2;
res=-J(u)\F(u);iter=0;
while norm(res,inf)>tol
    iter+=1;
    u+=res;
    res=-J(u)\F(u);
end
u+=res;

figure
plot(x,sin(x),'g',x,u,'r-x')
xlabel('x')
ylabel('f(x)')
legend('analytical','numerical')
title('Risoluzione numerica BVP non lineare')