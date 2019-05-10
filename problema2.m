%Es 2

clear all
close all

%3*f''-2*g'+f=0
%g''+2*f'-2*g=-cosy
%f(0)=0
%f(pi/2)=1
%g(0)=1
%g(pi/2)=0

mrange=2.^(4:10);
count=0;

for m=mrange
count=count+1;
y=pi/2*(((1:m)-1)/(m-1))'.^2;


h=diff(y);
d1=1./(h(1:m-2)+h(2:m-1));
D1=spdiags([[-d1;0;0],[0;0;d1]],[-1,1],m,m);

v1=2./h(1:m-2).*d1;
v2=-2./(h(2:m-1).*h(1:m-2));
v3=2./h(2:m-1).*d1;

D2=spdiags([[v1;0;0],[0;v2;0],[0;0;v3]],[-1,0,1],m,m);
I=speye(m);

Z=spdiags(zeros(m,1),0,m,m);

F1=3*D2+I;
F2=-2*D1;
F3=2*D1;
F4=D2-2*I;
F=[F1,F2;
   F3,F4];



b=[zeros(m,1);-cos(y)];
F(1,1)=1;
F(m,m)=1;
F(m+1,m+1)=1;
F(2*m,2*m)=1;
b(m)=1;
b(m+1)=1;
b(2*m)=0;
r=F\b;

Sol=[sin(y);cos(y)];

err(count)=norm(r-Sol,inf);

end


figure
plot(y,r(1:m),'or',y,r(m+1:2*m),'ok',y,sin(y),'r',y,cos(y),'b')
legend('f(x)','g(x)','sin(x)','cos(x)')
title('Soluzioni numeriche e analitiche')
figure
loglog(mrange,err,'*',mrange,mrange(1)*(mrange/mrange(1)).^(-2),'r')
title('Errore in norma infinito')