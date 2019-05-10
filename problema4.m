clear all
close all

mrange=2.^[3:10];
count=0;
for m=mrange
    count=count+1;
    y=pi*(((1:m)-1)/(m-1))'.^2;
    %y=linspace(0,1,m)';

    h=diff(y);
    d1=1./(h(1:m-2)+h(2:m-1));
    D1a=spdiags([[-d1;0;0],[0;0;d1]],[-1,1],m,m);
    D1b=spdiags([[-d1;0;0],[0;0;d1]],[-1,1],m,m);
    I=speye(m);

    %condizioni al bordo
    %u_1=0
    D1a(1,2)=0;
    %u_m=0
    D1a(m,m-1)=0;

    %v_1=1
    C=sparse(zeros(m,m));
    C(1,1)=1;

    %v_m=1
    C(m,m)=-1;

    %
    F=@(w) [I*w(1:m) + D1a*w(m+1:2*m);...
         diag(w(m+1:2*m))*(D1b*w(1:m))- diag(w(1:m))*(D1b*w(m+1:2*m))-1+C*[w(m+1:2*m)]]; %400x1

    %costruzione dello jacobiano
    r=@(w) diag(w(m+2:2*m-1))*d1;
    s=@(w) -diag(w(m+2:2*m-1))*d1; %g_i/(y_i+1 - y_i-1)
    B3=@(w) spdiags([[s(w);0;0],[0;0;r(w)]],[-1,1],m,m) - diag(D1b*w(m+1:2*m))

    r=@(w) diag(w(2:m-1))*d1; %f_i/(y(i+1)-y(i-1))
    s=@(w) -diag(w(2:m-1))*d1; 
    B4=@(w) -spdiags([[s(w);0;0],[0;0;r(w)]],[-1,1],m,m) + diag(D1b*w(1:m))+C %C deriva da condizioni al bordo

    J=@(w) [I,D1a;
            B3(w),B4(w)]; %400x400

     %Newton

    w0=[0;ones(m-2,1);0;1;ones(m-2,1);-1]; %guess iniziale
    w=w0
    res=-J(w0)\F(w0);
    tol=max(h)^2;
    while norm(res,inf)>tol
        w=w+res;
        res=-J(w)\F(w);
        norm(res,inf);
    end
    w=w+res;
    uesatta=[sin(y);cos(y)];
    errore(count)=norm(uesatta-w,inf)
end
plot(y,w(1:m),'ok',y,w(m+1:2*m),'ok',y,sin(y),'.r',y,cos(y),'.r')

%% grafico errore spaziale
figure;
loglog(mrange,errore,'*',mrange,errore(end)*(mrange/mrange(end)).^(-1),...
                         mrange,errore(end)*(mrange/mrange(end)).^(-2),...
                         mrange,errore(end)*(mrange/mrange(end)).^(-3),...
                         mrange,errore(end)*(mrange/mrange(end)).^(-4))
legend('errore','ordine 1','ordine 2','ordine 3','ordine 4')