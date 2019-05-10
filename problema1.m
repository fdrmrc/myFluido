clear all
close all

% Risoluzione del problema ai limiti 1. 
% Approssimazione della derivata seconda mediante FD non equispaziate. Ordine è O(max h_i^2, h_i-1^2)

sol=@(y) -sin(y); %soluzione analitica
mrange=2.^(4:9) +1;
count=0;

for m=mrange 
    count+=1;
    for n=1:m
         x(n,1)=1.5*pi*((n-1)/(m-1))^2; %genera distribuzione di nodi di tipo parabolico
    end
    h=diff(x);
    % Costruzione matrice discretizzante
    v=-1./(h(1:m-2).*h(2:m-1)); %diag princ 
    w=1./((h(1:m-2)).*(h(1:m-2) + h(2:m-1))); %sotto diag
    s=1./((h(2:m-1)).*(h(1:m-2) + h(2:m-1))); %sopra diag
 
    D2=2*spdiags([[w;0;0],[0;v;0],[0;0;s]],[-1,0,1],m,m);
    I=speye(m);
    %Funzionale da azzerare è lineare: non è necessario utilizzare il metodo di Newton.

    A=3*D2 + 2*I; %jacobiano lineare

    %%condizioni al bordo & risoluzione sistema lineare
    A(1,1:2)=[1,0];
    A(m,m-1:m)=[0,1];
    b=[sin(x(1:m-1));1];
    u=A\(b);

    err(count)=norm(u-sol(x),inf);
end

figure
plot(x,u,'*',x,sol(x),'g')
title('Risoluzione numerica')
xlabel('x')
ylabel('f(x)')
legend('numerical','analytical')
figure
loglog(mrange,err,'*',mrange,err(end) * (mrange/mrange(end)).^(-2),'g')
title('Errore in norma infinito')
