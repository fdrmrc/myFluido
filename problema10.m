clear all
close all

%% Risoluzione del problema di Riemann non-lineare con dato iniziale discontinuo
%% utilizzando lax-wendroff prima e lax-friedrichs conservativo poi
%% L'utente fornisce l'ampiezza del passo h, e il numero CFL. Il passo k nel tempo viene
%% quindi scelto di conseguenza.

h=0.01; tf=1;
sx=-2;dx=2;
mx=floor((dx-sx)/h);
xx=linspace(sx,dx,mx+1);
u0=@(x) 1.2*(x<0) + 0.4*(x>=0);
S=0.5*(1.2+0.4); %Speed dell'onda d'urto
CFL=0.9;
k=(CFL*h)/max(u0(xx)); %time steps, massimo della velocità iniziale=(tf)/k;
U=zeros(2,mx+1); %prima riga lax wendroff, seconda riga lax friedrichs

U=[u0(xx);u0(xx)];
f=@(u)0.5*u.^2;


    t=0; 
while (t+k)<tf
    %lax-Wendroff.
    U(1,2:end-1)= U(1,2:end-1)-0.5*(k/h)*U(1,2:end-1).*(U(1,3:end)-U(1,1:end-2))+...
        0.5*((k/h)^2)*(U(1,2:end-1).^2).*(U(1,3:end)-2*U(1,2:end-1)+U(1,1:end-2));    
    U(1,1)=U(1,2);U(1,end)=U(1,end-1);
    %lax-Friedrichs
    U(2,2:end-1)=0.5*(U(2,3:end) + U(2,1:end-2)) - (k/(2*h))*(f(U(2,3:end)) - f(U(2,1:end-2)));
    U(2,1)=U(2,2);U(2,mx+1)=U(2,mx);
    k=(CFL*h)/max(U(1,:));
    t=t+k;
    plot(xx,U(1,:),'r-o',xx,U(2,:),'b-o',xx,1.2*(xx-S*t<0)+0.4*(xx-S*t>0),'k')
    xlabel('x')
    ylabel('u(t,x)')
    title(sprintf('Problema di Riemann non-lineare a t=%0.1f',t));
    legend('Lax-Wendroff','Lax-Friedrichs','sol. analitica')
    pause(0.01)
end
