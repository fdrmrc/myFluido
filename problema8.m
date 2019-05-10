clear all
close all

% Problema 8
% Risoluzione della PDE
%                       |u_t + au_x=0
%                       |
%                                   1.2 se x<0
%                       |u(x,0)=
%                       |           0.4 se x>0
%mediante Lax Wendroff (soluzione numerica non buona) e Lax Friedrichs conservativo.
% CFL= numero caratteristico del problema

% L'utente definisce h, CFL e k viene determinato di conseguenza.


h=0.01;
ax=-6; bx=6;
Mx = (bx-ax)/h; %nodi in x. Mx=dx-sx/h;
x = linspace(ax,bx,Mx+1);
CFL = 0.95; % Courant number
a = -0.5; %a positivo
tend = 2; %tempo finale
u0 = @(x) 1.5*(x<=0) + 0.5*(x>0);

k = (CFL*h)/abs(a); %passo temporale
t = 0;
U = zeros(2,Mx+1);
U = [u0(x); u0(x)];


while (t+k)<tend    
    %Lax Wendroff155
    U(1,2:end-1)= U(1,2:end-1)-0.5*(k/h)*a*(U(1,3:end)-U(1,1:end-2))+...
        0.5*((k/h)^2)*(a^2).*(U(1,3:end)-2*U(1,2:end-1)+U(1,1:end-2));    
    U(1,1)=U(1,2);U(1,end)=U(1,end-1);
    %lax-Friedrichs
    U(2,2:end-1)=0.5*(U(2,3:end) + U(2,1:end-2)) - (k/(2*h))*a*(U(2,3:end) - U(2,1:end-2));
    U(2,1)=U(2,2);U(2,end)=U(2,end-1);
    
    % Aggiorno
    t = t + k;
    plot(x,U(1,:),'b-o',x,U(2,:),'r-o',x,1.5*(x-a*t<0)+0.5*(x-a*t>=0),'k')
    title(sprintf('Soluzioni a t = %0.2f',t));
    legend('Lax Wendroff','Lax Friedrichs','sol analitica')
    axis([-6,6,0.3,1.6])
    pause(0.01)
end
figure  %zoom sui due metodi
plot(x,U(1,:),'b-o',x,U(2,:),'r-o',x,1.5*(x-a*t<0)+0.5*(x-a*t>=0),'k')
title('Soluzioni al tempo finale')
xlabel('x') 
ylabel('sol')
legend('LW','LF','Sol analitica')
axis([-1.2,0.6,0.3,1.5])
