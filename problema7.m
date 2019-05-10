clear all
close all


% Problema 8
% Risoluzione della PDE
%                       |u_t + au_x=0
%                       |
%
%                       |u(x,0)= 1.5*max(0,1-abs(x));
%mediante Lax Wendroff

Mx = 300; %nodi in x
CFL = 0.95; % Courant numero
a = 1.5;%a positivo

tend = 1.1;%tempo finale

u0 = @(x) 1.5*max(0,1-abs(x));
x = linspace(-2,2,Mx+1);
h = 4/Mx;%passo spaziale
k = CFL*h/abs(a); %passo temporale

t = 0;
U = zeros(1,Mx+1);
U = [u0(x)];

while (t+k)<tend
    % Lax-Wendroff
    %Utemp(1,j) = U(1,j) - 0.5*k/h*a*(U(1,j+1)-U(1,j-1)) + 0.5*(k/h*a)^2*(U(1,j+1)-2*U(1,j)+U(1,j-1));
    U(1,2:end-1)= U(1,2:end-1)-0.5*(k/h)*a*(U(1,3:end)-U(1,1:end-2))+...
        0.5*(((k*a)/h)^2)*(U(1,3:end)-2*U(1,2:end-1)+U(1,1:end-2));
    % Nodi estremi: considero il valore adiacente.
    U(:,1)    = U(:,2);
    U(:,Mx+1) = U(:,Mx);
    % Aggiorno
    t = t + k;
    %Area=trapz(x,Utemp(1,:));
    plot(x,u0(x-a*t),'r',x,U,'og')
    xlabel('x')
    ylabel('u(t,x)')
    title(sprintf('t = %0.2f',t))
    legend('sol analitica','Lax Wendroff')
    pause(0.01)
end