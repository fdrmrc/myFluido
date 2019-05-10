clear all
close all

%% Risoluzione dell'equazione di Burgers con metodo delle
%% caratteristiche e Lax - Friedrichs conservativo


tb=2/3; %tb=-1/(min(du0/dx))
tf=1.5*tb;
sx=-10;dx=10;
h=0.01;
mx=floor((dx-sx)/h); %numero nodi 
xx=linspace(sx,dx,mx+1); 
CFL=0.8;
u0=@(x) 1.5*max(0,1-abs(x)); %dato iniziale
k=(CFL*h)/max(u0(xx)); %time steps, massimo della velocità iniziale=(tf)/k;
f=@(u) 0.5*u.^2; %funzione che rappresenta il termine non-lineare dell'equazione di Burgers

U=zeros(1,mx+1);
U=[u0(xx)];

%% Metodo delle caratteristiche

AreaChar=quad(u0,sx,dx); %ovviamente fa 1.5 visto che l'altezza e la base non sono cambiate
%uso quad perché è implementato sia in octave che in MatLab

sprintf('L''area sottesa con caratteristiche è %0.5f',AreaChar)

%% Lax-Friedrichs conservativo
t=0;
i=0; %contatore tempo, per metodo caratteristiche
T=0:k:tf; %linspace dei tempi
while (t+k)<tf
  %caratteristiche
    i=i+1;
    %Lax-Friedrichs conservativo
    U(1,2:end-1)=0.5*(U(1,3:end) + U(1,1:end-2)) - (k/(2*h))*(f(U(1,3:end)) - f(U(1,1:end-2)));
    U(1,1)=U(1,2); %nodi interni
    U(1,mx+1)=U(1,mx);
    k=(CFL*h)/(max(U(1,:))); %il k va cambiato poiché la condizione su CFL va modificata
    t=t+k; %aggiorno il tempo e la soluzione
    AreaLax=trapz([sx:h:dx],U(1,:))
    plot(xx,U(1,:),'r-o',xx+T(i)*u0(xx),u0(xx),'b',ones(2),[0,1.5],'r--')
    xlabel('x')
    ylabel('u(t,x)')
    legend('Sol. numerica','Caratteristiche','tempo sol multivalued')
    title(sprintf('Lax-Friedrichs L''area sottesa è %0.8f per t=%0.2f',AreaLax,t));
    pause(0.0001)
end