clear all
%close all

%% Problema 6
% f*u' + 2*u'' = 0
% f' - u = 0
% f(0) = 0        u(0) = 0        u(inf) = 1
% f = f(eta) 
% u = f'

N = 5000; 
aa = 0;
bb = 8;
eta_i = (bb-aa) * (((1:N)-1)/(N-1)).^2;
eta_i = eta_i(:);
h = diff(eta_i);

    
%% discretizzazione 
d1 = 1./(h(1:N-2)+h(2:N-1)); % sopradiagonale D1 = - sottodiagonale
D1 = spdiags([[-d1;0;0],[0;0;d1]],[-1,1],N,N); % già zeri in prima e ultima riga

sopra_diag = [0; 0; 2./(h(2:end).*(h(1:end-1)+h(2:end)))]; 
sotto_diag = [2./(h(1:end-1).*(h(1:end-1)+h(2:end))); 0; 0]; 
diago = [0; -2./(h(1:end-1).*h(2:end)); 0]; 
D2 = spdiags([sotto_diag, diago, sopra_diag], [-1,0,1], N,N);


%% sistema
C = zeros(N,N);
C(1,1) = 1;
alpha = h(N-1)/(h(N-2)*(h(N-1)+h(N-2)));
beta = -(h(N-1)+h(N-2))/(h(N-1)*h(N-2)); 
gamma = - alpha - beta; % combinazione lineare per differenze non centrate non equisp all'indietro ordine 2
C(N,N-2:N) = [alpha, beta, gamma];
c = [zeros(N-1,1);1]; % per u(infinito) con + c, per f(infinito) con -c
F = @(f,u) [diag(f)*(D1*u) + 2*D2*u + C*f - c; % diag([f(1:N-1);h(N-2)+h(N-1)+f(N-1)])
    D1*f - u + c];

I = eye(N,N);
I(N,N) = 0;
I(N,N-1) = 1;
JF = @(f,u) [diag(D1*u) + C, diag(f)*D1 + 2*D2;
    D1, -eye(N,N)];

tol = max(h)^2;
% u0 = [zeros(N-1,1); h(N-2)+h(N-1); zeros(N-1,1); 1];
u0 = [0;ones(N-1,1);0;ones(N-1,1)];
delta = -JF(u0(1:N),u0(N+1:end))\F(u0(1:N),u0(N+1:end));
while norm(delta,inf)>tol
    u0 = u0 + delta;
    delta = -JF(u0(1:N),u0(N+1:end))\F(u0(1:N),u0(N+1:end));
end
u0 = u0 + delta;

f = u0(1:N);
fSecondo=D2*f;
figure
plot(u0(1:N),eta_i,'b',u0(N+1:end),eta_i,'r',fSecondo(10:end),eta_i(10:end),'k')
axis([0,6.5,0,8])
title('soluzioni')
ylabel('eta')
legend('f','f''','f''''')

    

%% u(2,y) e v(2,y)
U = 3;
nu = 1.5;
u_sol = U*(D1*f); % formula u(x,y) con Blasuis, u(2,y)beta
v_sol = 1.5*((0.5*eta_i).*(D1*f) - 0.5*f); % formula v(2,y) con Blasuis


figure
plot(u_sol(1:end-1),eta_i(1:end-1),'b-o',v_sol(1:end-1),eta_i(1:end-1),'r-o')
ylabel('eta')
legend('u','v')
% l'ultima componenete è stata tolta poiché la matrice D1 ha come
%ultima riga zeri e quindi l'ultimo valore delle soluzioni risulta negativo.


%% H fattore forma per u(2,y)
spessore_spostamento_numerico = trapeziNonEq(eta_i,1-(u_sol/U));
spessore_quantita_moto_numerico = trapeziNonEq( eta_i,u_sol/U .* (1-u_sol/U));
spessore_spostamento_teorico = 1.7208*(sqrt(nu*2/U));
spessore_quantita_moto_teorico = 0.664*sqrt(nu*2/U);
H_teorico = spessore_spostamento_teorico/spessore_quantita_moto_teorico;
H_numerico = spessore_spostamento_numerico/spessore_quantita_moto_numerico;
disp(sprintf('Il fattore di forma numerico vale %d',H_numerico));
disp(sprintf('Il fattore di forma teorico vale %d',H_teorico));