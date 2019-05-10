clear all
%close all

M_x = 800; % nodi in x
M_y = 200; % nodi in y
U = 3; % parametri
nu = 1.5;
ay = 0; % estr sinistro in y
by = 10; % estr destro in y
ax = 0;% estr sinistro in x
bx = 2;% estr destro in x
y = (by-ay)*(((1:M_y)-1)/(M_y-1))'.^2;
h = diff(y);
x = (bx-ax)*(((1:M_x)-1)/(M_x-1))'.^2;
k = diff(x);



maxit = 50; %numero massimo di iterazioni
% dati iniziali
u0 = U*ones(M_y,1); % u(0,y)=U
u0(1) = 0; %no-splip condition (primo punto della parete)
v0 = zeros(M_y,1); %prima della lamina la corrente è parallela alla direz x
w0 = [u0;v0]; % argomento da passare alla funzione ([u1,u2,..,un,v1,v2,...vn])
uOld = u0; % come guess iniziali di Newton
vOld = v0; 
tic 
for i = 1:length(x)-1 % arrivo a M_y
    % prima equazione sistema
    delta_x = x(i+1) - x(i); % k(i); %x(2)-x(1) = k(1) = secondo punto - 0 alla prima iterata
    C = (0.5/delta_x) * spdiags([ones(M_y,1), ones(M_y,1)], [0,-1], M_y,M_y); % moltiplica u_i
    delta_y = h ; % vettore diff(y) rinominato per comodit�
    D = (-0.5/delta_x) * spdiags([ones(M_y,1), ones(M_y,1)], [0,-1], M_y,M_y); % moltiplica uOld
    E = spdiags([[0; 1./(y(2:end-1)-y(1:end-2)); 0], [-1./(y(2:end-1)-y(1:end-2)); 0; 0]], [0,-1], M_y,M_y); % moltiplica v_i
    % seconda equazione
    F = (1/delta_x) * speye(M_y); % moltiplica (u_i)^2
    G = spdiags([[0;0;1./(y(3:end)-y(1:end-2))], [-1./(y(3:end)-y(1:end-2));0;0]],[1,-1],M_y,M_y); % moltiplica u.*v
    H = spdiags([[0;-2./((y(3:end)-y(1:end-2)).*((y(3:end)-y(2:end-1))));0],...
        [0;0;2./((y(3:end)-y(1:end-2)).*((y(3:end)-y(2:end-1))))]],[0,1],M_y,M_y); % moltiplica la prima parte u_i
    I = spdiags([ [0;-2./((y(3:end)-y(1:end-2)).*((y(2:end-1)-y(1:end-2))));0],...
        [2./((y(3:end)-y(1:end-2)).*((y(2:end-1)-y(1:end-2))));0;0]], [0,-1], M_y,M_y); % moltiplica la seconda parte u_i
    J = -(1/delta_x) * speye(M_y); % moltiplica uOld
    K = sparse(M_y,M_y); % aggiungo per condizioni al bordo su v
    
    % CONDIZIONI AI BORDI
    % u(x,0) = 0
    C(1,1) = 1;
    C(1,2:M_y) = zeros(1,M_y-1);
    D(1,:) = zeros(1,M_y);
    E(1,:) = zeros(1,M_y);
    % u(x,inf) = U
    C(M_y,M_y) = 1;
    C(M_y,1:M_y-1) = zeros(1,M_y-1);
    D(M_y,:) = zeros(1,M_y);
    E(M_y,:) = zeros(1,M_y);
    noto = [zeros(M_y-1,1);-U];
    % v(x,0) = 0 v'(x,inf) = 0 DERIVATA NULLA perchè lontano dalla lamina
    % non varia più effetto di lamina non si sente su v
    K(1,1) = 1;
    F(1,:) = zeros(1,M_y);
    G(1,:) = zeros(1,M_y);
    H(1,:) = zeros(1,M_y);
    I(1,:) = zeros(1,M_y);
    J(1,:) = zeros(1,M_y);
    K(M_y,M_y-1:M_y) = [-1,1]/(y(M_y)-y(M_y-1)); %condizione sulla derivata prima
    F(M_y,:) = zeros(1,M_y);
    G(M_y,:) = zeros(1,M_y);
    H(M_y,:) = zeros(1,M_y);
    I(M_y,:) = zeros(1,M_y);
    J(M_y,:) = zeros(1,M_y);
    
    fun = @(u,v) [C*u + D*uOld + E*v + noto;
        F*u.^2 + J*uOld.^2 + G*(u.*v) - nu*(H+I)*u + K*v];
    % i meno sono gia' inizializzati nella matrice
    % metto K2*ones cosi nell'ultima riga della seconda equazione do la
    % condizione v_n = vOld_n (v al nodo end vale quanto vOld, al passo
    % precedente, al nodo end), cioè praticamente do condizione iniziale
    % diversa per v ogni volta che la calcolo
    
    %jacobiano
    jfun = @(u,v) [               C              ,       E       ;
        2*F*diag(u) + G*diag(v) - nu*(H+I),  G*diag(u) + K  ];
    
    tol = max(h)^2;
    % guess iniziale
    w0 = [u0;v0];
    j0 = jfun(w0(1:M_y),w0(M_y+1:2*M_y));
    res = -j0\fun(w0(1:M_y),w0(M_y+1:2*M_y));
    iter = 0;
    while(norm(res,inf)>tol)
        w0 = w0+res;
        iter = iter+1;
        res = -jfun(w0(1:M_y),w0(M_y+1:2*M_y))\fun(w0(1:M_y),w0(M_y+1:2*M_y));
        norm(res,inf); % controllo decada quadraticamente
    end
    w0 = w0+res;
    
    uOld = w0(1:M_y);  
end
xx = y; % plot su nodi non equisp
plot(w0(1:M_y),xx,'b',w0(M_y+1:2*M_y),xx,'r')
xlabel('soluzioni')
ylabel('y')
toc

%% Fattore forma, coerente con la teorai
delta = trapz(y,1-w0(1:M_y)/U);
theta = trapz(y,(w0(1:M_y)/U).*(1-w0(1:M_y)/U));
H = delta/theta