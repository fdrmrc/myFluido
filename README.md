# myFluido

ESERCIZIO 1:

3*f''+2f=sin(y)
f(0)=0, f(3*pi/2)=1

Verificare numericamente ordine 2 per FD centrate con distribuzione parabolica dei nodi y_i. 


ESERCIZIO 2:
3f''-2g'+f=0
g''+2f'-2g=-cos(y)
f(0)=0, f(pi/2)=1
g(0)=1, g(pi/2)=0

Verificare ordine 2 (la coppia di soluzioni è (f(y),g(y))=(sin(y),cos(y)) su una distribuzione parabolica dei nodi


ESERCIZIO 3, ESERCIZIO 4,


ESERCIZIO 5:

Risolvere equazione di Blasius (strato limite su lamina piana) in forma conservativa con corrente esterna U costante, nu=1.5 m/s^2 e condizioni al contorno

u(x,0)=v(x,0)=0 
u(x,infty)=U

e condizioni iniziali 

u(0,y)=U, v(0,y)=0


ESERCIZIO 6:

Risolvere equazioni di Blasius scritta come problema ai limiti

fu'+2u''=0
f'-u=0

e condizioni al contorno
f(0)=0, u(0)=0, u(infty)=1

Calcolare fattore di forma H per u(2,y) e confrontarlo con il valore previsto dalla teoria


ESERCIZIO 7:

Risolvere il problema
u_t+au_x =0
u(x,0)=1.5*max(0,1-abs(x))

con un metodo conservativo e confrontarlo con la soluzione ottenuta con il metodo delle caratteristiche


ESERCIZIO 8:

u_t+au_x=0
u(x,0)=1.2 if x<0, 0.4 if x>0

utilizzando un metodo conservativo e confrontarlo con la soluzione ottenuta con il metodo delle caratteristiche. Utilizzando un metodo conservativo, si ottengono risultati più o meno corretti?
Quali conclusioni si possono trarre sulla dipendenza del metodo numerico da eventuali discontinuità del dato iniziale? Come sono spiegabili numericamente?


ESERCIZIO 9:

Risolvere con il metodo delle caratteristiche l'equazione di Burgers

u_t+uu_x=0

con dato iniziale u(x,0)=1.5*max(0,1-abs(x)) fino ad un tempo tstar=1.5*T, dove T è il tempo per il quale si osserva una soluzione a più valori (non fisica).
Quindi si risolva lo stesso problema
numericamente (con un metodo conservativo), ottenendo una soluzione ad un sol valore. 
Si mostri che l’area racchiusa tra l’asse x e la curva viene conservata, ovvero è la stessa per entrambi i metodi.


ESERCIZIO 10:

Risolvere numericamente il problema di Riemann non lineare

u_t + uu_x=0
u(x,0)=u_left if x<0, u_right if x>0

confrontando tra loro le soluzioni ottenute con un metodo non conservativo, uno conservativo e la soluzione esatta
