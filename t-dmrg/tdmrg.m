
function tdmrg(Phi,h_int,Bl,Br,d,m,tf,dt)

global dispR

% Numeros de pasos de tiempo
time_steps = tf/dt;
% Estado inicial:
Psi = Phi; 

% Medir Sz para todos los sitios en el estado Psi0
Szprof = zeros(numel(d),time_steps+1);
Szprof(:,1) = medirSz(Psi,Bl,Br,d);
% Sxprof = zeros(numel(d),time_steps+1);
% Sxprof(:,1) = medirSx(Psi,Bl,Br,d);
SvonN = zeros(time_steps+1,2);
SvonN(1,1:2) = [0,vonNeumannS(Phi,0.5*numel(d),Bl,Br,d)];

if dispR, fprintf('-> t-DMRG, evolucion temporal con (dt= %.3f)\n',dt), end
for n = 1:time_steps
    % Tiempo transcurrido t = n*dt
    t = n*dt;
    % Aplicamos e^(-i*delta*H(t)) a Psi(t-dt) para obtener Psi(t).
    [Psi,Bl,Br,uPm] = evolve(Psi,dt,h_int,Bl,Br,m(end),d);
    % Medir Sz para todos los sitios en el estado Psi(t)
    Szprof(:,n+1) = medirSz(Psi,Bl,Br,d);
    %Sxprof(:,n+1) = medirSx(Psi,Bl,Br,d);
    SvonN(n+1,1:2) = [t,vonNeumannS(Psi,0.5*numel(d),Bl,Br,d)];
    if dispR, fprintf('tiempo = %.3f\n',t), end
    save SzDM Szprof SvonN
end

%========================================================================
function [Psi,Bl,Br,uPm] = evolve(Psi,dt,h_int,Bl,Br,m,d)
% EVOLVE - Ejecuta un paso temporal del algoritmo tdmrg

L = length(d);
% Funcion para identificar los sitios de la red
b = @(l) sitio(l,d);
% Hamiltoniano de interaccion de sitios vecinos x,y
h = @(x,y) feval(h_int,b(x),b(y));

Psi = A1v(expm(-0.5*i*h(1,2)*dt),Psi,Br,d);
Psi = Alv(expm(-0.5*i*h(2,3)*dt),Psi,1,Bl,Br,d);
Bl{2} = sweepR(Bl{1},d(2),d(3),Br{L-3},m,Psi);

for l = 2:L-3
    Psi = chbasis(Psi,l-1,l,Bl,Br,d);
    Psi = Alv(expm(-0.5*i*h(l+1,l+2)*dt),Psi,l,Bl,Br,d);
    Bl{l+1} = sweepR(Bl{l},d(l+1),d(l+2),Br{L-l-2},m,Psi);
end

Psi = ALv(expm(-i*h(L-1,L)*dt),Psi,Bl,d);
Psi = Alv(expm(-0.5*i*h(L-2,L-1)*dt),Psi,L-3,Bl,Br,d);
Br{2} = sweepL(Bl{L-3},d(L-2),d(L-1),Br{1},m,Psi);

for l = L-4:-1:1
    Psi = chbasis(Psi,l+1,l,Bl,Br,d);
    Psi = Alv(expm(-0.5*i*h(l+1,l+2)*dt),Psi,l,Bl,Br,d);
    [Br{L-l-1},uPm] = sweepL(Bl{l},d(l+1),d(l+2),Br{L-l-2},m,Psi);
end

Psi = A1v(expm(-0.5*i*h(1,2)*dt),Psi,Br,d);
%========================================================================
function Av = A1v(A,v,Br,d)
% A1v - Multiplicacion matriz-vector A{1,2} * v

L  = length(d);
dA = d(1)*d(2);
dR = d(3)*Br{L-3}.dim;

w0 = reshape(v,[dR dA])*A.';
Av = reshape(w0,[dA*dR 1]);
%========================================================================
function Av = Alv(A,v,l,Bl,Br,d) 
% Alv - Multiplicacion matriz-vector A{l+1,l+2} * v

L  = length(d);
dA = d(l+1)*d(l+2);
dL = Bl{l}.dim;
dR = Br{L-l-2}.dim;

w0  = reshape(v, [dR dL*dA]).';
w1  = A*reshape(w0, [dA dL*dR]);
w2  = reshape(w1, [dL*dA dR]).';
Av  = reshape(w2, [dL*dA*dR 1]);
%========================================================================
function Av = ALv(A,v,Bl,d)
% ALv - Multiplicacion matriz-vector A{L-1,L} * v

L  = length(d);
dA = d(L-1)*d(L);
dL = Bl{L-3}.dim*d(L-2);

w0 = A*reshape(v,[dA dL]);
Av = reshape(w0,[dL*dA 1]);
%=========================================================================
function b = sitio(l,d)
% SITIO - Determina el sitio etiquetado por l

global b0 b1

if d(l) == b0.dim
    b = b0;
else
    b = b1;
end
%=========================================================================
