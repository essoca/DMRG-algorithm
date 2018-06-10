
function main()

% Parametros del sistema
global b0 b1 
S1      = 0.5;
S0      = 0.5;
celdaU  = [1 0];
Nrp     = 16;
% Parametros del modelo
global Jxy Jz
h_int = 'hXXZ';
Jxy   = 1;
Jz    = 1;
% Parametros del DMRG
global Ne dispR Msamp Lt
err     =  1e-6;
Msamp   = 100;
Lt      = 50;
Ne      = 1;
dispR   = 1;

clc
% Sitios b0 y b1 iniciales
[Sz0 Sp0] = matrices_spin(S0);
[Sz1 Sp1] = matrices_spin(S1);
b0        = bloque(1,2*S0+1,zeros(2*S0+1),Sp0,Sz0,eye(2*S0+1));
b1        = bloque(1,2*S1+1,zeros(2*S1+1),Sp1,Sz1,eye(2*S1+1));

% Bloques iniciales 
[Bl,Br,d,L] = bloques_ini(S0,S1,celdaU,Nrp);

% Funcion para identificar los sitios de la red
b = @(l) sitio(l,d);
% Energia de 2 sitios
[Phi2,E2] = eig(full(feval(h_int,b(1),b(L))));
E_ = E2(1);

% Temperatura
T = 0.1:0.1:2; beta = 1./T; ET = zeros(length(T),2);
for iT = 1:length(T)
    % Algoritmo infinito
    if dispR, fprintf('-> Algoritmo Infinito\n'), end
    for l = 1:L/2-1
        fprintf('==========================')
        fprintf(' LBL = %d, LBR = %d\n',l+1,l+1)
        [Bl{l+1},Br{l+1},Emin,Emax,uPm] = ...
            dmrg_infinito(Bl{l},b(l+1),b(L-l),Br{l},err,h_int,beta(iT));
        fprintf('E0 = %.12f\n',0.5*(Emin-E_)); E_ = Emin;
        fprintf('Emin = %.12f, Emax = %.12f\n',Emin,Emax)
        fprintf('err = %.0E, m = %d\n',err,uPm)
    end
    % Bloques del superbloque
    Bll = Bl{L/2-1}; bl = b(L/2);
    Brl = Br{L/2-1}; br = b(L/2+1);
    % Dimension del superbloque
    dimSB = Bll.dim*d(L/2)*d(L/2+1)*Brl.dim;
    % Alargar bloques
    Ble = alargar(Bll,bl,'haciaR',h_int);
    Bre = alargar(Brl,br,'haciaL',h_int);
    Z = 0; U = 0;
    % Random sampling and averaging
    fprintf('-> Random sampling and averaging (M = %d)...\n',Msamp)
    for k = 1:Msamp
        % Crear xi aleatorio
        xi = 2*rand(dimSB,1)-1; xi = xi/norm(xi);
        % Obtener el target state xim
        xim = getxim(xi,Emin,Emax,Ble,Bre,h_int,beta(iT));
        % Medir valores esperados
        Z = Z + norm(xim)^2;
        U = U + real(xim'*Av(feval(h_int,bl,br),xim,Bll,Brl));
    end
    fprintf('\n****************************************************\n')
    fprintf('-> Energy = %.12f at T/J = %.2f\n',U/Z,T(iT));
    fprintf('****************************************************\n\n')
    ET(iT,1:2) = [T(iT),U/Z]; save ET ET
end

%=========================================================================
function [Bl,Br,d,L] = bloques_ini(S0,S1,celdaU,Nrp)
% BLOQUES_INI - Construye los bloques iniciales 

global b0 b1

% Construir patron de dimensiones de los sitios
L      = Nrp*length(celdaU);
d0     = 2*S0+1;
d1     = 2*S1+1;
vecP   = repmat(celdaU,[1 Nrp]);
d      = zeros(1,L);
d(1:L) = vecP(1:L)*d1 + (1-vecP(1:L))*d0; % dim_sitios

Bl = cell(1,L-3);
Br = cell(1,L-3);

% Determinar bloques extremos
if d(1) == b0.dim
    Bl{1} = b0;
else
    Bl{1} = b1;
end
if d(L) == b0.dim
    Br{1} = b0;
else
    Br{1} = b1;
end
%=========================================================================
function b = sitio(l,d)
% SITIO - Determina el sitio etiquetado por l

global b0 b1

if d(l) == b0.dim
    b = b0; 
else
    b = b1; 
end
%========================================================================
function Av = Av(A,v,Bl,Br) 
% Alv - Multiplicacion matriz-vector A{l+1,l+2} * v

dA  = size(A,1);
dL  = Bl.dim;
dR  = Br.dim;
w0  = reshape(v, [dR dL*dA]).';
w1  = A*reshape(w0, [dA dL*dR]);
w2  = reshape(w1, [dL*dA dR]).';
Av  = reshape(w2, [dL*dA*dR 1]);
%========================================================================

