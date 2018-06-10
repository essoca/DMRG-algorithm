
function [Szprof,E0,uPm] = gsdmrg(S0,S1,celdaU,Nrp,err,Nsweeps,h_int)
%function [E,uPm] = gsdmrg(S0,S1,celdaU,Nrp,err,Nsweeps,h_int)

global dispR 

% Bloques iniciales 
[Bl,Br,d,L] = bloques_ini(S0,S1,celdaU,Nrp);

% Funcion para identificar los sitios de la red
b = @(l) sitio(l,d);

[Phi2,E2] = eig(full(feval(h_int,b(1),b(L))));
E_ = E2(1);

% Sweep inicial
if dispR, fprintf('-> Algoritmo Infinito\n'), end
for l = 1:L/2-1
    [Bl{l+1},Br{l+1},Phi,E,uPm] = ...
        dmrg_infinito(Bl{l},b(l+1),b(L-l),Br{l},err(1),h_int);   
    E0 = 0.5*(E(1)-E_); E_ = E(1);
    fprintf('LBL = %d, LBR = %d, E0 = %.8f, m = %d\n',l+1,l+1,E0,uPm)
end

if dispR, fprintf('-> Sweep inicial\n'), end
for l = L/2:L-3
    Phi0  = chbasis(Phi(:,1),l-1,l,Bl,Br,d);
    [Bl{l+1},Phi,E,uPm] = ...
        sweepR(Bl{l},b(l+1),b(l+2),Br{L-l-2},err(1),h_int,Phi0);
    %fprintf('LBL = %d, LBR = %d, E = %.8f, uPm = %.0E\n',l+1,L-l-1,E,uPm)
end
for l = L-4:-1:1
    Phi0  = chbasis(Phi(:,1),l+1,l,Bl,Br,d);
    [Br{L-l-1},Phi,E,uPm] = ...
        sweepL(Bl{l},b(l+1),b(l+2),Br{L-l-2},err(1),h_int,Phi0);
end

% Sweeps de convergencia
for n = 1: Nsweeps-1
    if dispR, fprintf('-> Sweep %d \n',n+1), end
    for l = 2:L-3
        Phi0  = chbasis(Phi(:,1),l-1,l,Bl,Br,d);
        [Bl{l+1},Phi,E,uPm] = ...
            sweepR(Bl{l},b(l+1),b(l+2),Br{L-l-2},err(n+1),h_int,Phi0);
        fprintf('LBL = %d, LBR = %d, E = %.8f, uPm = %.0E\n',l+1,L-l-1,E,uPm)
    end
    for l = L-4:-1:1
        Phi0  = chbasis(Phi(:,1),l+1,l,Bl,Br,d);
        [Br{L-l-1},Phi,E,uPm] = ...
            sweepL(Bl{l},b(l+1),b(l+2),Br{L-l-2},err(n+1),h_int,Phi0);
        fprintf('LBL = %d, LBR = %d, E = %.8f, uPm = %.0E\n',l+1,L-l-1,E,uPm)
    end
end
%Szprof = medirSz(Phi(:,1),Bl,Br,d);
%E = vesperadoE(Phi,h_int,d);
%=========================================================================
function [Bl,Br,d,L] = bloques_ini(S0,S1,celdaU,Nrp)
% BLOQUES_INI - Construye los bloques iniciales 

global b0 b1

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
%=========================================================================
