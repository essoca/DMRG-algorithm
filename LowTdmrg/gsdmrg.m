
function [Phi0,uPm,d] = gsdmrg(S0,S1,celdaU,Nrp,m)

global b0

% Bloques iniciales 
[d,L] = bloques_ini(S0,S1,celdaU,Nrp);

% Funcion para identificar los sitios de la red
b = @(l) sitio(l,d);

% Construir estado inicial
Bl = rbloque('l',1); Br = rbloque('r',1);
for l = 1:L/2-1
    [Bl,Br,uPm] = dmrg_infinito(Bl,b(l+1),b(L-l),Br,m(1));
    wbloque('l',l+1,Bl); wbloque('r',l+1,Br);
end

% Estado a beta = 0
Bl = rbloque('l',L/2-1); Br = rbloque('r',L/2-1);
Phi0 = kron(kron(Bl.v0,b0.v0),kron(b0.v0,Br.v0));  

% Adaptar la base a Phi0
for l = L/2:L-3
    Phi0 = chbasis(Phi0,l-1,l,d);
    Br = rbloque('r',L-l-2);
    [Bl,uPm] = sweepR(Bl,d(l+1),d(l+2),Br,m(1),Phi0);
    wbloque('l',l+1,Bl);
end
Br = rbloque('r',2);
for l = L-4:-1:1
    Phi0 = chbasis(Phi0,l+1,l,d);
    Bl = rbloque('l',l); 
    [Br,uPm] = sweepL(Bl,d(l+1),d(l+2),Br,m(1),Phi0);
    wbloque('r',L-l-1,Br);
end
%=========================================================================
function [d,L] = bloques_ini(S0,S1,celdaU,Nrp)
% BLOQUES_INI - Construye los bloques iniciales 

global b0 b1

L      = Nrp*length(celdaU);
d0     = 2*S0+1;
d1     = 2*S1+1;
vecP   = repmat(celdaU,[1 Nrp]);
d      = zeros(1,L);
d(1:L) = vecP(1:L)*d1 + (1-vecP(1:L))*d0; d = d.^2;

% Determinar bloques extremos
if d(1) == b0.dim
    Bl = b0;
else
    Bl = b1;
end
if d(L) == b0.dim
    Br = b0;
else
    Br = b1;
end
wbloque('l',1,Bl)
wbloque('r',1,Br)
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
function wbloque(strb,l,B)
% Escribe al disco el bloque B de l sitios. Aqui strb es: 'l' o 'r'

fname = strcat('tmp/B',strb,'(',num2str(l),').mat');
save(fname,'B')
%=========================================================================
function B = rbloque(strb,l)
% Lee del disco el bloque B de l sitios. Aqui strb es: 'l' o 'r'

fname = strcat('tmp/B',strb,'(',num2str(l),').mat');
load(fname,'B')
%=========================================================================
