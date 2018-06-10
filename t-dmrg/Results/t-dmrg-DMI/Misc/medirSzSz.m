
function [SzSz,SzSzr] = medirSzSz(Phi,Bl,Br,d)

global dirop

L  = length(d); Lm = L/2;
Op = @(l) operador(l,d);

if ~exist(dirop,'dir')
    mkdir(dirop)
end
Oi = kron(Op(1),speye(d(2)));
save(strcat(dirop,'/op1'),'Oi')
Oi = kron(speye(d(1)),Op(2));
save(strcat(dirop,'/op2'),'Oi')
OOij = kron(Op(1),Op(2));
save(strcat(dirop,'/oop12'),'OOij')

for l = 2:L-3
    Phi = chbasis(Phi(:,1),l-1,l,Bl,Br,d);
    if l <= Lm-1
        wcorreL(l,Bl,d,Op)
    end
end

Oi = kron(speye(d(L-1)),Op(L));
save(strcat(strcat(dirop,'/op'),int2str(L)),'Oi')
Oi = kron(Op(L-1),speye(d(L)));
save(strcat(strcat(dirop,'/op'),int2str(L-1)),'Oi')
OOij = kron(Op(L-1),Op(L));
save(strcat(strcat(dirop,'/oop'),int2str(L-1),int2str(L)),'OOij')

for l = L-4:-1:Lm-1
    Phi = chbasis(Phi(:,1),l+1,l,Bl,Br,d);
    wcorreR(l,Br,d,Op)
end

dimBL = Bl{Lm-1}.dim * d(Lm);
dimBR = d(Lm+1) * Br{Lm-1}.dim;
[rhoL,rhoR] = calcular_rho(Phi,dimBL,dimBR);
SzSz = zeros(L-1,1);
SzSzr = zeros(Lm,1);

for i = 1:Lm-1 
    % sitios vecinos en bloque izquierdo
    j = i+1;
    load(strcat(strcat(dirop,'/oop'),int2str(i),int2str(j)));
    load(strcat(strcat(dirop,'/op'),int2str(i))); Oii = Oi;
    load(strcat(strcat(dirop,'/op'),int2str(j))); Ojj = Oi;
    SzSz(i) = real(trace(rhoL*OOij)-trace(rhoL*Oii)*trace(rhoL*Ojj)); 
    clear OOij Oii Ojj Oi;
    % sitios vecinos en bloque derecho
    j = L-i;
    load(strcat(strcat(dirop,'/oop'),int2str(j),int2str(L-i+1)));
    load(strcat(strcat(dirop,'/op'),int2str(j))); Ojj = Oi;
    load(strcat(strcat(dirop,'/op'),int2str(L-i+1))); Oii = Oi;
    SzSz(L-i) = real(trace(rhoR*OOij)-trace(rhoR*Ojj)*trace(rhoR*Oii)); 
    clear OOij Oii Ojj Oi;    
    % par de sitios simetricos al centro 
    load(strcat(strcat(dirop,'/op'),int2str(Lm-i+1))); Oii = Oi;
    load(strcat(strcat(dirop,'/op'),int2str(Lm+i))); Ojj = Oi;
    OiOjPhi = OOv(Oii,Ojj,Phi,dimBL,dimBR);
    SzSzr(i) = real(Phi'*OiOjPhi-trace(rhoL*Oii)*trace(rhoR*Ojj));
end
% Par de sitios extremos 
load(strcat(strcat(dirop,'/op'),int2str(1))); Oii = Oi;
load(strcat(strcat(dirop,'/op'),int2str(L))); Ojj = Oi;
OiOjPhi = OOv(Oii,Ojj,Phi,dimBL,dimBR);
SzSzr(Lm) = real(Phi'*OiOjPhi-trace(rhoL*Oii)*trace(rhoR*Ojj));
SzSz(Lm) = SzSzr(1);

delete(strcat(dirop,'/*.mat'))
%=========================================================================
function wcorreL(l,Bl,d,Op)

global dirop
lmas1 = l+1;
for i = 1:l
    % Actualizar operadores i dentro del bloque l
    load(strcat(strcat(dirop,'/op'),int2str(i))); Oii = Oi;
    Oi = kron(Bl{l}.O' *Oi* Bl{l}.O,speye(d(lmas1)));
    save(strcat(strcat(dirop,'/op'),int2str(i)),'Oi')
    clear Oi
    % Dentro del mismo bloque, solo corre entre vecinos
    j = i+1;
    if j ~=lmas1
        % Operadores i con j dentro del bloque l
        load(strcat(strcat(dirop,'/oop'),int2str(i),int2str(j)));
        OOij = kron(Bl{l}.O' *OOij* Bl{l}.O,speye(d(lmas1)));
        save(strcat(strcat(dirop,'/oop'),int2str(i),int2str(j)),'OOij')
        clear OOij
    else
        % Operadores i dentro del bloque l con j=l+1
        OOij = kron(Bl{l}.O' *Oii* Bl{l}.O,Op(lmas1));
        save(strcat(strcat(dirop,'/oop'),int2str(i),int2str(j)),'OOij')
        clear OOij Oii
    end
end
% Actualizar operador l+1
Oi = kron(speye(Bl{l}.dim),Op(lmas1));
save(strcat(strcat(dirop,'/op'),int2str(lmas1)),'Oi')
%=========================================================================
function wcorreR(l,Br,d,Op)

global dirop
L = size(Br,2)+2; lmas2 = l+2; Lmlm2 = L-l-2;
for i = L:-1:l+3
    % Actualizar operadores i dentro del bloque L-l-2
    load(strcat(strcat(dirop,'/op'),int2str(i))); Oii = Oi;
    Oi = kron(speye(d(lmas2)),Br{Lmlm2}.O' *Oi* Br{Lmlm2}.O);
    save(strcat(strcat(dirop,'/op'),int2str(i)),'Oi')
    clear Oi
    % Dentro del mismo bloque, solo corre entre vecinos
    j = i-1;
    if j ~= lmas2
        % Operadores i con j dentro del bloque L-l-2
        load(strcat(strcat(dirop,'/oop'),int2str(j),int2str(i)));
        OOij = kron(speye(d(lmas2)),Br{Lmlm2}.O' *OOij* Br{Lmlm2}.O);
        save(strcat(strcat(dirop,'/oop'),int2str(j),int2str(i)),'OOij')
        clear OOij
    else
        % Operadores i dentro del bloque L-l-2 con j=l+2
        OOij = kron(Op(lmas2),Br{Lmlm2}.O' *Oii* Br{Lmlm2}.O);
        save(strcat(strcat(dirop,'/oop'),int2str(j),int2str(i)),'OOij')
        clear OOij Oii
    end
end
% Actualizar operador l+2
Oi = kron(Op(lmas2),speye(Br{Lmlm2}.dim));
save(strcat(strcat(dirop,'/op'),int2str(lmas2)),'Oi')
%=========================================================================
function w = OOv(Oii,Ojj,v,dL,dR)
% OOV - Realiza la multiplicacion matriz-vector

w0 = reshape(v,[dR,dL]).';
w1 = Oii* w0* (Ojj.');                          
w  = real(reshape(w1.',[dL*dR,1]));
%=========================================================================
function Op = operador(l,d) 
% OPERADOR - Determina el operador a medir en el sitio l

global b0 b1

if d(l) == b0.dim
    Op = b0.Szb;
else
    Op = b1.Szb;
end
%=========================================================================

