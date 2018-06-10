
function [SzSz,SzSzr] = correSzSz(Phi,Bl,Br,d)

L  = length(d);
Op = @(l) operador(l,d);

if ~exist('tmp','dir')
    mkdir('tmp')
else
    delete('tmp/*.mat')
end
Oi = kron(Op(1),speye(d(2)));
save('tmp/op1','Oi')
Oi = kron(speye(d(1)),Op(2));
save('tmp/op2','Oi')
OOij = kron(Op(1),Op(2));
save('tmp/oop12','OOij')

for l = 2:L-3
    Phi = chbasis(Phi(:,1),l-1,l,Bl,Br,d);
    if l <= L/2-1
        wcorreL(l,Bl,d,Op)
    end
end

Oi = kron(speye(d(L-1)),Op(L));
save(strcat('tmp/op',int2str(L)),'Oi')
Oi = kron(Op(L-1),speye(d(L)));
save(strcat('tmp/op',int2str(L-1)),'Oi')
OOij = kron(Op(L-1),Op(L));
save(strcat('tmp/oop',int2str(L-1),int2str(L)),'OOij')

for l = L-4:-1:L/2-1
    Phi = chbasis(Phi(:,1),l+1,l,Bl,Br,d);
    wcorreR(l,Br,d,Op)
end

dimBL = Bl{Lm-1}.dim * d(Lm);
dimBR = d(Lm+1) * Br{Lm-1}.dim;
[rhoL,rhoR] = calcular_rho(Phi,dimBL,dimBR);
SzSz = SzSzVeci(Phi,rhoL,rhoR,L);

% Connected spin-spin correlations
Gzz = zeros(L-1,1);
for r = 1:L-1
    termGzz = 0; 
    for k = 1:L-r
        load(strcat('tmp/oop',int2str(k),int2str(k+r)));
        load(strcat('tmp/op',int2str(k))); Okk = Oi;
        load(strcat('tmp/op',int2str(k+r))); Okr = Oi;
        % revisar corre entre sitios en B diferentes o dentro del mismo B
        SzkSzkr = real(trace(rhoL*OOij)-trace(rhoL*Okk)*trace(rhoL*Okr)); 
        termGzz = termGzz + SzkSzkr;
    end
    Gzz(r) = (1/(L-r))*termGzz;
end

delete('tmp/*.mat')
%=========================================================================
function SzSz = SzSzVeci(Phi,rhoL,rhoR,L)
% Correlaciones entre sitios vecinos a lo largo de la red

SzSz = zeros(L-1,1);
dimBL = size(rhoL,1);
dimBR = size(rhoR,1);
for i = 1:Lm-1 
    % sitios en bloque izquierdo
    j = i+1;
    load(strcat('tmp/oop',int2str(i),int2str(j)));
    load(strcat('tmp/op',int2str(i))); Oii = Oi;
    load(strcat('tmp/op',int2str(j))); Ojj = Oi;
    SzSz(i) = real(trace(rhoL*OOij)-trace(rhoL*Oii)*trace(rhoL*Ojj)); 
    clear OOij Oii Ojj Oi;
    % sitios en bloque derecho
    j = L-i;
    load(strcat('tmp/oop',int2str(j),int2str(L-i+1)));
    load(strcat('tmp/op',int2str(j))); Ojj = Oi;
    load(strcat('tmp/op',int2str(L-i+1))); Oii = Oi;
    SzSz(L-i) = real(trace(rhoR*OOij)-trace(rhoR*Ojj)*trace(rhoR*Oii)); 
    clear OOij Oii Ojj Oi;    
end
% Par de sitios del centro (bloques diferentes)
load(strcat('tmp/op',int2str(L/2))); Oii = Oi;
load(strcat('tmp/op',int2str(L/2+1))); Ojj = Oi;
OiOjPhi = OOv(Oii,Ojj,Phi,dimBL,dimBR);
SzSz(L/2) = real(Phi'*OiOjPhi-trace(rhoL*Oii)*trace(rhoR*Ojj));
%=========================================================================
function wcorreL(l,Bl,d,Op)

lmas1 = l+1;
for i = 1:l
    % Actualizar operadores i dentro del bloque l
    load(strcat('tmp/op',int2str(i))); Oii = Oi;
    Oi = kron(Bl{l}.O' *Oi* Bl{l}.O,speye(d(lmas1)));
    save(strcat('tmp/op',int2str(i)),'Oi')
    clear Oi
    for j = i+1:lmas1
        if j ~=lmas1
            % Operadores i con j dentro del bloque l
            load(strcat('tmp/oop',int2str(i),int2str(j)));
            OOij = kron(Bl{l}.O' *OOij* Bl{l}.O,speye(d(lmas1)));
            save(strcat('tmp/oop',int2str(i),int2str(j)),'OOij')
            clear OOij
        else
            % Operadores i dentro del bloque l con j=l+1
            OOij = kron(Bl{l}.O' *Oii* Bl{l}.O,Op(lmas1));
            save(strcat('tmp/oop',int2str(i),int2str(j)),'OOij')
            clear OOij Oii
        end
    end
end
% Actualizar operador l+1
Oi = kron(speye(Bl{l}.dim),Op(lmas1));
save(strcat('tmp/op',int2str(lmas1)),'Oi')
%=========================================================================
function wcorreR(l,Br,d,Op)

L = size(Br,2)+2; lmas2 = l+2; Lmlm2 = L-l-2;
for i = L:-1:l+3
    % Actualizar operadores i dentro del bloque L-l-2
    load(strcat('tmp/op',int2str(i))); Oii = Oi;
    Oi = kron(speye(d(lmas2)),Br{Lmlm2}.O' *Oi* Br{Lmlm2}.O);
    save(strcat('tmp/op',int2str(i)),'Oi')
    clear Oi
    for j = i-1:-1:lmas2
        if j ~= lmas2
            % Operadores i con j dentro del bloque L-l-2
            load(strcat('tmp/oop',int2str(j),int2str(i)));
            OOij = kron(speye(d(lmas2)),Br{Lmlm2}.O' *OOij* Br{Lmlm2}.O);
            save(strcat('tmp/oop',int2str(j),int2str(i)),'OOij')
            clear OOij           
        else
            % Operadores i dentro del bloque L-l-2 con j=l+2
            OOij = kron(Op(lmas2),Br{Lmlm2}.O' *Oii* Br{Lmlm2}.O);
            save(strcat('tmp/oop',int2str(j),int2str(i)),'OOij')
            clear OOij Oii            
        end
    end
end
% Actualizar operador l+2
Oi = kron(Op(lmas2),speye(Br{Lmlm2}.dim));
save(strcat('tmp/op',int2str(lmas2)),'Oi')
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
