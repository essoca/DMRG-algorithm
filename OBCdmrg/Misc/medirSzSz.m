
function [SzSz11,SzSz00,SzSz10] = medirSzSz(Phi,Bl,Br,d,m,h_int)

L  = length(d);
Op = @(l) operador(l,d);
b  = @(l) sitio(l,d);

Oi = kron(Op(1),speye(d(2)));
save('tmp/op1','Oi')
Oi = kron(speye(d(1)),Op(2));
save('tmp/op2','Oi')
OOij = kron(Op(1),Op(2));
save('tmp/oop12','OOij')

for l = 2:L-3
    Phi0  = chbasis(Phi(:,1),l-1,l,Bl,Br,d);
    [Bl{l+1},Phi] = sweepR(Bl{l},b(l+1),b(l+2),Br{L-l-2},m,h_int,Phi0);
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
    Phi0  = chbasis(Phi(:,1),l+1,l,Bl,Br,d);
    [Br{L-l-1},Phi] = sweepL(Bl{l},b(l+1),b(l+2),Br{L-l-2},m,h_int,Phi0);
    wcorreR(l,Br,d,Op)
end

% Medir correlaciones entre sitios dentro de un mismo bloque
profSzSz = zeros(L,L);
dimBL = Bl{L/2-1}.dim * d(L/2);
dimBR = d(L/2+1) * Br{L/2-1}.dim;
[rhoL,rhoR] = calcular_rho(Phi,dimBL,dimBR);

for i = 1:L/2 
    for j = i+1:L/2
        load(strcat('tmp/oop',int2str(i),int2str(j))); 
        profSzSz(i,j) = trace(rhoL*OOij); clear OOij;
        load(strcat('tmp/oop',int2str(L-j+1),int2str(L-i+1)));
        profSzSz(L-j+1,L-i+1) = trace(rhoR*OOij); clear OOij;
    end
end

% Medir correlaciones entre sitios dentro de bloques diferentes
for i = 1:L/2 
    for j = L/2+1:L
        load(strcat('tmp/op',int2str(i))); Oii = Oi;
        load(strcat('tmp/op',int2str(j))); Ojj = Oi; 
        profSzSz(i,j) = Phi(:,1)'* kron(Oii,Ojj) *Phi(:,1);
        clear OOij Oi;
    end
end

SzSz11 = zeros(L/2-1,1);
SzSz00 = zeros(L/2-1,1);
SzSz10 = zeros(L/2,1);

% Correlaciones entre sitios [1,1]
n = 0;
for k = 3:2:L-1
    n = n+1;
    SzSz11(n) = profSzSz(1,k) - profSzSz(1,1)*profSzSz(k,k);
end
% Correlaciones entre sitios [0,0]
n = 0;
for k = 4:2:L
    n = n+1;
    SzSz00(n) = profSzSz(2,k) - profSzSz(2,2)*profSzSz(k,k);
end
% Correlaciones entre sitios [1,0]
n = 0;
for k = 2:2:L
    n = n+1;
    SzSz10(n) = profSzSz(1,k) - profSzSz(1,1)*profSzSz(k,k);
end

%=========================================================================
function wcorreL(l,Bl,d,Op)

lmas1 = l+1;
for i = 1:lmas1
    for j = i:lmas1 
        if (j ~= lmas1) && (j ~= i) 
            load(strcat('tmp/oop',int2str(i),int2str(j)));
            OOij = kron(Bl{l}.O' *OOij* Bl{l}.O,speye(d(lmas1)));
            save(strcat('tmp/oop',int2str(i),int2str(j)),'OOij')
            clear OOij
        elseif (j ~= lmas1) && (j == i)
            load(strcat('tmp/op',int2str(j))); Oii = Oi;
            Oi = kron(Bl{l}.O' *Oi* Bl{l}.O,speye(d(lmas1)));
            save(strcat('tmp/op',int2str(j)),'Oi')
            clear Oi
        elseif (j == lmas1) && (j ~= i)
            OOij = kron(Bl{l}.O' *Oii* Bl{l}.O,Op(lmas1));
            save(strcat('tmp/oop',int2str(i),int2str(j)),'OOij')
            clear OOij Oii
        elseif (j == lmas1) && (j == i)
            Oi = kron(speye(Bl{l}.dim),Op(lmas1));
            save(strcat('tmp/op',int2str(j)),'Oi')
            clear Oi
        end
    end
end
%=========================================================================
function wcorreR(l,Br,d,Op)

L = size(Br,2)+2; lmas2 = l+2; Lmlm2 = L-l-2;
for i = L:-1:lmas2
    for j = i:-1:lmas2
        if (j ~= lmas2) && (j ~= i)
            load(strcat('tmp/oop',int2str(j),int2str(i)));
            OOij = kron(speye(d(lmas2)),Br{Lmlm2}.O' *OOij* Br{Lmlm2}.O);
            save(strcat('tmp/oop',int2str(j),int2str(i)),'OOij')
            clear OOij
        elseif (j ~= lmas2) && (j == i)
            load(strcat('tmp/op',int2str(j))); Oii = Oi;
            Oi = kron(speye(d(lmas2)),Br{Lmlm2}.O' *Oi* Br{Lmlm2}.O);
            save(strcat('tmp/op',int2str(j)),'Oi')
            clear Oi
        elseif (j == lmas2) && (j ~= i)
            OOij = kron(Op(lmas2),Br{Lmlm2}.O' *Oii* Br{Lmlm2}.O);
            save(strcat('tmp/oop',int2str(j),int2str(i)),'OOij')
            clear OOij Oii
        elseif (j == lmas2) && (j == i)
            Oi = kron(Op(lmas2),speye(Br{Lmlm2}.dim));
            save(strcat('tmp/op',int2str(j)),'Oi')
            clear Oi
        end
    end
end
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
function b = sitio(l,d)
% SITIO - Determina el sitio etiquetado por l

global b0 b1

if d(l) == b0.dim
    b = b0;
else
    b = b1;
end
%=========================================================================
