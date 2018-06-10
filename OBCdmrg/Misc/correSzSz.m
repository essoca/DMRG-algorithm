
function SzSzr = correSzSz(Phi,Bl,Br,d)
% Mide correlaciones SzSz con respecto al centro

if ~exist('tmp','dir')
    mkdir('tmp')
end
% Identificar operadores
Op = @(l) operador(l,d);
L = length(d);

Oi = kron(Op(1),speye(d(2)));
save('tmp/op1','Oi')
Oi = kron(speye(d(1)),Op(2));
save('tmp/op2','Oi')
clear Oi

for l = 2:L-3
    Phi = chbasis(Phi,l-1,l,Bl,Br,d);
    if l <= L/2-1
        lmas1 = l+1;
        for n = 1:lmas1
            if n <= l
                load(strcat('tmp/op',int2str(n)));
                Oi = kron(Bl{l}.O' *Oi* Bl{l}.O,speye(d(lmas1)));
                save(strcat('tmp/op',int2str(n)),'Oi')
                clear Oi
            else
                Oi = kron(speye(Bl{l}.dim),Op(lmas1));
                save(strcat('tmp/op',int2str(n)),'Oi')
                clear Oi
            end
        end
    end
end

Oi = kron(speye(d(L-1)),Op(L));
save(strcat('tmp/op',int2str(L)),'Oi')
Oi = kron(Op(L-1),speye(d(L)));
save(strcat('tmp/op',int2str(L-1)),'Oi')

for l = L-4:-1:L/2-1
    Phi = chbasis(Phi,l+1,l,Bl,Br,d);
    lmas2 = l+2; Lmlm2 = L-l-2;
    for n = L:-1:lmas2
        if n >= l+3
            load(strcat('tmp/op',int2str(n)));
            Oi = kron(speye(d(lmas2)),Br{Lmlm2}.O' *Oi* Br{Lmlm2}.O);
            save(strcat('tmp/op',int2str(n)),'Oi')
            clear Oi
        else
            Oi = kron(Op(lmas2),speye(Br{Lmlm2}.dim));
            save(strcat('tmp/op',int2str(n)),'Oi')
            clear Oi
        end
    end
end

dL = Bl{L/2-1}.dim*d(L/2); dR = d(L/2+1)*Br{L/2-1}.dim;
clear Bl Br d

% Corre entre sitios equidistantes al centro de la cadena
SzSzr = zeros(L-1,1);
for r = 1:L-1
    if rem(r,2) == 0
        j =0.5*(L-r); k = 0.5*(L+r);
    else
        j =0.5*(L-r+1); k = 0.5*(L+r+1);
    end
    fprintf('r = %d\n',r)
    load(strcat('tmp/op',int2str(j))); Oii = Oi; clear Oi
    load(strcat('tmp/op',int2str(k))); Ojj = Oi; clear Oi
    OiOjPhi = OOv(Oii,Ojj,Phi,dL,dR);
    SzSzr(r) = abs(Phi'*OiOjPhi); save SzSzr SzSzr
end
delete('tmp/*.mat')

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
function w = OOv(Oii,Ojj,v,dL,dR)
% OOV - Realiza la multiplicacion matriz-vector

w0 = reshape(v,[dR,dL]).';
w1 = Oii* w0* (Ojj.');                          
w  = real(reshape(w1.',[dL*dR,1]));
%=========================================================================