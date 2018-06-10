
function chiT = correSzSz(Phi,d)
% Medir correlaciones para calcular chi*T = sum_i <Sz0*Szi>

if ~exist('tmp/Oper','dir')
    mkdir('tmp/Oper')
end
% Identificar operadores del sitio l
Op = @(l) operador(l,d); 
L = length(d);

Oi = kron(Op(1),speye(d(2)));
save('tmp/Oper/op1','Oi')
Oi = kron(speye(d(1)),Op(2));
save('tmp/Oper/op2','Oi')
clear Oi

for l = 2:L-3
    Phi = chbasis(Phi,l-1,l,d);
    if l <= L/2-1
        lmas1 = l+1;
        for n = 1:lmas1
            if n <= l
                load(strcat('tmp/Oper/op',int2str(n)));
                Bl = rbloque('l',l);
                Oi = kron(Bl.O' *Oi* Bl.O,speye(d(lmas1)));
                save(strcat('tmp/Oper/op',int2str(n)),'Oi')
                clear Oi Bl
            else
                Bl = rbloque('l',l);
                Oi = kron(speye(Bl.dim),Op(lmas1));
                save(strcat('tmp/Oper/op',int2str(n)),'Oi')
                clear Oi Bl
            end
        end
    end
end

Oi = kron(speye(d(L-1)),Op(L));
save(strcat('tmp/Oper/op',int2str(L)),'Oi')
Oi = kron(Op(L-1),speye(d(L)));
save(strcat('tmp/Oper/op',int2str(L-1)),'Oi')

for l = L-4:-1:L/2-1
    Phi = chbasis(Phi,l+1,l,d);
    lmas2 = l+2; Lmlm2 = L-l-2;
    for n = L:-1:lmas2
        if n >= l+3
            load(strcat('tmp/Oper/op',int2str(n)));
            Br = rbloque('r',Lmlm2);
            Oi = kron(speye(d(lmas2)),Br.O' *Oi* Br.O);
            save(strcat('tmp/Oper/op',int2str(n)),'Oi')
            clear Oi Br
        else
            Br = rbloque('r',Lmlm2);
            Oi = kron(Op(lmas2),speye(Br.dim));
            save(strcat('tmp/Oper/op',int2str(n)),'Oi')
            clear Oi Br
        end
    end
end
Bl = rbloque('l',L/2-1); Br = rbloque('r',L/2-1);
dL = Bl.dim*d(L/2); dR = d(L/2+1)*Br.dim;
[rhoL,rhoR] = calcular_rho(Phi,dL,dR);
clear Bl Br d

% Calcular la suma: chi*T = sum_i <Sz0*Szi>
% Corre entre sitios equidistantes al centro de la cadena
SzSzr = zeros(L-1,1);
for r = 1:L-1
    if rem(r,2) == 0
        j =0.5*(L-r); k = 0.5*(L+r);
    else
        j =0.5*(L-r+1); k = 0.5*(L+r+1);
    end
    load(strcat('tmp/Oper/op',int2str(j))); Oii = Oi; clear Oi
    load(strcat('tmp/Oper/op',int2str(k))); Ojj = Oi; clear Oi
    OiOj = trace(rhoL*Oii)*trace(rhoR*Ojj);
    OiOjPhi = OOv(Oii,Ojj,Phi,dL,dR);
    SzSzr(r) = Phi'*OiOjPhi-OiOj;
end
chiT = abs(sum(SzSzr));
delete('tmp/Oper/*.mat')
%=========================================================================
function Op = operador(l,d)
% Operador Sz del sitio l

global Sz0 Sz1

if d(l) == size(Sz0,1)
    Op = Sz0;
else
    Op = Sz1;
end
%=========================================================================
function B = rbloque(strb,l)
% Lee del disco el bloque B de l sitios. Aqui strb es: 'l' o 'r'

fname = strcat('tmp/B',strb,'(',num2str(l),').mat');
load(fname,'B')
%=========================================================================
function w = OOv(Oii,Ojj,v,dL,dR)
% OOV - Realiza la multiplicacion matriz-vector

w0 = reshape(v,[dR,dL]).';
w1 = Oii* w0* (Ojj.');                          
w  = real(reshape(w1.',[dL*dR,1]));
%=========================================================================
