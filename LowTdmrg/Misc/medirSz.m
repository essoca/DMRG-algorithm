
function Szprof = medirSz(xim,Bl,Br,d)

L  = length(d);
Op = @(l) operador(l,d); 

Szprof    = zeros(L,1);
Szprof(1) = real(xim'*A1v(kron(Op(1),eye(d(2))),xim,Br,d));
Szprof(2) = real(xim'*A1v(kron(eye(d(1)),Op(2)),xim,Br,d));
for l = 2:L-3
    xim         = chbasis(xim,l-1,l,Bl,Br,d);
    Szprof(l+1) = real(xim'*Alv(kron(Op(l+1),eye(d(l+2))),xim,l,Bl,Br,d));
end
Szprof(L-1) = real(xim'*ALv(kron(Op(L-1),eye(d(L))),xim,Bl,d));
Szprof(L)   = real(xim'*ALv(kron(eye(d(L-1)),Op(L)),xim,Bl,d));

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