
function E = vesperadoE(Phi,Bl,Br,h_int,d)

global b0
L = length(d);
% Hamiltoniano de interaccion de sitios vecinos x,y
h = @(x,y) feval(h_int,b0,b0);

% Calcular <Phi| H |Phi> descomponiendo H en operadores de 2 sitios
% Esto es: <H> = Sum_n <Phi| h_int(n,n+1)+Hloc(n)+Hloc(n+1) |Phi>
E = real(Phi'*A1v(h(1,2),Bl,Br,Phi,d));
for l = 1:L-3
    E = E + real(Phi'*Alv(h(l+1,l+2),Bl,Br,Phi,l,d));
    if l ~= L-3, Phi = chbasis(Phi,l,l+1,Bl,Br,d); end
end
E = E + real(Phi'*ALv(h(L-1,L),Bl,Br,Phi,d));

%=========================================================================
function Av = A1v(A,Bl,Br,v,d)
% A1v - Multiplicacion matriz-vector A{1,2} * v

L  = length(d);
dA = d(1)*d(2);
dR = d(3)*Br{L-3}.dim;

w0 = reshape(v,[dR dA])*A.';
Av = reshape(w0,[dA*dR 1]);
%========================================================================
function Av = Alv(A,Bl,Br,v,l,d) 
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
function Av = ALv(A,Bl,Br,v,d)
% ALv - Multiplicacion matriz-vector A{L-1,L} * v

L  = length(d);
dA = d(L-1)*d(L);
dL = Bl{L-3}.dim*d(L-2);

w0 = A*reshape(v,[dA dL]);
Av = reshape(w0,[dL*dA 1]);
%=========================================================================