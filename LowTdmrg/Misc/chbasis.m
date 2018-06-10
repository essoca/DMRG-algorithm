
function Phi = chbasis(Phi,lb,la,Bl,Br,d)
% CHBASIS - Cambia la representacion de Phi

% Phi: B(lb)b(lb+1)b(lb+2)B(L-lb-2) --> Phi: B(la)b(la+1)b(la+2)B(L-la-2)

L   = length(d);

if la == lb
    % Transf. identidad
    return
elseif la-lb == 1
    % Transf. de White de izquierda a derecha -->
    Phi = A1Bv( Bl{la}.O', d(la+1), Br{L-lb-2}.O, Phi); 
elseif lb-la == 1
    % Transf. de White de derecha a izquierda <--
    Phi = A1Bv( Bl{lb}.O, d(lb+1), Br{L-la-2}.O', Phi);
else
   error('Las l-configuraciones no son consecutivas')
end

%========================================================================
function w = A1Bv (A,dn,B,v)
% kronA1Bv - Calcula w = kron(A,kron(eye(dn),B))*v
  
[dAi,dAj] = size(A);
[dBi,dBj] = size(B);

W0 = B * reshape(v, [dBj,dAj*dn]);
w0 = reshape(W0, [dAj*dn*dBi,1]);
W1 = reshape(w0, [dn*dBi,dAj]) * A.';
w1  = reshape(W1, [dAi*dn*dBi,1]);
w  = w1/norm(w1);

%========================================================================