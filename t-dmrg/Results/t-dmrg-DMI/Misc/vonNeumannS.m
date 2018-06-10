
function S = vonNeumannS(Phi,l,Bl,Br,d)
% VONNEUMANNS - Calcula la entropia de Von Neumann del bloque con l sitios

dBle = Bl{l}.dim*d(l+1);
dBre = d(l+2)*Br{length(d)-l-2}.dim;

if size(Phi,1) ~= dBle*dBre
    % Se asume que Phi esta en la configuracion Bl{1}b(2)b(3)Br{L-3}
    for k = 2:l
        Phi = chbasis(Phi,k-1,k,Bl,Br,d);
    end
    % Asi expresamos Phi en la configuracion  Bl{l}b(l+1)b(l+2)Br{L-l-2}
end

[VL,DL]  = eig(calcular_rho(Phi,dBle,dBre,'L')); 
rho      = diag(DL);

% Calcular la entropia de von Neumann
S = -real(sum(rho.*log2(rho)));
