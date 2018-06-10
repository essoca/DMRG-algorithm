
function [Bl_new,uPm] = sweepR(Bl,dl,dr,Br,m,Phi)

% Construir matriz densidad con el estado Phi
dBle = Bl.dim * dl;
dBre = dr * Br.dim;
rhoL = calcular_rho(Phi,dBle,dBre,'L');

% Construir la base reducida
m           = min(dBle,m);        
[vrho wrho] = eig(rhoL); 
[wrho irho] = sort(diag(wrho),'descend'); 
vrho        = vrho(:,irho);
mnew        = check_m(wrho,m);
O           = vrho(:,1:mnew); 

% Calculmos uPm = 1-Pm
uPm = 1-sum(wrho(1:mnew));

% Construir el nuevo bloque izquierdo
Bl_new = bloque_simple(Bl.l+1, mnew, O);
