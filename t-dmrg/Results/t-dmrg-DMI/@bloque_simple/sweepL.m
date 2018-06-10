
function [Br_new,uPm] = sweepL(Bl,dl,dr,Br,m,Phi)

% Construir matriz densidad con el estado Phi
dBle = Bl.dim * dl;
dBre = dr * Br.dim;
rhoR = calcular_rho(Phi,dBle,dBre,'R');

% Construir la base reducida
m           = min(dBre,m);        
[vrho wrho] = eig(rhoR);
[wrho irho] = sort(diag(wrho),'descend'); 
vrho        = vrho(:,irho);      
mnew        = check_m(wrho,m);
O           = vrho(:,1:mnew); 

% Calculmos uPm = 1-Pm
uPm = 1-sum(wrho(1:mnew));

% Construir el nuevo bloque izquierdo
Br_new = bloque_simple(Br.l+1, mnew, O);
