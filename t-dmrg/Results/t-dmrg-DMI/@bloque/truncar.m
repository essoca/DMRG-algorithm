
function [Be_new,uPm] = truncar(Be,rho,m0)
% TRUNCAR - Trunca un bloque DMRG
%---------------------------------------------------------------
%    Trunca un bloque alargado existente a la base reducida
%    La operacion se simboliza como B(l+1,mxd) --> B(l+1,m)
%---------------------------------------------------------------
%  Sintaxis:
%    [Be_new,uPm] = truncar(Be,rho);
%  Entrada:
%    Be     : Bloque alargado a truncar
%    rho    : Matriz densidad reducida del bloque alargado
%    m      : Numero de autoestados de rho a retener para formar la base
%             truncada de los bloques alargados
%  Salida:
%    Be_new : Bloque truncado
%    uPm    : Estimativo del error en el truncamiento (1-Pm)

m             = min(Be.dim,m0);              
[vrho wrho]   = eig(rho);                   
[wrho irho]   = sort(diag(wrho),'descend'); 
vrho          = vrho(:,irho);   
mnew          = check_m(wrho,m);
O             = vrho(:,1:mnew); 

% Truncamos los operadores del bloque alargado
HB  = O' * Be.HB  * O; 
Spb = O' * Be.Spb * O;
Szb = O' * Be.Szb * O;

% Calculmos uPm = 1-Pm
uPm = 1-sum(wrho(1:mnew));

% Construir bloque truncado
Be_new = bloque(Be.l, mnew, HB, Spb, Szb, O);
