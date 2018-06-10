
function [Bl_new,Phi,E,uPm] = sweepR(Bl,bl,br,Br,m,h_int,Phi0)
% SWEEPR - Ejecuta un paso del algoritmo DMRG finito (izquieda-derecha)
%-------------------------------------------------------------------------
%    Calcula el bloque izquierdo que se va a usar en el siguiente paso del 
%    algoritmo, haciendo que las propiedades del estado base converjan al
%    estado base real de una cadena de espines mixtos.
%-------------------------------------------------------------------------
%  Sintaxis:
%    [Bl_new,Phi,E,uPm] = sweepR(Bl,bl,br,Br,m,h_int,xdens,Phi0);
%  Entrada:
%    Bl     : Bloque izquierdo de longitud l
%    bl     : Sitio agregado a Bl
%    Br     : Bloque derecho de longitud L-l-2
%    br     : Sitio agregado a Br
%    m      : Numero de autoestados de rho a retener para formar la base
%             truncada de los bloques alargados
%    h_int  : Hamiltoniano de interaccion entre sitios vecinos (@)
%    Phi0   : Estado inicial para la rutina de diagonalizacion (ARPACK)
%  Salida:
%    Bl_new : Bloque Bl alargado hacia la derecha y truncado
%    Phi    : Estado blanco del superbloque de longitud L
%    E      : Energia correpondiente
%    uPm    : Estimativo del error de truncamiento (1-Pm)

% Paso sistema finito
Ble          = alargar(Bl,bl,'haciaR',h_int);
Bre          = alargar(Br,br,'haciaL',h_int);
[Phi,E]      = superbloque(Ble,Bre,h_int,Phi0); 
rhoL         = calcular_rho(Phi,Ble.dim,Bre.dim,'L');
[Bl_new,uPm] = truncar(Ble,rhoL,m);