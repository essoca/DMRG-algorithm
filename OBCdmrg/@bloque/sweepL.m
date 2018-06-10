
function [Br_new,Phi,E,mnew] = sweepL(Bl,bl,br,Br,err,h_int,Phi0)
% SWEEPL - Ejecuta un paso del algoritmo DMRG finito (derecha-izquieda)
%-------------------------------------------------------------------------
%    Calcula el bloque derecho que se va a usar en el siguiente paso del 
%    algoritmo, haciendo que las propiedades del estado base converjan al
%    estado base real de una cadena de espines mixtos.
%-------------------------------------------------------------------------
%  Sintaxis:
%    [Br_new,Phi,E,uPm] = sweepL(Bl,bl,br,Br,m,h_int,xdens,Phi0);
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
%    Br_new : Bloque Br alargado hacia la izquierda y truncado
%    Phi    : Estado blanco del superbloque de longitud L
%    E      : Energia correpondiente
%    uPm    : Estimativo del error de truncamiento (1-Pm)

% Paso sistema finito
Ble           = alargar(Bl,bl,'haciaR',h_int);
Bre           = alargar(Br,br,'haciaL',h_int);
[Phi,E]       = superbloque(Ble,Bre,h_int,Phi0); 
rhoR          = calcular_rho(Phi,Ble.dim,Bre.dim,'R');
[Br_new,mnew] = truncar(Bre,rhoR,err);
