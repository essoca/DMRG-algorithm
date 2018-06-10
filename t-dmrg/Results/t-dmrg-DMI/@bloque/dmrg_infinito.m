
function [Bl_new,Br_new,Phi,E,uPm] = dmrg_infinito(Bl,bl,br,Br,m,h_int)
% DMRG_INFINITO - Ejecuta un paso del algoritmo DMRG infinito
%-------------------------------------------------------------------------
%    Agrega dos sitios en el medio (OBC) de una cadena de espines mixtos
%    calcula la energia de esa configuracion y prepara el mismo calculo 
%    para el siguiente paso.
%-------------------------------------------------------------------------
%  Sintaxis:
%    [Bl_new,Br_new,E,uPm] = dmrg_infinito(Bl,dl,dr,Br,m);
%  Entrada:
%    Bl     : Bloque izquierdo de longitud l
%    bl     : Sitio agregado a Bl
%    Br     : Bloque derecho de longitud l
%    br     : Sitio agregado a Br
%    m      : Numero de autoestados de rho a retener para formar la base
%             truncada de los bloques alargados
%    h_int  : Hamiltoniano de interaccion entre sitios vecinos (@)
%  Salida:
%    Bl_new : Bloque Bl alargado hacia la derecha y truncado
%    Br_new : Bloque Br alargado hacia la izquierda y truncado
%    E      : Energia del superbloque de longitud 2*l+2
%    uPm    : Estimativo del error de truncamiento (1-Pm)

% Paso sistema infinito
Ble          = alargar(Bl,bl,'haciaR',h_int);
Bre          = alargar(Br,br,'haciaL',h_int);
[Phi,E]      = superbloque(Ble,Bre,h_int);
[rhoL,rhoR]  = calcular_rho(Phi,Ble.dim,Bre.dim);
[Bl_new,uPm] = truncar(Ble,rhoL,m);
Br_new       = truncar(Bre,rhoR,m);
