
function [Blnew,Brnew,Emin,Emax,uPm] = dmrg_infinito(Bl,bl,br,Br,err,h_int,beta)
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

% Alargar bloques
Ble = alargar(Bl,bl,'haciaR',h_int);
Bre = alargar(Br,br,'haciaL',h_int);
% Encontrar Emin y Emax del superbloque
[Emin,Emax] = superbloque(Ble,Bre,h_int);
% Crear xi aleatorio
xi = 2*rand(Ble.dim*Bre.dim,1)-1; xi = xi/norm(xi);
% Obtener el target state xim
xim = getxim(xi,Emin,Emax,Ble,Bre,h_int,beta);
% Calcular matriz densidad con xim
[rhoL,rhoR] = calcular_rho(xim/norm(xim),Ble.dim,Bre.dim);
% Truncar bloques
[Blnew,uPm] = truncar(Ble,rhoL,err);
Brnew = truncar(Bre,rhoR,err);
