
function [Phi,E] = superbloque(varargin)
% SUPERBLOQUE - Calcula autoestados del superbloque
%--------------------------------------------------------------------------
%    Construimos el superbloque del sistema al conectar el bloque izquierdo  
%    alargado Ble con el bloque derecho alargado Bre. Esta conexion se hace
%    en el Hamiltoniano, el cual se diagonaliza finalmente para obtener
%    los autoestados que se buscan con las autoenergias correspondientes.
%--------------------------------------------------------------------------
%  Sintaxis:
%    [Phi,E] = superbloque(Ble,Bre,h_int);      
%  Entrada:
%     Ble   : Bloque izquierdo alargado
%     Bre   : Bloque derecho alargado
%     h_int : Hamiltoniano de interaccion entre sitios vecinos (@)
%  Salida:
%    Phi    : Estados que se buscan
%     E     : Energias correspondientes

[Phi,E] = diagIrbl(varargin{:});
%[Phi,E] = diagEigs(varargin{:});


