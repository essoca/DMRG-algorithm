
function B = bloque_simple(varargin)
% BLOQUE_SIMPLE - Constructor de un bloque DMRG simple
%---------------------------------------------------------------------
%     Bloque B(l,m) : Estructura fundamental de un sistema DMRG
%     Se caracteriza por tener un numero l de sitios internos descritos
%     por un espacio de Hilbert (en general truncado) de dimension m
%---------------------------------------------------------------------
%  Sintaxis:
%     B = bloque(l,dim,HB,SZ,Spb,Szb,O); 
%  Entrada: 
%    l   : Numero de sitios en el bloque
%    dim : Dimension del espacio de Hilbert del bloque = min(m,d^l)
%    O   : Matriz de que transforma B(l,m) a B(l+1,m)
%  Salida:
%    B   : Estructura con los campos de entrada

switch nargin    
    case 0  % Constructor vacio
        B.l   = 0;
        B.dim = 0;
        B.O   = [];
        B     = class(B,'bloque_simple');  
    case 1  % Constructor de copia
        if(isa(varargin{1},'bloque_simple'))
            B = varargin{1}; 
        else
            error('El argumento no es un bloque_simple')
        end
    case 3 % Constructor con argumentos
        B.l   = varargin{1};
        B.dim = varargin{2};
        B.O   = varargin{3};
        B     = class(B,'bloque_simple'); 
    otherwise
        error('Numero errado de argumentos')
end