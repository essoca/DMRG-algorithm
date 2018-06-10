
function B = bloque(varargin)
% BLOQUE - Constructor de un bloque DMRG
%---------------------------------------------------------------------
%     Bloque B(l,m) : Estructura fundamental de un sistema DMRG
%     Se caracteriza por tener un numero l de sitios internos descritos
%     por un espacio de Hilbert (en general truncado) de dimension m
%---------------------------------------------------------------------
%  Sintaxis:
%     B = bloque(l,dim,HB,Spb,Szb,O); 
%  Entrada: 
%    l   : Numero de sitios en el bloque
%    dim : Dimension del espacio de Hilbert del bloque = min(m,d^l)
%    HB  : Hamiltoniano del bloque
%    Spb : S+ del spin del extremo derecho del bloque 
%    Szb : Sz del spin del extremo derecho del bloque
%    O   : Matriz que transforma la base DMRG de B(l-1,m) a la de B(l,m)
%  Salida:
%    B   : Estructura con los campos de entrada

switch nargin    
    case 0  % Constructor vacio
        B.l   = 0;
        B.dim = 0;
        B.HB  = [];
        B.Spb = [];
        B.Szb = [];
        B.O   = [];
        B     = class(B,'bloque');  
    case 1  % Constructor de copia
        if(isa(varargin{1},'bloque'))
            B = varargin{1}; 
        else
            error('El argumento no es un bloque')
        end
    case 6 % Constructor con argumentos
        B.l   = varargin{1};
        B.dim = varargin{2};
        B.HB  = varargin{3};
        B.Spb = varargin{4};   
        B.Szb = varargin{5};
        B.O   = varargin{6};
        B     = class(B,'bloque');
    otherwise
        error('Numero errado de argumentos')
end
        


