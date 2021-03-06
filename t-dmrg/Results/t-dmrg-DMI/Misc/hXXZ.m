
function H = hXXZ(BL,BR,varargin)

global Jxy Jz

switch nargin
    case 2
        SpSm = kron(BL.Spb,BR.Spb');
        SzSz = kron(BL.Szb,BR.Szb);
        H    = 0.5*Jxy*(SpSm + SpSm') + Jz*SzSz;
    case 3 
        % Multiplicacion matriz-vector
        w0 = varargin{1};
        H  = 0.5*Jxy* BL.Spb * w0 * (BR.Spb').' + ...
             0.5*Jxy* BL.Spb' * w0 * (BR.Spb).' + ...
             Jz* BL.Szb * w0 *(BR.Szb).';
    otherwise
        error('Numero invalido de argumentos')
end
        