
function H = hdimer(BL,BR,varargin)

global dd

ind = find(max(BL.dim,BR.dim)==[BL.dim,BR.dim]);
if ind(1) == 1, lB = BL.l; else lB = BR.l; end
% Determinar magnitud de interaccion
if rem(lB,2) ~= 0
    Jxy = 1; Jz = 1;
else
    Jxy = dd; Jz = dd;
end
        
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
        