
function H = hdimer(BL,BR,varargin)

global delta dd

ind = find(max(BL.dim,BR.dim)==[BL.dim,BR.dim]);
if ind(1) == 1, lB = BL.l; else lB = BR.l; end
% Determinar magnitud de interaccion
if rem(lB,2) ~= 0
    Jxy = 2*(1-delta); Jz = 2*delta;
else
    Jxy = 2*dd*(1-delta); Jz = 2*dd*delta;
end

switch nargin
    case 2
        % Interaccion bloques
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