
function H = h_ini(BL,BR,varargin)
% Interaccion nula

switch nargin
    case 2
        H  = zeros(BL.dim*BR.dim);
    case 3 
        % Multiplicacion matriz-vector
        w0 = varargin{1};
        H  = zeros(BL.dim) * w0;
    otherwise
        error('Numero invalido de argumentos')
end
        