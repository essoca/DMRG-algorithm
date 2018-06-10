
function H = hDM(BL,BR,varargin)

global Dx Dy Dz J

Sxl = 0.5*(BL.Spb+(BL.Spb)');
Sxr = 0.5*(BR.Spb+(BR.Spb)');
Syl = -0.5*i*(BL.Spb-(BL.Spb)');
Syr = -0.5*i*(BR.Spb-(BR.Spb)');
Szl = BL.Szb;
Szr = BR.Szb;

switch nargin
    case 2
        % Productos entre dos sitios
        SySz = kron(Syl,Szr); SzSy = kron(Szl,Syr);
        SzSx = kron(Szl,Sxr); SxSz = kron(Sxl,Szr);
        SxSy = kron(Sxl,Syr); SySx = kron(Syl,Sxr);
        % Interaccion D-M
        H = Dx*(SySz-SzSy)+Dy*(SzSx-SxSz)+Dz*(SxSy-SySx) + ...
            J*(kron(Sxl,Sxr)+kron(Syl,Syr));
    case 3 
        % Multiplicacion matriz-vector
        w0 = varargin{1};
        H  = Dx* Syl * w0 * Szr.' - Dx* Szl * w0 * Syr.' + ...
            Dy* Szl * w0 * Sxr.' - Dy* Sxl * w0 * Szr.' + ...
            Dz* Sxl * w0 * Syr.' - Dz* Syl * w0 * Sxr.' + ...
            J* Sxl * w0 * Sxr.' + J* Syl * w0 * Syr.';
    otherwise
        error('Numero invalido de argumentos')
end
        