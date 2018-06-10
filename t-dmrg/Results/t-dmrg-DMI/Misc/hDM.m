
function H = hDM(BL,BR)

global Dx Dy Dz J

Sxl = 0.5*(BL.Spb+(BL.Spb)');
Sxr = 0.5*(BR.Spb+(BR.Spb)');
Syl = -0.5*i*(BL.Spb-(BL.Spb)');
Syr = -0.5*i*(BR.Spb-(BR.Spb)');
Szl = BL.Szb;
Szr = BR.Szb;

% Productos entre dos sitios
SySz = kron(Syl,Szr); SzSy = kron(Szl,Syr);
SzSx = kron(Szl,Sxr); SxSz = kron(Sxl,Szr);
SxSy = kron(Sxl,Syr); SySx = kron(Syl,Sxr);
SzSz = kron(Szl,Szr); SxSx = kron(Sxl,Sxr);
SySy = kron(Syl,Syr);
% Interaccion D-M 
%H = J*SzSz + Dx*(SySz-SzSy)+Dy*(SzSx-SxSz)+Dz*(SxSy-SySx);
%H = J*(SxSx+SySy+SzSz) + Dx*(SySz-SzSy)+Dy*(SzSx-SxSz)+Dz*(SxSy-SySx);
H = J*(SxSx+SySy) + Dx*(SySz-SzSy)+Dy*(SzSx-SxSz)+Dz*(SxSy-SySx);
    

        