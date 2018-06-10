
function varargout = calcular_rho(varargin)

Phi = varargin{1};
dL  = varargin{2};
dR  = varargin{3};

if nargin == 3
    varargout{1} = matrizrho(Phi,dL,dR,1);
    varargout{2} = matrizrho(Phi,dL,dR,2);
else
    switch varargin{4}
        case 'L'
            varargout{1} = matrizrho(Phi,dL,dR,1);
        case 'R'
            varargout{1} = matrizrho(Phi,dL,dR,2);
        otherwise
            error('Opciones para Arg4 son los caracteres L o R')
    end
end

%=========================================================================
function rho = matrizrho(Phi,dL,dR,sp)
% MATRIRHO - Evalua la matriz rho del espacio indicado por sp

Ne = size(Phi,2);
p  = (1/Ne)*ones(Ne,1);

if sp == 1
    rho = zeros(dL);
    for k = 1:Ne
        PHI = reshape(Phi(:,k),[dR dL]).';
        rho = rho + p(k)*(PHI*PHI');
    end
else 
    rho = zeros(dR);
    for k = 1:Ne
        PHI = reshape(Phi(:,k),[dR dL]);
        rho = rho + p(k)*(PHI*PHI');
    end
end
if ~isequal(rho,rho')
    rho = 0.5*(rho+rho');
end
%=========================================================================


