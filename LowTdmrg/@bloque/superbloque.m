
function [Emin,Emax] = superbloque(varargin)
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

global Ne

Ble   = varargin{1};
Bre   = varargin{2};
h_int = varargin{3};

% Opciones para IRBLEIGS
opts.k      = Ne;
opts.funpar = {Ble,Bre,h_int};
opts.nbls   = 3;
opts.maxit  = 200;
opts.blsz   = 1;
opts.sigma  = 'se';
opts.tol    = 1e-8;
if nargin == 4
    opts.v0 = varargin{4};
end

% Dimension del superbloque
n = Ble.dim*Bre.dim;

conv = 1; ct = 0;
while conv(1)
    % Diagonalizar 'implicitly restarted block-Lanczos'
    [Phid,Emin,conv] = irbleigs('Hv',n,opts);
    if conv(1)
        opts.v0 = Phid; ct = ct+1; 
        fprintf('Buscando Emin: %d veces\n',ct+1)
    else
        break
    end
    % Si no converge aun, duplicar numero de iteraciones
    if ct > 3
        fprintf('Refinando la busqueda: duplicar maxit.\n')
        if opts.maxit > 5000, return, end
        opts.maxit = 2*opts.maxit;
        opts.nbls = 2*opts.nbls;
    end
end

opts.tol = 1e-6;
opts.sigma = 'le';
conv = 1; ct = 0;
while conv(1)
    % Diagonalizar 'implicitly restarted block-Lanczos'
    [Phiu,Emax,conv] = irbleigs('Hv',n,opts);
    if conv(1)
        opts.v0 = Phid; ct = ct+1; 
        fprintf('Buscando Emax: %d veces\n',ct+1)
    else
        return
    end
    % Si no converge aun, duplicar numero de iteraciones
    if ct > 3
        fprintf('Refinando la busqueda: duplicar maxit.\n')
        if opts.maxit > 5000, return, end
        opts.maxit = 2*opts.maxit;
        opts.nbls = 2*opts.nbls;
    end
end


