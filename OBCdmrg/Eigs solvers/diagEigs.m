
function [Phi,E] = diagEigs(varargin)

global Ne

Ble   = varargin{1};
Bre   = varargin{2};
h_int = varargin{3};

% Opciones para ARPACK 
%opts.isreal = 0;              % <---
opts.disp  = 0;
opts.tol   = 1e-12; 
opts.maxit = 5000; 
if nargin == 4
    opts.v0 = varargin{4};
end
% Diagonalizar con ARPACK
Hv = @(v) Hfun(v,Ble,Bre,h_int);
[Phi,E] = eigs(Hv,Ble.dim*Bre.dim,Ne,'SR',opts);

%=========================================================================
function w = Hfun(v,Ble,Bre,h_int)
% HFUN - Realiza la multiplicacion matriz-vector

w0 = reshape(v,[Bre.dim,Ble.dim]).';
w1 = Ble.HB * w0 + w0 * Bre.HB.' + feval(h_int,Ble,Bre,w0);                           
w  = real(reshape(w1.',[Ble.dim*Bre.dim,1]));
%w  = reshape(w1.',[Ble.dim*Bre.dim,1]);           %<----
%=========================================================================
