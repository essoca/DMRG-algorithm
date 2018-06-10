
function xim = getxim(xi,Emin,Emax,Ble,Bre,h_int,beta)

global Lt 

% Escalar Hamiltoniano
lambda = 0.5*(Emax+Emin);
wH = 0.7*(2/(Emax-Emin));
Hsv = @(v) Hs(v,Ble,Bre,h_int,lambda,wH);

% Definicion de funciones de Bessel esfericas modificadas de 1ra clase
ibess = @(l,x) sqrt(0.5*pi/x)*besseli(l+0.5,x);
% Crear el estado target:
xim = zeros(size(xi));
% Ejecutar relacion de recurrencia
PlHs = CoaliRecursion(xi,Hsv,Lt,Ble.dim*Bre.dim);
for l = 1:Lt+1
    xim = xim + (l-0.5)*ibess(l-1,-0.5*beta/wH)*PlHs(:,l);
end
xim = real(xim);
%=========================================================================
function PlHs = CoaliRecursion(xi,Hsv,Lt,dimSB)
% Relacion de recurrencia para calcular xim

% Polinomios de legendre regulados aplicados en xi
PlHs = zeros(dimSB,Lt+1);
PlHs(:,1) = xi;
PlHs(:,2) = Hsv(xi);
% Derivadas correspondientes aplicadas en xi
PplHs = zeros(dimSB,Lt+1);
PplHs(:,2) = xi;
for l = 2:Lt
    % Recurrencia en las derivadas
    PplHs(:,l+1) = (2*l-1)*PlHs(:,l) + PplHs(:,l-1);
    % Recurrencia principal
    cl = (2*l-1)/l;
    PlHs(:,l+1) = cl * Hsv(PlHs(:,l)) - (l-1)/l* PlHs(:,l-1) + ...
                  cl * (2*pi/Lt)^2 * PplHs(:,l);
end
%=========================================================================
function w = Hs(v,Ble,Bre,h_int,lambda,wH)
% HFUN - Realiza la multiplicacion matriz-vector

w0 = reshape(v,[Bre.dim,Ble.dim]).';
w1 = wH * (Ble.HB * w0 + w0 * Bre.HB.' + feval(h_int,Ble,Bre,w0));  
w2 = wH * lambda * w0;
w  = real(reshape((w1-w2).',[Ble.dim*Bre.dim,1]));
%=========================================================================