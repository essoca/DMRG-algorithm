
function w = Hv(v,n,blsz,Param)
% HFUN - Realiza la multiplicacion matriz-vector

Ble   = Param{1};
Bre   = Param{2};
h_int = Param{3};

w0 = reshape(v,[Bre.dim,Ble.dim]).';
w1 = Ble.HB * w0 + w0 * Bre.HB.' + feval(h_int,Ble,Bre,w0);                           
%w  = real(reshape(w1.',[n,1]));
w  = reshape(w1.',[n,1]);

