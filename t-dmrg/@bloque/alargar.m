
function Be = alargar(B,b,haciaDonde,h_int)
% ALARGAR - Alarga un bloque DMRG al agregar un sitio
%-----------------------------------------------------------------------
%     Alarga un bloque existente mediante la conexion con un nuevo sitio.
%     La operacion se simboliza como B(l,m)B(1,d) --> Ble(l+1,mxd) si se 
%     alarga el bloque izquierdo o   B(1,d)B(l,m) --> Bre(l+1,mxd) si se
%     alarga el bloque derecho. Aqui, d representa la dimension del espacio
%     de Hilbert del sitio agregado.
%-----------------------------------------------------------------------
%  Sintaxis:
%     Ble = alargar(B,b,'haciaR',h_int);
%     Bre = alargar(B,b,'haciaL',h_int);  
%  Entrada:
%     B          : Bloque que se va a alargar
%     b          : Bloque (sitio) que se agrega
%     haciaDonde : Se alarga un bloque 'haciaL' o 'haciaR'
%     h_int      : Hamiltoniano de interaccion entre sitios vecinos (@)
%  Salida:
%     Be         : Bloque izquierdo (Ble) o derecho (Bre) alargado

IB = speye(B.dim);
Ib = speye(b.dim);

switch haciaDonde
    case 'haciaR'  % alarga bloque izquierdo hacia la derecha
        HB  = kron(B.HB,Ib) + kron(IB,b.HB) + feval(h_int,B,b);
        Spb = kron(IB,b.Spb);
        Szb = kron(IB,b.Szb);
    case 'haciaL'  % alarga bloque derecho hacia la izquierda
        HB  = kron(Ib,B.HB) + kron(b.HB,IB) + feval(h_int,b,B);
        Spb = kron(b.Spb,IB);
        Szb = kron(b.Szb,IB);
    otherwise
        error(' Alargar haciaL o haciaR ')
end

% Construimos el bloque alargado
Be = bloque(B.l+1, B.dim*b.dim, HB, Spb, Szb, []);

