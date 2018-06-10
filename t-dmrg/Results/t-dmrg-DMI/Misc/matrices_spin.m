%-------------------------------------------------------------------------
%     Calcula las representaciones matriciales de los operadores de spin 
%     en el sistema.
%-------------------------------------------------------------------------
%  Sintaxis:
%     [Sp Sz] = matrices_spin(J);
%  Entrada:
%     S  : número cuántico de spin
%  Salida:
%     Sp : operador de subida -> S+ = Sx + iSy. Notar que S- = Sp'
%     Sz : componente z del operador de spin

function [Sz Sp] = matrices_spin(S)

N = 2*S+1;
I = 1:N-1;
J = 2:N;
Sz = sparse(1:N, 1:N, S-(0:2*S), N, N);
Sp = sparse(I, J, sqrt((N-I).*I), N, N);
%  m = -S:S-1;
%  Sp2 = sparse(I, J, sqrt(S*(S+1)-m.*(m+1)), N, N);
