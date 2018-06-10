
function m_new = check_m (wrho,m)
% CHECK_M - Garantiza la correcta seleccion de los autoestados de rho
%-----------------------------------------------------------------------
%    Determina si m es la suma de las multiplicidades de los autoestados
%    de rho que se seleccionan para la nueva base del espacio de Hilbert
%    efectivo en el siguiente paso del DMRG.
%-----------------------------------------------------------------------
%  Sintaxis:
%    m_new = check_m (wrho,m)
%  Entrada:
%    wrho  : Vector con los autovalores de rho ordenados descendentemente
%    m     : Numero de autoestados de rho a retener en el truncamiento
%  Salida:
%    m_new : Nuevo m corregido para incluir multiplicidades completas

k = 0; wrho(abs(wrho)<eps) = 0;
wm = wrho(m);
while m+k+1 <= length(wrho)
    wnext = wrho(m+k+1); 
    if wnext ~= 0 && abs(wnext-wm) <= eps
        k = k + 1;
    else
        break
    end
end
m_new = m+k;


