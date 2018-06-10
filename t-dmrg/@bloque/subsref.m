%--------------------------------------------------------------------------
%  Define la operacion 'subscripted reference', por ejemplo, cuando
%  hacemos asignaciones como Hbloque = B.HB, o ElemArray = B(3) o  
%  o ElemCell = B{n}. Estas dos ultimas asignaciones no seran soportadas.
%--------------------------------------------------------------------------

function val = subsref(B,S)

switch S.subs % valor del subscript
    case 'l'
        val = B.l;
    case 'dim'
        val = B.dim;
    case 'HB'
        val = B.HB;   
    case 'Spb'
        val = B.Spb;
    case 'Szb'
        val = B.Szb;
    case 'O'
        val = B.O;
    otherwise
        error('Campo invalido')
end