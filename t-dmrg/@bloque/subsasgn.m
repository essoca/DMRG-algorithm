
function B = subsasgn(B,S,val)

switch S(1).type % tipo de subscript
    case '()' 
        error('() no esta soportado')
    case '{}'
        switch S(2).subs
            case 'l'
                B(S(1).subs{:}).l   = val;
            case 'dim'
                B(S(1).subs{:}).dim = val;
            case 'HB'
                B(S(1).subs{:}).HB  = val;
            case 'Spb'
                B(S(1).subs{:}).Spb = val;
            case 'Szb'
                B(S(1).subs{:}).Sz  = val;
            case 'O'
                B(S(1).subs{:}).O   = val;
            otherwise
                error('Campo invalido')
        end
    case '.'
        switch S(1).subs % valor del subscript
            case 'l'
                B.l   = val;
            case 'dim'
                B.dim = val;
            case 'HB'
                B.HB  = val;
            case 'Spb'
                B.Spb = val;
            case 'Szb'
                B.Szb = val;
            case 'O'
                B.O   = val;
            otherwise
                error('Campo invalido')
        end
end