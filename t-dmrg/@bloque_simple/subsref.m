
function val = subsref(B,S)

switch S(1).type % tipo de subscript
    case '()' 
        error('() No esta soportado como subscript')
    case '{}'
        switch S(2).subs
            case 'l'
                val = B(S(1).subs{:}).l;
            case 'dim'
                val = B(S(1).subs{:}).dim;
            case 'O'
                val = B(S(1).subs{:}).O;
            otherwise
                error('Estas utilizando un campo invalido')
        end
    case '.'
        switch S(1).subs % valor del subscript
            case 'l'
                val = B.l;
            case 'dim'
                val = B.dim;
            case 'O'
                val = B.O;
            otherwise
                error('Estas utilizando un campo invalido')
        end
end