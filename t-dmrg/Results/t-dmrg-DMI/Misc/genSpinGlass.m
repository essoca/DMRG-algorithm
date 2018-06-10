
function h = genSpinGlass(Bn,L)
% Genera un vector h con L entradas Bn o -Bn aleatorias

P  = [0.5,0.5];
vd = [Bn,-Bn];
h  = zeros(L,1);
for i = 1:L
    h(i) = genRandom(vd,P);
end

function xout = genRandom(vd,P)
% Muestrea los numeros en el vector vd con probabilidad P

ld = length(vd);
pd = zeros(ld-1,1); pd(1) = P(1);
for k = 1:ld-2
    pd(k+1) = pd(k)+P(k+1);
end

r = rand;
if r < pd(1)
    xout = vd(1);
    return
else
    for k = 1:ld-2
        if r > pd(k) && r < pd(k+1)
            xout = vd(k+1);
            return
        end
    end
end
xout = vd(ld);
