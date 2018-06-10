
function nmUtiles()

clc, clear all
J = 1;
m = 1:20; n = 1:20; 
D = @(m,n) sqrt((n/m)^2-1); % J=1
T = @(m,n) (2^p(m,n))*m*pi/(0.5*J);

ct = 0;
for k = 1:numel(m)
    for l = 1:numel(m)
        if m(k)/n(l) < 1
            if gcd(m(k),n(l)) == 1
                if D(m(k),n(l)) < 1
                    ct = ct+1; 
                    ind(ct,1:4)=[k,l,D(m(k),n(l)),T(m(k),n(l))];
                end
            end
        end
    end
end

x = [];
y = ind;
ct = 0;
zr = 1000*ones(size(ind));
% Escoger valores de m y n en donde el menor D defina Tmn univocamente
[Dord,ord]=sort(ind(:,3));
ind = ind(ord,1:4);
pause(1)


function vp = p(m,n)
% Funcion paridad simultanea de m y n

if (rem(n,2)==0 && rem(m,2)==0)||(rem(n+1,2)==0 && rem(m+1,2)==0)
    vp = 0;
else
    vp = 1;
end