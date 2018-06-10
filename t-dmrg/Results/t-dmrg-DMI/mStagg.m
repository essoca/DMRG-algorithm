
load D100FE
t_steps = 90;
%t_steps = size(Szprof,2);
ms = zeros(t_steps,2);
for k = 1:t_steps
    pSz = Szprof(:,k); N = size(Szprof,1);
    mags = (1/N)*sum(((-1).^(1:N)).*(pSz.'));
    ms(k,1:2) = [(k-1),mags];
end
plot(ms(:,1),ms(:,2),'r')
clear

