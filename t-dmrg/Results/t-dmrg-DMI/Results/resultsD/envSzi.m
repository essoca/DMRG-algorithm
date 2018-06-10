
function envSzi()

global D ts L
D = 1; ts = 1/D; t = 80; L = 200;
load(strcat('D',num2str(D),'-00FE.mat'))

for k = 1:L/2
    pf(k) = Sz(k,t);
    pf(L-k+1) = pf(k); 
end
% subplot(2,1,1), plot(Szprof(1:L,t),'ko-','MarkerFaceColor','k')
% hold on
% plot(Syprof(1:L,t),'ko-','MarkerFaceColor','g')
% plot(Sxprof(1:L,t),'ko-','MarkerFaceColor','r')
% plot(pf,'k--')
% subplot(2,1,2), plot(SzSz(:,t),'ko-','MarkerFaceColor','k')
plot(Szprof(1:L,t),'ko-','MarkerFaceColor','k')
hold on
plot(pf,'k--')
% Solucion clasica
trange = [0,80];
theta0 = zeros(1,size(Szprof,1));
[t,theta] = ode45('ecdif',trange,theta0);
plot(1:200,0.5*cos(theta(end,1:200)),'r')

%===========================================
function szi = Sz(i,t)
% Sz en cada sitio

global D ts 
ti = (i-1)*ts; 
if t < ti
    szi = 0.5;
else   
    c = D*(t-ti);
    szi = (4+c)^2/((4+c)^2+c^2)-0.5;    
end
%===========================================