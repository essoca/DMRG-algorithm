
function spinEdge()

clc, clear all
D = 1;
load(strcat('D',num2str(D),'-00FE.mat'))
tf = 30; t = 0:0.1:tf; 
subplot(2,1,1), plot(0:tf,Szprof(1,1:tf+1),'ko-','MarkerFaceColor','k','LineWidth',0.1)
%plot(0:tf,Szprof(1,1:tf+1),'ko','MarkerFaceColor','k','LineWidth',0.1)
hold on
Sz1 = 0.5*(1-(tanh(D*t/4)).^2)./(1+(tanh(D*t/4)).^2);
plot(t,Sz1,'g','LineWidth',1)
% Solucion clasica
trange = [t(1),t(end)];
theta0 = zeros(1,size(Szprof,1));
%theta0 = zeros(1,2);
[t,theta] = ode45('ecdif',trange,theta0);
plot(t,0.5*cos(theta(:,1)),'b')
% Segundo sitio
subplot(2,1,2), plot(0:tf,Szprof(2,1:tf+1),'ko-','MarkerFaceColor','k','LineWidth',0.1)
hold on
plot(t,0.5*cos(theta(:,2)),'b')

% clc, clear all
% D = 0.2:0.2:2.0;
% for m = 1:length(D)
%     TD(m,1:3) = [D(m),PeriodoD(D(m))];
% end
% axes('pos',[0.5 0.6 0.25 0.3]);
% errorbar(TD(:,1),TD(:,2),TD(:,3),'ko','MarkerFaceColor','k')
% Dt = 0.2:0.01:2;
% hold on, plot(Dt,2*pi./Dt,'k')

%=============================================================
function T = PeriodoD(D)
% Mide periodo (temporal o espacial) de la curva Sz

st = 3; %sitio
load(strcat('D',num2str(D),'-00FE.mat'))
%plot(Szprof(st,:),'go-')
% Derivada numerica
for k = 1:size(Szprof,2)-1    
    dSz(k) = Szprof(st,k+1)-Szprof(st,k);
end
% Rastrear cambios de signo para hallar minimos y maximos
n=0;
for k = 1:length(dSz)-1
    if dSz(k+1)*dSz(k) < 0
        n = n+1;
        ind(n) = k+1;
    end
end
% Escoger valles en la curva Szprof
tv = ind(1:2:end);
if length(tv) >= 2
    for k = 1:length(tv)-1
        % Hallar periodo
        tT(k) = tv(k+1)-tv(k);
    end
    T(1:2) = [mean(tT),std(tT)];
elseif length(ind) == 2
    T(1:2) = [mean(2*ind(1),ind(2)),std(2*ind(1),ind(2))];
elseif length(ind) == 1
    T(1:2) = [2*ind(1),0];
else
    error('No se puede definir el periodo')
end
%=============================================================
