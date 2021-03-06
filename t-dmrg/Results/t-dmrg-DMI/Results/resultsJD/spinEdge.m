
function spinEdge()

% clc, clear all
% D = 1.2;
% load(strcat('D',num2str(D),'-00FE.mat'))
% tf = 20; t = 0:0.1:tf;
% % Comparacion con dinamica de espines clasica
% Dc = 0.25*D; Dm = 0.5*D;
% Sz1 = 0.5*(1-(tanh(Dc*t)).^2)./(1+(tanh(Dc*t)).^2);
% Sz11 = 0.5+0.5*(1-(Dm*t+1).^2)./(1+(Dm*t+1).^2);
% %subplot(4,1,4), plot(0:tf,abs(Szprof(1,1:tf+1)),'ko-','LineWidth',0.2,'MarkerFaceColor','k')
% plot(0:tf,Szprof(1,1:tf+1),'ko-','MarkerFaceColor','k')
% hold on
% %plot(t,Sz1.*cos(Dm*t),'k','LineWidth',1.5,'LineStyle','--')
% plot(t,Sz1,'c','LineWidth',1.5)
% plot(t,Sz11,'r','LineWidth',1.5)


clc, clear all
Tu = @(m,n) (2^(p(m,n)+2))*m*pi;
Td = @(D) 2*pi./sqrt((0.5*D).^2+0.0625);
%m = 1:2:9; n = 3:2:11; r = m./n;
%m = 2:2:6; n = 3:2:7; r = m./n;
%m = [5,3,7,4,9,5,11]; n = [9,8,11,7,13,8,17]; r = m./n;
%m = [1,1,3,2,3]; n = [7,8,19,13,20]; r = m./n;
%m = [1,1]; n = [11,12]; r = m./n;
m = [3,5,1,1,5,4,11,3,1,1,3,2,3,3,3,2,3,3,3,3,5]; 
n = [5,9,2,3,8,7,17,8,7,8,19,13,20,16,17,15,13,14,7,11,19]; 
r = m./n; D = 0.5*sqrt(1-r.^2)./r; [D,orden]=sort(D);
m = m(orden); n = n(orden);
for k = 1:length(m)
    TJ(k,1:3) = [D(k),2*PeriodoJ(m(k),n(k),D(k))];
    %Ttd(k,1:2) = [D(k),Td(D(k))];
    %Ttu(k,1:2) = [D(k),Tu(m(k),n(k))];
end
errorbar(TJ(:,1),TJ(:,2),TJ(:,3),'ko-','MarkerFaceColor','k')
%plot(TJ(:,1),TJ(:,2),'ko:','MarkerFaceColor','k')
hold on
Dt = 0.5:0.01:4; Ttd = [Dt',Td(Dt)'];
plot(Ttd(:,1),Ttd(:,2),'r','LineWidth',1)
% Teorico J/D >> 1
m = [11,4,5,1,1]; n = [17,7,9,2,3];
r = m./n; D = 0.5*sqrt(1-r.^2)./r; [D,orden]=sort(D);
m = m(orden); n = n(orden);
for k = 1:length(m)
    Ttu(k,1:2) = [D(k),Tu(m(k),n(k))];
end
plot(Ttu(:,1),Ttu(:,2),'ko','MarkerFaceColor','b')
plot(Dt,2*pi./Dt,'c-')
plot(Dt,4*pi./Dt,'g-')

%=============================================================
function T = PeriodoJ(m,n,D)
% Mide periodo (temporal o espacial) de la curva Sz

load(strcat('r',int2str(m),'_',int2str(n),'D',num2str(D),'.mat'))
Sz = Szprof(1,1:80);
%plot(Sz,'ko-')
% Derivada numerica
for k = 1:length(Sz)-1
    dSz(k) = Sz(k+1)-Sz(k);
end
% Rastrear cambios de signo para hallar minimos y maximos
c=0;
for k = 1:length(dSz)-1
    if dSz(k+1)*dSz(k) < 0
        c = c+1;
        ind(c) = k+1;
    end
end
if ~exist('ind','var'), T=[]; return; end
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
function vp = p(m,n)
% Funcion paridad simultanea de m y n

if (rem(n,2)==0 && rem(m,2)==0)||(rem(n+1,2)==0 && rem(m+1,2)==0)
    vp = 0;
else
    vp = 1;
end
%=============================================================