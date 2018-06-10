
function diagExacta()

global Dx Dy Dz J
Dx = 0; 
Dy = 0;
Dz = 1;
J  = 0;

Sx = 0.5*[0,1;1,0];
Sy = 0.5*[0,-i;i,0];
Sz = 0.5*[1,0;0,-1];
h = @(x,y) hDM(Sx,Sy,Sz);
% Numero de sitios
L = 6; clc
% Sectores de spin
ss = spinsectors(L);
% Hamiltoniano y estado inicial
[H,Psi0] = Hamilt(L,h);
[Phi,E] = eig(full(H)); 

% Evolucionar estado inicial
dt = 0.01; tf = 3; tm = 1; 
ntm = tf/tm+1;
tr = 0:dt:tf;
Psi = Psi0; k=0;
lt = length(tr);
sx = zeros(L,ntm);
sy = zeros(L,ntm);
sz = zeros(L,ntm);
sxa = zeros(L,ntm);
sya = zeros(L,ntm);
sza = zeros(L,ntm);
for n = 1:lt
    t = (n-1)*dt;    
    Psi = expm(-i*H*dt)*Psi;
    if rem(t,tm) == 0
        k = k+1;
        fprintf('t = %d\n',t)
        % Medir perfiles de mag
        [sx(:,k),sy(:,k),sz(:,k)] = vEspe(Psi,Sx,Sy,Sz,L);
        [sxa(:,k),sya(:,k),sza(:,k)] = vEspeA(t,Sx,Sy,Sz,L);
    end    
end
save profs sx sy sz sxa sya sza

% t = 0:0.01:30; L = 200;
% for k = 1:length(t)
%     [sxa(:,k),sya(:,k),sza(:,k)] = vEspeA(t(k),Sx,Sy,Sz,L);
% end
% save sa sxa sya sza

%load D001FE.mat
%load sa
%subplot(1,2,2), plot(t,Sxprof(1,:),'k','MarkerFaceColor','k','MarkerSize',2)
%hold on
%plot(t,sxa(1,:),'b','MarkerFaceColor','b','MarkerSize',2)

% %Graficar
% clc
% load profs
% subplot(1,2,1), plot(0:6,sx(1,:),'ko-','MarkerFaceColor','k')
% hold on
% plot(0:6,sxa(1,:),'ko-','MarkerFaceColor','r')


% % Graficar
% tp = 3; clc
% load profs
% subplot(2,2,3), plot(sx(:,tp+1),'ko-','MarkerFaceColor','k')
% hold on
% plot(sy(:,tp+1),'ro-','MarkerFaceColor','r')
% plot(sz(:,tp+1),'go-','MarkerFaceColor','g')
% subplot(2,2,4), plot(sxa(:,tp+1),'ko-','MarkerFaceColor','k')
% hold on
% plot(sya(:,tp+1),'ro-','MarkerFaceColor','r')
% plot(sza(:,tp+1),'go-','MarkerFaceColor','g')

E = diag(E); %E = round(2.*E).*0.5;
% Producir lista de E sin multiplicidad
x = [];
y = E;
while (~isempty(y))
  x(end+1) = max(y);
  y = y(y~=max(y));
end
% Valores de E con degeneracy
sE = zeros(length(x),2);
for k = 1:length(x)
    dS = numel(find(x(k)==E));
    sE(k,1:2) = [x(k),dS];
end
sE

%=======================================================================
function [sxa,sya,sza] = vEspeA(t,Sx,Sy,Sz,L)
% Valores esperados con Psi propuesto como ansatz

W = @(tau,lor) getRho(tau,Sx,Sz,lor);
sxa = zeros(L,1);
sya = zeros(L,1);
sza = zeros(L,1);
for j = 1:0.5*L   
    tj = j-1;
    sxa(j) = real(trace(Sx*W(t-tj,'l')));
    sya(j) = real(trace(Sy*W(t-tj,'l')));
    sza(j) = real(trace(Sz*W(t-tj,'l')));
    sxa(L-j+1) = real(trace(Sx*W(t-tj,'r')));
    sya(L-j+1) = real(trace(Sy*W(t-tj,'r')));
    sza(L-j+1) = real(trace(Sz*W(t-tj,'r')));
end
%=======================================================================
function rho = getRho(tau,Sx,Sz,lor)
% Matriz densidad estatica o dinamica de un sitio

rhos = 0.5*eye(2)+Sx; D = 1;
switch lor
    case 'l'
        rhod = 0.5*eye(2)+cos(D*tau/2)*Sx+0.5*sin(D*tau)*Sz;
    case 'r'
        rhod = 0.5*eye(2)+cos(D*tau/2)*Sx-0.5*sin(D*tau)*Sz;
end
if tau<0, rho = rhos; else rho = rhod; end
%=======================================================================
function [sx,sy,sz] = vEspe(Psi,Sx,Sy,Sz,L)
% Valores esperados con Psi exacto

sx = zeros(L,1);
sy = zeros(L,1);
sz = zeros(L,1);
for k = 1:L
    nl = k-1; nr = L-(nl+1);
    sx(k) = real(Psi'*kron(kron(speye(2^nl),Sx),speye(2^nr))*Psi);
    sy(k) = real(Psi'*kron(kron(speye(2^nl),Sy),speye(2^nr))*Psi);
    sz(k) = real(Psi'*kron(kron(speye(2^nl),Sz),speye(2^nr))*Psi);
end
%=======================================================================
function ss = spinsectors(L)
% Sectores de spin (con dimension) para L sitios

S = 0.5*L;
s = 0:S; ss = zeros(numel(s),2);
for k = 1:numel(s)
    ss(k,1:2) = [s(k),2*(S-s(k)+1)];
end
%=======================================================================
function [H,vx] = Hamilt(L,h)
% Hamiltoniano de L sitios

dH = 2^L;
H = sparse(zeros(dH));
for k = 1:L-1
    nl = k-1; nr = L-(nl+2);
    H = H + kron(kron(speye(2^nl),h(k,k+1)),speye(2^nr));
end
% Estado inicial ->->->
vx = ones(dH,1)/norm(ones(dH,1));
%=======================================================================
function h = hDM(Sx,Sy,Sz)
% Hamiltoniano de 2 sitios

global Dx Dy Dz J

% Productos entre dos sitios
SySz = kron(Sy,Sz); SzSy = kron(Sz,Sy);
SzSx = kron(Sz,Sx); SxSz = kron(Sx,Sz);
SxSy = kron(Sx,Sy); SySx = kron(Sy,Sx);
SzSz = kron(Sz,Sz); SxSx = kron(Sx,Sx);
SySy = kron(Sy,Sy);
% Interaccion D-M 
h = J*(SxSx+SySy+SzSz) + Dx*(SySz-SzSy)+Dy*(SzSx-SxSz)+Dz*(SxSy-SySx);
%=======================================================================