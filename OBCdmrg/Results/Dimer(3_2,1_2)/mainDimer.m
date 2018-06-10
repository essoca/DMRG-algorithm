
function mainDimer()

% Parametros del sistema
global b0 b1 
S1      = 1.;
S0      = 0.5;
celdaU  = [1 0];
Nrp     = 50; 
% Parametros de interaccion experimentales (55cm-1 -> 80K)
Jexp = 101.2; jexp = -2.8; Dexp = 0.11;
global dd
h_int = 'hdimer'; 
dd = jexp/Jexp; 
D = Dexp/Jexp; 
% Parametros del DMRG
global Ne dispR
m       =  [80];
Nsweeps = numel(m);
Ne      = 2;
dispR   = 1;
D=0;

clc
% Imprimir parametros usadosr
fprintf('\n\n')
fprintf('\t\t\t DMRG cadena de espines mixtos (%.1f,%.1f)\n',S1,S0)
fprintf('\t\t\t E. Solano-Carrillo, Fecha: %s\n\n',date)
fprintf('Parametros usados:\n\n')
fprintf('Tamanho de la cadena\nL = %d sitios\n',numel(celdaU)*Nrp)
fprintf('Hamiltoniano\nJ = %.1f\nj = %.1f\n',1+dd,1-dd)
fprintf('Autoestados de rho a retener en %d sweeps\n',Nsweeps)
disp(['m = [',num2str(m),']'])
fprintf('\n')

% Sitios b0 y b1 iniciales
[Sz0 Sp0] = matrices_spin(S0);
[Sz1 Sp1] = matrices_spin(S1);
b0 = bloque(1,2*S0+1,zeros(2*S0+1),Sp0,Sz0,eye(2*S0+1));
b1 = bloque(1,2*S1+1,D*(Sz1^2),Sp1,Sz1,eye(2*S1+1));

%Egap = gap(S0,S1,celdaU,Nrp,m,Nsweeps,h_int);
%Esectores(S0,S1,celdaU,Nrp,m,Nsweeps,h_int)
gsdmrg(S0,S1,celdaU,Nrp,m,Nsweeps,h_int);


% % Campo magnetico
% h = 0.05:0.2:8.05; lh = length(h);
% % Anisotropia magnetica
% D1 = 1:0.1:3; lD1 = length(D1);
% D0 = 0;
% 
% hm  = zeros(lh,2,lD1);
% he  = zeros(lh,2,lD1);
% psz = zeros(numel(celdaU)*Nrp,lh,lD1);
% for k = 1:lD1
%     fprintf('\n\n D1 = %.2f\n\n',D1(k))
%     for n = 1:lh
%         b0.HB = -h(n)*Sz0 + D0*(Sz0^2);
%         b1.HB = -h(n)*Sz1 + D1(k)*(Sz1^2);
%         [profSz,E,uPm] = gsdmrg(S0,S1,celdaU,Nrp,m,Nsweeps,h_int);
%         fprintf('-> h = %.2f   E = %.8f    1-Pm = %.0E\n',h(n),E,uPm)
%         hm(n,1:2,k) = [h(n),sum(profSz(50:150))/(Nrp+1)];
%         he(n,1:2,k) = [h(n),E/(2*Nrp)];
%         psz(:,n,k)  = profSz;
%         save HeisD1/hm hm
%         save HeisD1/he he
%         save HeisD1/psz psz
%     end
% end

