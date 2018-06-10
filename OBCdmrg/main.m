
function main()

% Parametros del sistema
global b0 b1 S1 S0
S1      = 0.5;
S0      = 0.5;
celdaU  = [1 0];
Nrp     = 50; 
% Parametros del modelo
h_int = 'hXXZ';
global Jxy Jz
Jxy = 1; 
Jz = 0; 
%h   = 1e-5;  
% Parametros del DMRG
global Ne dispR
err     =  [1e-6,1e-6,1e-6,1e-8,1e-8];
Nsweeps = numel(err);
Ne      = 1;
dispR   = 1;

clc
% Imprimir parametros usados
fprintf('\n\n')
fprintf('\t\t\t DMRG cadena de espines mixtos (%.1f,%.1f)\n',S1,S0)
fprintf('\t\t\t E. Solano-Carrillo, Fecha: %s\n\n',date)
fprintf('Parametros usados:\n\n')
fprintf('Tamanho de la cadena\nL = %d sitios\n',numel(celdaU)*Nrp)
fprintf('Hamiltoniano\nJxy = %.8f, Jz = %.8f\n',Jxy,Jz)
fprintf('Autoestados de rho a retener en %d sweeps\n',Nsweeps)
disp(['err = [',num2str(err),']'])
fprintf('\n')

% Sitios b0 y b1 iniciales
[Sz0 Sm0] = matrices_spin(S0);
[Sz1 Sm1] = matrices_spin(S1);
b0 = bloque(1,2*S0+1,zeros(2*S0+1),Sm0,Sz0,eye(2*S0+1));
b1 = bloque(1,2*S1+1,zeros(2*S1+1),Sm1,Sz1,eye(2*S1+1));

[profSz,E,uPm] = gsdmrg(S0,S1,celdaU,Nrp,err,Nsweeps,h_int);
%save profSz profSz

% % Campo magnetico
% h = 0.05:0.1:4.05; lh = length(h);
% % Anisotropia magnetica
% D1 = 0:0.1:4; lD1 = length(D1);
% D0 = 0;
% 
% hm  = zeros(lh,2,lD1);
% psz = zeros(numel(celdaU)*Nrp,lh,lD1);
% for k = 1:lD1
%     fprintf('\n\n D1 = %.2f\n\n',D1(k))
%     for n = 1:lh
%         b0.HB = -h(n)*Sz0 + D0*(Sz0^2);
%         b1.HB = -h(n)*Sz1 + D1(k)*(Sz1^2);
%         [profSz,E,uPm] = gsdmrg(S0,S1,celdaU,Nrp,err,Nsweeps,h_int);
%         fprintf('-> h = %.2f   E = %.8f    m = %d\n',h(n),E,uPm)
%         hm(n,1:2,k) = [h(n),sum(profSz)/(2*Nrp)];
%         psz(:,n,k)  = profSz;
%     end
%     save hm hm
%     save psz psz
% end


