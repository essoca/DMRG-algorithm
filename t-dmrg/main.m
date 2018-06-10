
function main()

% Parametros del sistema
global b0 b1 
S1      = 0.5;
S0      = 0.5;
celdaU  = [1 0];
Nrp     = 100;
% Parametros del modelo
% global Jxy Jz 
% h_int = 'hXXZ';
% Jxy   = 1;
% Jz    = 0;
global Dx Dy Dz J
h_int = 'hDM';
J  = 0;
Dx = 1;
Dy = 1;
Dz = 0;
% Hamiltoniano inicial
global h
Bn = 2;
% hL = repmat(Bn,[1,Nrp]);
% hR = repmat(-Bn,[1,Nrp]);
% h  = [hL,hR]; % Domain wall
%h = repmat([Bn,-Bn],[1,2*Nrp]); % Neel state
h = repmat(-Bn,[1,2*Nrp]);
% Parametros del DMRG
global Ne dispR 
m       = [20,40,60];
Nsweeps = numel(m);
Ne      = 1;
dispR   = 1;
% Parametros del t-DMRG
tf = 100; 
dt = 0.05;

clc

% Sitios b0 y b1 iniciales
[Sz1 Sp1] = matrices_spin(S1); d1 = 2*S1+1;
[Sz0 Sp0] = matrices_spin(S0); d0 = 2*S0+1;
b1 = bloque(1,d1,zeros(d1),Sp1,Sz1,eye(d1));
b0 = bloque(1,d0,zeros(d0),Sp0,Sz0,eye(d0));

[Bl,Br,d,Phi] = gsdmrg(S0,S1,celdaU,Nrp,m,Nsweeps,'h_ini');

% No se necesitan mas operadores de bloques -> construir bloques simples
for l = 1:length(Bl)
    Bl{l} = bloque_simple(Bl{l}.l, Bl{l}.dim, Bl{l}.O);
    Br{l} = bloque_simple(Br{l}.l, Br{l}.dim, Br{l}.O);
end

% Empezar evolucion temporal
tdmrg(Phi,h_int,Bl,Br,d,m,tf,dt)