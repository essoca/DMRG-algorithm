
clc, clear all
load D1-00FE.mat
dw = 1:10;
t = 2; 
% constante de red
a = 1;
% Posiciones de los espines
x = (dw(1):a:a*numel(dw)).'; y = zeros(numel(dw),1); 
u = Szprof(dw,t+1); v = Sxprof(dw,t+1); 
pos = [x,y]; vec = [u,v];
arrow3(pos,pos+vec,'1.5_r',0.9)

% dw = 1:20; t = 4; 
% te= atan(Sxprof(dw,t+1)./Szprof(dw,t+1));
% for k = 1:length(te)-1    
%     teap(k,1:3) = [k,te(k),te(k)-te(k+1)];   
% end 
% teap

