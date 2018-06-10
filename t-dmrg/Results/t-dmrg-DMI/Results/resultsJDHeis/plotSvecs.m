
clc, clear all
load r9_7JDz.mat
dw = 1:10;
t = 2; 
% constante de red
a = 1;
% Posiciones de los espines
x = (dw(1):a:a*numel(dw)).'; y = zeros(numel(dw),1); 
u = Sxprof(dw,t+1); v = Szprof(dw,t+1); 
pos = [x,y]; vec = [u,v];
arrow3(pos,pos+vec,'1.5_r',0.9)

te= atan(Szprof(dw,t+1)./Sxprof(dw,t+1));
for k = 1:length(te)-1    
    teap(k,1:3) = [k,te(k),te(k)-te(k+1)];   
end 
teap

