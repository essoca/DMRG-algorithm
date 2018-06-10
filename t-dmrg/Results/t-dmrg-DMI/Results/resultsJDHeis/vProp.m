
function vProp()

clc
tf = 150; n = 9; m = 7; 
load(strcat('r',int2str(n),'_',int2str(m),'JDz.mat'))
for k = 1:0.5*size(Sxprof,1)
    Sxk = Sxprof(k,1:tf);
    %ind = find(abs(Szk-0.5) < 0.01);
    ind = find(abs(Sxk-0.5) < 0.1);
    td(k,1:2) = [k,numel(ind)-1];
end
plot(td(:,1),td(:,2),'ko-')

% D = 0.2:0.2:2;
% v = [4.87,2.463,1.649,1.237,0.993,0.825,0.707,0.617,0.549,0.496];
% plot(D,v,'ko-')
