
function densplot()

clc
yrange = 0:0.05:100;
%xrange = 40:160;
xrange = 1:200;
load SzDM

[X,Y] = meshgrid(xrange,yrange); szt = [];
for k = 1:length(yrange)
    szt = [szt;abs(Szprof(xrange,k)')];
end
colormap(gray(256))
% Large values will be near black and small ones will be near white 
mapc = colormap; smp = size(mapc,1); ct = 0;
for k = smp:-1:1
    ct = ct+1;
    map(ct,1:3) = mapc(k,:);
end
colormap(map)
pcolor(X,Y,szt)
shading interp 