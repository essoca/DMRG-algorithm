
function crossovCorre()

n = [10,6,9,11,10,8,7,17,19,9,17,13,14,3,19]; 
m = [9,5,7,8,7,5,4,9,9,4,7,5,5,1,6];
D = sqrt((n./m).^2-1);
t = [100,80:-5:15];

for i = 1:length(n)
    load(strcat('r',int2str(n(i)),'_',int2str(m(i)),'JHDz.mat'))
    subplot(2,1,1), plot(SxSx(:,t(i)),'ko-')
    subplot(2,1,2), plot(SxSxr(:,t(i)),'ro-') 
end



plot(mxx(:,1),mxx(:,2),'ko-','MarkerFaceColor','k','MarkerSize',5)
hold on
plot(mxxr(:,1),mxxr(:,2),'ro-','MarkerFaceColor','r','MarkerSize',5)