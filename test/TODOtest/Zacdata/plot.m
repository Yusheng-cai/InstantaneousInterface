%% read ply 
tri = pcread("zac.ply");
pos = tri.Location;
curve  = readtable("c.dat");
meanC = curve.Average;
scatter3(pos(1:end,1),pos(1:end,2), pos(1:end,3), 100,abs(meanC),"filled");
colormap jet;
colorbar;
caxis([0,1.5]);