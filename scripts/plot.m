%% First read the ply file
pc = pcread("I30.ply");
pos = pc.Location;
px = pos(pos(1:end,1)>5,1:end);
curvature = tdfread("c.out"," ");
cx = curvature.Average(pos(1:end,1)>5);
curve = abs(curvature.Average);
coeff = tdfread("coeff.out"," ");
dir=tdfread("dir.out"," ");
%% Plot principal dir zy
value = sqrt(dir.p1y.^2 + dir.p1z.^2);
pcshow(pos, value)
colorbar;
caxis([0,1])
colormap jet
%% Plot principal dir x
value = abs(dir.p1x);
pcshow(pos, value)
colorbar;
caxis([0,1])
colormap jet
%% Plot curvature 
value = abs(curve);
pcshow(px, cx)
colorbar;
caxis([0,0.5])
colormap jet
%% Plot derivative of curvature in principal dir 1
value = coeff.coefficient3;
pcshow(pos, value)
colorbar;
caxis([0,0.5])
colormap jet