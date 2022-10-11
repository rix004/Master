function IntensityMap(c,v,x)
logX = log10(x);
maxX = max(logX);
minX = min(logX);
z = (logX-minX)/(maxX-minX);
cmap = colormap(gray);
num = linspace(0,1,256);
col = zeros(numel(z),3);
for i1 = 1 : 3
    col(:,i1) = interp1(num, cmap(:,i1), z);
end
for i2 = 1 : numel(z)
    coords=v(c{i2},:);
    patch(coords(:,1),coords(:,2),col(i2,:));
    hold on
end