function IntensityMap(c,v,x,Name)
    x_sorted = sort(x,'descend');
    cm = colormap(gray(length(x)));
    figure('Name',Name)
    for i = 1:size(c,1)
        coords=v(c{i},:);
        pgon = polyshape(coords(:,1),coords(:,2));
        pg = plot(pgon);
        x_here = x(i);
        ind = find(x_sorted==x_here);
        ind = ind(1);
        pg.FaceColor = cm(ind,:);
        hold on
    end
end