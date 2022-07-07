function IntensityMap(c,v,x,Name)
    x_unique = unique(x);
    x_sorted = sort(x_unique,'descend');
    cm = colormap(gray(length(x_sorted)));
    %cm = cm(1:length(x_sorted),:);
    figure('Name',Name)
    for i = 1:size(c,1)
        coords=v(c{i},:);
        pgon = polyshape(coords(:,1),coords(:,2));
        pg = plot(pgon);
        x_here = x(i);
        ind = find(x_sorted==x_here);
        ind = ind(1);
        if x_here == 0
            pg.FaceColor = [1 1 1];
        else
            pg.FaceColor = cm(ind,:);
        end
        hold on
    end
end