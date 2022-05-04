function[O]=Overlap(cells,vertices,grid_elements,grid_vertices)
for i = 1:size(cells,1)
%     NeighbourNodes = Tn;
%     if size(Tn,1)>1
%         NeighbourNodes(i,:)=[];
%         [~,~,~,rad]=NearestNode(Tn(i,1),Tn(i,2),NeighbourNodes);
%         rad = rad/3*2;
%     elseif size(Tn,1)==1
%         rad = 1;
%     end
%     x = Tn(i,1);
%     y = Tn(i,2);
%     th = 0:pi/50:2*pi;
%     xunit = rad *cos(th) + x;
%     yunit = rad *sin(th) + y;
%     xcoords = xunit(1:end-1);
%     ycoords = yunit(1:end-1);
%     B(i) = polyshape(xcoords,ycoords);
    poly_coords = vertices(cells{i},:);
    pgon(i) = polyshape(poly_coords(:,1),poly_coords(:,2));
    for j = 1:size(grid_elements,1)
        coords=grid_vertices(grid_elements{j},:);
        triangle(j) = polyshape(coords(:,1),coords(:,2));
        O(i,j)=area(intersect(pgon(i),triangle(j)));
    end
    figure(3)
    Bplot = plot(pgon(i));
    Bplot.FaceColor = [0.8500 0.3250 0.0980];
    Bplot.EdgeColor = 'b';
    Bplot.LineWidth = 1;
    hold on
    axis equal
end
for k = 1:size(grid_elements,1)
    PolyPlot = plot(triangle(k));
    PolyPlot.EdgeColor = [0.8500 0.3250 0.0980];
    PolyPlot.FaceColor = [1 1 1];

%     [xp,yp]=centroid(pgon(k));
%     text(xp,yp,[num2str(k)])
    hold on
    axis equal
end
    
O = sparse(O');
end