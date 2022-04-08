% Plots

% IntensityMap(cells,vertices,p_tn,'Voronoi diagram: Pressure in terminal network nodes')
% plot(Tn(:,1),Tn(:,2),'b.')
% hold on
% plot(RootNode(1),RootNode(2),'r*');
% hold on
% axis(D)



% figure('Name','Pressure plot for Darcy domain')
% plot3(cell_center(:,1),cell_center(:,2),p_darcy,'.','MarkerSize',20);
% grid on
% xlabel('x','FontSize',15)
% ylabel('y','FontSize',15)
% zlabel('p','FontSize',15)




% figure('Name','Pressure plot')
% for i = 1:Ne
%     x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
%     y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
%     pressure = [p_network(edges(i,2)) p_network(edges(i,3))];
%     plot3(x,y,pressure,'.-','Color',[0 0 0],'LineWidth',edges(i,4)*10)
%     hold on
%     grid on
% end
% p_map = autumn(Nn+Ntn*3);
% network_map = p_map(Ntn*2:end,:);
% for i = 1:Nn
%     x = nodes(i,1);
%     y = nodes(i,2);
%     pressure = p_network(i);
%     p_sorted = sort(p_network);
%     ind = find(p_sorted==pressure);
%     ind = ind(1);
%     plot3(x,y,pressure,'.','Color',network_map(ind,:),'MarkerSize',20);
%     hold on
% end
% zlabel('Node pressure')
% 
% darcy_map = p_map(1:3*Ntn,:);
% 
% [xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
% p_points = [Tn(:,1:2);boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
% p_values = [p_darcy;Bv_darcy*ones(size(boundary_cells,1)+4,1)];
% vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
% surf(xq,yq,vq)
% colormap(darcy_map)



% IntensityMap(cells,vertices,p_darcy,'Pressure in Darcy domain');
% axis(D)




% figure('Name','Peaceman pressure function')
% nx=180;
% ny=180;
% points = GridGeneration(nx,ny,D);
% map = gray(2);
% for i = 1:size(cells,1)
%     coords=vertices(cells{i},:);
%     pgon = polyshape(coords(:,1),coords(:,2));
%     pg = plot(pgon);
%     pg.FaceColor=(map(end,:));
%     hold on
%     for j = 1:size(points,1)
%         in(j) = inpolygon(points(j,1),points(j,2),coords(:,1),coords(:,2));
%     end
%     p_in=points(in==1,:);
%     distToVessel = sqrt((p_in(:,1)-cell_center(i,1)).^2+(p_in(:,2)-cell_center(i,2)).^2);
%     p_pm=-darcy_source(i)*mu/(2*k*pi)*log(distToVessel/edges(Tn(i,3),4))+p_tn(i);
%     p_pm_sorted = sort(p_pm);
%     dp = p_pm_sorted(2)-p_pm_sorted(1);
%     p_axis=min(p_pm_sorted):dp:max(p_pm_sorted);
%     newmap = parula(length(p_axis));
%     for iter = 1:size(p_in,1)
%         p_here = p_pm(iter);
%         ind = find(abs(p_axis-p_here)==min(abs(p_axis-p_here)));
%         ind = ind(1);
%         plot(p_in(iter,1),p_in(iter,2),'.','Color',newmap(ind,:))
%         hold on
%     end
%     c=scircle1(cell_center(i,1),cell_center(i,2),0.2*sqrt(cell_area(i)));
%     plot(c(:,1),c(:,2),'--','LineWidth',0.4,'Color',newmap(1,:));
%     hold on
% end
% axis(D)
% 
% figure('Name','Peaceman pressure at 0.2*deltaX')
% for i = 1:size(cells,1)
%     coords=vertices(cells{i},:);
%     pgon = polyshape(coords(:,1),coords(:,2));
%     pg = plot(pgon);
%     pg.FaceColor=(map(end,:));
%     hold on
%     for j = 1:size(points,1)
%         in(j) = inpolygon(points(j,1),points(j,2),coords(:,1),coords(:,2));
%     end
%     p_in=points(in==1,:);
%     distToVessel = sqrt((p_in(:,1)-cell_center(i,1)).^2+(p_in(:,2)-cell_center(i,2)).^2);
%     p_pm=-darcy_source(i)*mu/(2*k*pi)*log(distToVessel/edges(Tn(i,3),4))+p_tn(i);
%     p_pm2 = -darcy_source(i)*mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4))+p_tn(i);
%     p_pm_sorted = sort(p_pm);
%     dp = p_pm_sorted(2)-p_pm_sorted(1);
%     p_axis=min(p_pm_sorted):dp:max(p_pm_sorted);
%     newmap = parula(length(p_axis));
%     ind2 = find(abs(p_axis-p_pm2)==min(abs(p_axis-p_pm2)));
%     ind2 = ind2(1);
%     for iter = 1:size(p_in,1)
%         plot(p_in(iter,1),p_in(iter,2),'.','Color',newmap(ind2,:))
%         hold on
%     end
% end
% plot(cell_center(:,1),cell_center(:,2),'.','MarkerSize',7,'Color',newmap(end,:));
% hold on
% axis(D)