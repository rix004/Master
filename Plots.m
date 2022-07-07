% Plots
% IntensityMap(cells,vertices,p_darcy,'Pressure in Darcy domain');
% axis(D)
% 
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
% for i = 1:size(edges,1)
%     x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
%     y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
%     pressure = [p_network(edges(i,2)) p_network(edges(i,3))];
%     plot3(x,y,pressure,'.-','Color',[0 0 0],'LineWidth',edges(i,4)*10)
%     hold on
%     grid on
% end
% p_map = autumn(size(nodes,1)+size(TNinfo,1)*3);
% network_map = p_map(size(TNinfo,1)*2:end,:);
% for i = 1:size(nodes,1)
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

% darcy_map = p_map(1:3*size(TNinfo,1),:);
% 
% [xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
% p_points = [TNinfo(:,1:2);boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
% p_values = [p_darcy;Bv_darcy*ones(size(boundary_cells,1)+4,1)];
% vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
% surf(xq,yq,vq)
% colormap(darcy_map)

%     figure(15)
%     a = num2str(iter);
%     subplot(2,5,iter)
%     IntensityMap(cells,vertices,kT(:,iter),a)
%     plot(nodes(MacroTermIndexes(TN),1),nodes(MacroTermIndexes(TN),2),'r*')
%     hold on
%     title(a)
%     axis(D)

figure
loglog([0.25 0.25/2 0.25/4 0.25/8],[0.1 0.05 0.025 0.025/2],'-','LineWidth',3,'Color','r')
hold on
loglog([0.2408 0.1796 0.1231 0.0885 0.0615 0.0437 0.0310],[0.0991 0.059 0.0463 0.0431 0.0399 0.0249 0.0136],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2408 0.1747 0.1262 0.0877 0.0618 0.0439 0.0310],[0.0656 0.081 0.0495 0.0586 0.0356 0.0245 0.0164],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2408 0.1747 0.1262 0.0877 0.0618 0.0439 0.0311],[0.1136 0.0576 0.0864 0.0597 0.0413 0.0218 0.0170],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2462 0.1721 0.1236 0.0868 0.0622 0.0438 0.0311],[0.1703  0.0557  0.0700 0.0532 0.0249 0.0266 0.0155],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2481 0.1684 0.1262 0.0875 0.0627 0.0439 0.0311],[0.1001 0.1056 0.0667 0.0310 0.0292 0.0219 0.0150],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2390 0.1741 0.1233 0.0870 0.0620 0.0439 0.0310],[0.0761 0.0544 0.0757 0.0333 0.0327 0.0181 0.0222],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2425 0.1754 0.1238 0.0875 0.0622 0.0441 0.0311],[0.0866 0.1182 0.0607 0.0449 0.0239 0.0226 0.0232],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2374 0.1728 0.1245 0.0886 0.0622 0.0438 0.0310],[0.0885 0.0722 0.0684 0.0437 0.0297 0.0207 0.0146],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2481 0.1728 0.1208 0.0861 0.0620 0.0440 0.0312],[0.0755 0.0748 0.0503 0.0371 0.0239 0.0281 0.0146],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on
loglog([0.2425 0.1768 0.1208 0.0875 0.0618 0.0439 0.0309],[0.0975 0.0521 0.0615 0.0450 0.0268 0.0329 0.0180],'-.','Color',[0    0.5686    0.7157],'LineWidth',1.5)
hold on

L2_matrix = [0.0656 0.081 0.0495 0.0586 0.0356 0.0245 0.0164;
            0.1136 0.0576 0.0864 0.0597 0.0413 0.0218 0.0170;
            0.1703  0.0557  0.0700 0.0532 0.0249 0.0266 0.0155;
            0.1001 0.1056 0.0667 0.0310 0.0292 0.0219 0.0150;
            0.0991 0.059 0.0463 0.0431 0.0399 0.0249 0.0136;
            0.0761 0.0544 0.0757 0.0333 0.0327 0.0181 0.0222;
            0.0866 0.1182 0.0607 0.0449 0.0239 0.0226 0.0232;
            0.0885 0.0722 0.0684 0.0437 0.0297 0.0207 0.0146;
            0.0755 0.0748 0.0503 0.0371 0.0239 0.0281 0.0146;
            0.0975 0.0521 0.0615 0.0450 0.0268 0.0329 0.0180];
h_matrix = [0.2408 0.1747 0.1262 0.0877 0.0618 0.0439 0.0310;
            0.2408 0.1747 0.1262 0.0877 0.0618 0.0439 0.0311;
            0.2462 0.1721 0.1236 0.0868 0.0622 0.0438 0.0311;
            0.2481 0.1684 0.1262 0.0875 0.0627 0.0439 0.0311;
            0.2408 0.1796 0.1231 0.0885 0.0615 0.0437 0.0310;
            0.2390 0.1741 0.1233 0.0870 0.0620 0.0439 0.0310;
            0.2425 0.1754 0.1238 0.0875 0.0622 0.0441 0.0311;
            0.2374 0.1728 0.1245 0.0886 0.0622 0.0438 0.0310;
            0.2481 0.1728 0.1208 0.0861 0.0620 0.0440 0.0312;
            0.2425 0.1768 0.1208 0.0875 0.0618 0.0439 0.0309];

mean_vect = zeros(1,size(L2_matrix,2));
mean_h = zeros(1,size(L2_matrix,2));
for iter4 = 1:size(L2_matrix,2)
    mean_vect(iter4)=mean(L2_matrix(:,iter4));
    mean_h(iter4)=mean(h_matrix(:,iter4));
end

std_dev = std(L2_matrix);
errorbar(mean_h,mean_vect,std_dev,'x-','LineWidth',2.5,'Color',[0    0.1882    0.9059])
hold on

xlabel('h','FontSize',18)
ylabel('Error','FontSize',18)
lgd = legend('1st order reference','','','','','','','','','','Error','Average error');
lgd.FontSize=(14);
legend('Location','northwest')








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

figure('Name','Pressure plot')
[xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
p_points = [cell_center(ok_cells,1:2);boundary_cells(:,1:2)];
p_values = [p_darcy(ok_cells);bv];
pEx_values = [p_exact(cell_center(ok_cells,1),cell_center(ok_cells,2));bv];
vq = griddata(p_points(:,1),p_points(:,2),p_values-pEx_values,xq,yq);
surf(xq,yq,vq)