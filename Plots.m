function[] = Plots(ThisFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TPFA TEST VORONOI GRID %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure, 'TPFAConvergenceVoronoi')
    % Load data
    indata = load('L2errorData');
    L2_error = indata.L2_error;
    h = indata.h;

    % Make figure
    figure

    % Reference line
    loglog([0.2 0.1 0.1/2 0.1/4 0.1/8],[0.01 0.005 0.0025 0.0025/2 0.0025/4],'-','LineWidth',3,'Color','r')
    hold on
    
    % All tests
    for iter3 = 1:size(L2_error,1)
        loglog(h(iter3,:),L2_error(iter3,:),'-.','LineWidth',1.5,'Color',[0    0.5686    0.7157])
        hold on
    end
    
    % Average error
    mean_vect = zeros(1,size(L2_error,2));
    mean_h = zeros(1,size(L2_error,2));
    for iter4 = 1:size(L2_error,2)
        mean_vect(iter4)=mean(L2_error(:,iter4));
        mean_h(iter4)=mean(h(:,iter4));
    end
    
    std_dev = std(L2_error);
    errorbar(mean_h,mean_vect,std_dev,'x-','LineWidth',2.5,'Color',[0    0.1882    0.9059])
    hold on
    
    xlabel('h','FontSize',16)
    xticks(flip(mean_h))
    yticks(flip(mean_vect))
    ylabel('e','FontSize',16,'Rotation',0)
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) = p(1)*0.98;
    set(yh,'position',p);
    lgd = legend('1st order reference','','','','','','','','','','Error','Average error');
    lgd.FontSize=(14);
    legend('Location','northwest')
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TOTAL SYSTEM CONVERGENCE TEST %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'TotalSystemTest')
    % Load data
    indata = load('L2errorDataSystemTest1');
    L2_error = indata.L2_error;
    h = indata.h;
    
    % Make figure
    figure

    % Reference line
    loglog([1 0.5 0.5/2 0.5/4 0.5/8 0.5/16 0.5/32],[0.4 0.2 0.1 0.05 0.025 0.025/2 0.024/4],'-','LineWidth',2.5,'Color','r')
    hold on
    
    % All tests
    for iter3 = 1:size(L2_error,1)
        loglog(h(iter3,:),L2_error(iter3,:),'-.','LineWidth',1.5,'Color',[0    0.5686    0.7157])
        hold on
    end
    
    % Average error
    mean_vect = zeros(1,size(L2_error,2));
    mean_h = zeros(1,size(L2_error,2));
    for iter4 = 1:size(L2_error,2)
        mean_vect(iter4)=mean(L2_error(:,iter4));
        mean_h(iter4)=mean(h(:,iter4));
    end
    
    std_dev = std(L2_error);
    errorbar(mean_h,mean_vect,std_dev,'x-','LineWidth',2.5,'Color',[0    0.1882    0.9059])
    hold on
    
    xlabel('h','FontSize',16)
    xticks(flip(mean_h));
    yticks(flip(mean_vect));
    ylabel('e','FontSize',16,'Rotation',0)
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) = p(1)*0.95;
    set(yh,'position',p);
    lgd = legend('1st order reference','','','','','','','','','','Error','Average error');
    lgd.FontSize=(14);
    legend('Location','northwest')
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ILLUSTRATION OF PEACEMAN CORRECTION %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'PeacemanCorrection')
    figure('Name','Peaceman pressure function')
    nx=180;
    ny=180;
    points = GridGeneration(nx,ny,D);
    map = gray(2);
    for i = 1:size(cells,1)
        coords=vertices(cells{i},:);
        pgon = polyshape(coords(:,1),coords(:,2));
        pg = plot(pgon);
        pg.FaceColor=(map(end,:));
        hold on
        for j = 1:size(points,1)
            in(j) = inpolygon(points(j,1),points(j,2),coords(:,1),coords(:,2));
        end
        p_in=points(in==1,:);
        distToVessel = sqrt((p_in(:,1)-cell_center(i,1)).^2+(p_in(:,2)-cell_center(i,2)).^2);
        p_pm=-darcy_source(i)*mu/(2*k*pi)*log(distToVessel/edges(Tn(i,3),4))+p_tn(i);
        p_pm_sorted = sort(p_pm);
        dp = p_pm_sorted(2)-p_pm_sorted(1);
        p_axis=min(p_pm_sorted):dp:max(p_pm_sorted);
        newmap = parula(length(p_axis));
        for iter = 1:size(p_in,1)
            p_here = p_pm(iter);
            ind = find(abs(p_axis-p_here)==min(abs(p_axis-p_here)));
            ind = ind(1);
            plot(p_in(iter,1),p_in(iter,2),'.','Color',newmap(ind,:))
            hold on
        end
        c=scircle1(cell_center(i,1),cell_center(i,2),0.2*sqrt(cell_area(i)));
        plot(c(:,1),c(:,2),'--','LineWidth',0.4,'Color',newmap(1,:));
        hold on
    end
    axis(D)
    
    figure('Name','Peaceman pressure at 0.2*deltaX')
    for i = 1:size(cells,1)
        coords=vertices(cells{i},:);
        pgon = polyshape(coords(:,1),coords(:,2));
        pg = plot(pgon);
        pg.FaceColor=(map(end,:));
        hold on
        for j = 1:size(points,1)
            in(j) = inpolygon(points(j,1),points(j,2),coords(:,1),coords(:,2));
        end
        p_in=points(in==1,:);
        distToVessel = sqrt((p_in(:,1)-cell_center(i,1)).^2+(p_in(:,2)-cell_center(i,2)).^2);
        p_pm=-darcy_source(i)*mu/(2*k*pi)*log(distToVessel/edges(Tn(i,3),4))+p_tn(i);
        p_pm2 = -darcy_source(i)*mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4))+p_tn(i);
        p_pm_sorted = sort(p_pm);
        dp = p_pm_sorted(2)-p_pm_sorted(1);
        p_axis=min(p_pm_sorted):dp:max(p_pm_sorted);
        newmap = parula(length(p_axis));
        ind2 = find(abs(p_axis-p_pm2)==min(abs(p_axis-p_pm2)));
        ind2 = ind2(1);
        for iter = 1:size(p_in,1)
            plot(p_in(iter,1),p_in(iter,2),'.','Color',newmap(ind2,:))
            hold on
        end
    end
    plot(cell_center(:,1),cell_center(:,2),'.','MarkerSize',7,'Color',newmap(end,:));
    hold on
    axis(D)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Tree and impact field on same figure %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'ImpactFieldDeterministic')
    indata = load('DeterministicIF');
        figure()
        IntensityMap(indata.cells,indata.vertices,indata.KT)
        hold on
        axis(indata.D)
        DrawTree(indata.Tree,150,'b',indata.D);
        DrawTree(indata.DETtree,150,[0.8500, 0.3250, 0.0980],indata.D);
        plot(indata.nodes(indata.MicroTermIndexes,1),indata.nodes(indata.MicroTermIndexes,2),'b.','MarkerSize',10)
        plot(indata.nodes(indata.MacroTermIndexes(indata.TN),1),indata.nodes(indata.MacroTermIndexes(indata.TN),2),'.','MarkerSize',30,'Color',[0.8500, 0.3250, 0.0980])
        hold on
        axis off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error plot for different K^D values %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'ErrorPlotDifferentKD')
    indata = load('CalculationsDetErrorPlot');
    figure()
    for j = 1:indata.iterations1
        plot(indata.delta_m(:,j),indata.error(:,j),'.-','MarkerSize',20)
        hold on
    end
    xlabel('Δm','FontSize',15)
    ylabel('RMS(p^{D}_{exact}-p^{D}_{coarse})','FontSize',14,'Rotation',90)
    set(gca,'YScale','log');
    ylim([1E-17 1E-6])
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1;
    set(yh,'position',p);
    xticks(indata.delta_m(:,1))
    lgd = legend('K^D = 1000','K^D = 100','K^D = 10','K^D = 1','K^D = 0.1','Location','northeast');
    lgd.FontSize = 10;
    grid on
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Analytically K^T vs. computationally calculated K^T value %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure, 'TrueVScodeKT')
    figure()
    indata = load('AnalyticVScompKT1');
    map = turbo(size(indata.kT,2));
    for iter = 1:size(indata.kT,2)
        plot(indata.r(indata.NotZeroVals(:,iter)==1,iter),indata.kT(indata.NotZeroVals(:,iter)==1,iter),'.','MarkerSize',20,'Color',map(iter,:));
        yline(indata.TrueKT(iter),'LineWidth',2.5,'Color',map(iter,:),'LineStyle','--')
        hold on
    end
    ylim([min(max(indata.kT))-0.1*min(max(indata.kT)) max(max(indata.kT))+0.1*max(max(indata.kT))])
    yyaxis left
    yticks(sort(indata.TrueKT));
    ylabel('K^{T}','FontSize',15,'Rotation',0)
    yyaxis right
    ylim([min(max(indata.kT))-0.1*min(max(indata.kT)) max(max(indata.kT))+0.1*max(max(indata.kT))])
    ylabel('Δm','FontSize',15,'Rotation',0)
    yticks(sort(indata.TrueKT));
    yticklabels({'8','7','6','5','4','3','2','1'}');
    lgd = legend('Computationally calculated value','Analytically calculated value');
    lgd.FontSize = 15;
    xlabel('R','FontSize',15)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% K^T varying with radius rate %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure, 'KTradiusrateFIG')
    indata = load('KTvsDeltaRData');
    map=turbo(size(indata.TrueKT,1));
    for i = 1:size(indata.TrueKT,2)
        for j = 1:size(indata.TrueKT,1)
        plot(indata.rr(j,i),indata.TrueKT(j,i),'.','MarkerSize',20,'Color',map(j,:))
        hold on
        end
    end
    xlabel('Δr','FontSize',15)
    ylabel('K^T','FontSize',15,'Rotation',0)
    lgd = legend('Δm = 1','Δm = 2','Δm = 3','Δm = 4','Δm = 5','Δm = 6','Δm = 7','Δm = 8','Location','southeast');
    set(gca,'Yscale','log')
    grid on
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*0.95;
    set(yh,'position',p);
end
%%%%%%%%%%%%%%%%%%
%%%% K^T VS R %%%%
%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'KTvsRFig')
    figure()
    indata = load('KTvsR');
    map = turbo(size(indata.kT,2));
for iter = 1:size(indata.kT,2)
    plot(indata.r(indata.NotZeroVals(:,iter)==1,iter),indata.kT(indata.NotZeroVals(:,iter)==1,iter),'.','MarkerSize',20,'Color',map(iter,:));
    hold on
end
set(gca,'YScale','log');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compare pressure solution %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'PressureAndDeviation')
indata = load('PressurePlotData');
h1 = figure('Name','Pressure in Darcy domain');
D = indata.D;
% Exact model
subplot(1,2,1)
[xq,yq]=meshgrid(D(1):0.1:D(2), D(3):0.1:D(4));
p_points = [indata.cell_center;indata.boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
p_values = [indata.p_darcyEx;indata.Bv_darcy*ones(size(indata.boundary_cells,1)+4,1)];
vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
surf(xq,yq,vq)
title('Exact model','FontSize',15)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
zlabel('P^{D}_h','FontSize',13,'Rotation',0)
zh = get(gca,'zlabel');
p = get(zh,'position');
p(1) = 0.5*p(1) ;
set(zh,'position',p)

% Deviation (%)
subplot(1,2,2)
p_points = [indata.cell_center;indata.boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
p_values = [100*(indata.p_darcyCoarse-indata.p_darcyEx)./indata.p_darcyEx;indata.Bv_darcy*zeros(size(indata.boundary_cells,1)+4,1)];
vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
surf(xq,yq,vq)
title('Deviation (%)','FontSize',15)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
colormap jet 
end


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


% figure('Name','Pressure plot')
% [xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
% p_points = [cell_center(ok_cells,1:2);boundary_cells(:,1:2)];
% p_values = [p_darcy(ok_cells);bv];
% pEx_values = [p_exact(cell_center(ok_cells,1),cell_center(ok_cells,2));bv];
% vq = griddata(p_points(:,1),p_points(:,2),p_values-pEx_values,xq,yq);
% surf(xq,yq,vq)