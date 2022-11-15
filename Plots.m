function[] = Plots(ThisFigure,file)
    axisfontsize = 16;
    set(gca,'FontSize',15)
    legendfontsize = 14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TPFA TEST VORONOI GRID %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure, 'TPFAConvergenceVoronoi')
    % Load data
    indata = load(file);
    L2_error = indata.L2_error;
    h = indata.h;

    % Make figure
    figure

    % Reference line
    loglog([0.2 0.1 0.1/2 0.1/4 0.1/8 0.1/16],[0.005 0.0025 0.0025/2 0.0025/4 0.0025/8 0.0025/16],'-','LineWidth',3,'Color','r')
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

    xlim([0.02 0.15])
    ylim([0.0005 0.004])
    
    xlabel('h [mm]','FontSize',axisfontsize)
    xticks([0.03 0.05 0.1 0.15 0.2])
    yticks([0.0005 0.001 0.002 0.003])
    ylabel('e [kPa\cdot mm]','FontSize',axisfontsize)
%     yh = get(gca,'ylabel');
%     p = get(yh,'position');
%     p(1) = p(1)*3.5;
%     set(yh,'position',p);
    lgd = legend('1st order reference','','','','','','','','','','Error','Average error');
    lgd.FontSize=legendfontsize;
    legend('Location','northwest')
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TOTAL SYSTEM CONVERGENCE TEST %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'TotalSystemTest')
    % Load data
    indata = load(file);
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
    
    ylim([0.009 0.3])
    xlim([0.03 0.3])
    xlabel('h [mm]','FontSize',16)
    xticks([0.05 0.1 0.2]);
    yticks([0.01 0.02 0.04 0.08 0.2]);
    ylabel('e [kPa \cdot mm]','FontSize',16)
%     yh = get(gca,'ylabel');
%     p = get(yh,'position');
%     p(1) = p(1)*0.95;
%     set(yh,'position',p);
    lgd = legend('1st order reference','','','','','','','','','','Error','Average error');
    lgd.FontSize=(14);
    legend('Location','northwest')
    grid on
end

% figure('Name','Pressure plot')
% [xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
% p_points = [cell_center(ok_cells,1:2);boundary_cells(:,1:2)];
% p_values = [p_darcy(ok_cells);bv];
% pEx_values = [p_exact(cell_center(ok_cells,1),cell_center(ok_cells,2));bv];
% vq = griddata(p_points(:,1),p_points(:,2),p_values-pEx_values,xq,yq);
% surf(xq,yq,vq)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DeterministicTree %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'DeterministicTree')
    indata = load('./Files/DeterministicTreeToDraw.mat');
    DrawTree(indata.Tree,200,'b',indata.D)
    plot(indata.Tree.nodes(1,1),indata.Tree.nodes(1,2),'.','MarkerSize',40,'Color',[0.8500,0.3250, 0.0980])
    axis off
end

%%%%%%%%%%%%%%%%%%
%%%% DLA tree %%%%
%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'DLATree')
    indata = load('./Files/DLATreeToDraw.mat');
    DrawTree(indata.Tree,150,'b',indata.D)
    plot(0.503,0.497,'.','MarkerSize',30,'Color',[0.8500,0.3250, 0.0980])
    axis off
end

%%%%%%%%%%%%%
%%%% RRT %%%%
%%%%%%%%%%%%%
if strcmp(ThisFigure,'RRT')
    indata = load('./Files/RRTToDraw.mat');
    DrawTree(indata.Tree,150,'b',indata.D)
    plot(indata.Tree.nodes(1,1),indata.Tree.nodes(1,2),'.','MarkerSize',40,'Color',[0.8500,0.3250, 0.0980])
    axis off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Tree and impact field on same figure %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'ImpactFieldAndDeterministicTree')
    indata = load(file);
        figure()
        IntensityMap(indata.cells,indata.vertices,indata.K_T(:,indata.TN))
        hold on
        axis(indata.D)
        DrawTree(indata.Tree,150,'b',indata.D);
        indata.VisibleTree
        DrawTree(indata.VisibleTree,150,[0.8500, 0.3250, 0.0980],indata.D);
        plot(indata.Tree.nodes(indata.MicroTermIndexes,1),indata.Tree.nodes(indata.MicroTermIndexes,2),'b.','MarkerSize',10)
        plot(indata.Tree.nodes(indata.MacroTermIndexes(indata.TN),1),indata.Tree.nodes(indata.MacroTermIndexes(indata.TN),2),'.','MarkerSize',30,'Color',[0.8500, 0.3250, 0.0980])
        hold on
        axis off
end

if strcmp(ThisFigure,'ImpactFieldAndUnstructuredTree')
    indata = load(file);
        figure()
        IntensityMap(indata.cells,indata.vertices,indata.K_T(:,indata.TN))
        hold on
        axis(indata.D)
        InvisibleTree.nodes = indata.Tree.nodes(2:end,:);
        InvisibleTree.edges = indata.Tree.edges(2:end,:);
        InvisibleTree.edges(:,2)=InvisibleTree.edges(:,2)-1;
        InvisibleTree.edges(:,3)=InvisibleTree.edges(:,3)-1;
        DrawTree(InvisibleTree,150,'b',indata.D);
        plot(indata.Tree.nodes(indata.MicroTermIndexes,1),indata.Tree.nodes(indata.MicroTermIndexes,2),'b.','MarkerSize',8)
        plot(indata.Tree.nodes(indata.MacroTermIndexes(indata.TN),1),indata.Tree.nodes(indata.MacroTermIndexes(indata.TN),2),'.','MarkerSize',30,'Color',[0.8500, 0.3250, 0.0980])
        hold on
        axis off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error plot for different values, deterministic %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'ErrorPlot_KD_deltam')
    indata = load(file);
    figure()
    for j = 1:indata.iterations1
        if j == 3
            plot(indata.delta_m(:,j),indata.error(:,j)/2,'.-','MarkerSize',25,'LineWidth',3.5)
            hold on
        else
            plot(indata.delta_m(:,j),indata.error(:,j)/2,'.-','MarkerSize',25,'LineWidth',1.5)
            hold on
        end
    end
    set(gca,'FontSize',15);
    xlabel('L_{M} - L_{M''}','FontSize',axisfontsize)
    ylabel('e_{rel}','FontSize',axisfontsize,'Rotation',90)
    set(gca,'YScale','log');
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1;
    set(yh,'position',p);
    xticks(indata.delta_m(:,1))
%     lgd = legend('K^D = 1000','K^D = 100','K^D = 10','K^D = 1','K^D = 0.1','Location','southeast');
%     lgd.FontSize = legendfontsize;
    grid on
end

if strcmp(ThisFigure,'ErrorPlot_alfaR_KD')
    indata = load(file);
    figure()
    for j = 1:indata.iterations1
        if j == 3
            plot(indata.delta_m(:,j),indata.error(:,j)/2,'.-','MarkerSize',25,'LineWidth',3.5)
            hold on
        else
            plot(indata.delta_m(:,j),indata.error(:,j)/2,'.-','MarkerSize',25,'LineWidth',1.5-0*j)
            hold on
        end
    end
    set(gca,'FontSize',15);
    xlabel('α_r','FontSize',axisfontsize)
    ylabel('e_{rel}','FontSize',axisfontsize,'Rotation',90)
    set(gca,'YScale','log');
    %ylim([1E-19 1E-1])
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1;
    set(yh,'position',p);
    xticks(indata.delta_m(:,1))
    lgd = legend('K^D = 1000','K^D = 100','K^D = 10','K^D = 1','K^D = 0.1','Location','southeast');
    lgd.FontSize = legendfontsize;
    grid on
end
if strcmp(ThisFigure,'ErrorPlot_M_deltaM')
    indata = load(file);
    figure()
    for j = 1:indata.iterations1
        plot(indata.delta_m(:,j),indata.error(:,j)/2,'.-','MarkerSize',27,'LineWidth',4-0.5*j)
        hold on
    end
    set(gca,'FontSize',15);
    xlabel('L_{M} - L_{M''}','FontSize',axisfontsize)
    ylabel('e_{rel}','FontSize',axisfontsize,'Rotation',90)
    set(gca,'YScale','log');
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1;
    set(yh,'position',p);
    xticks(indata.delta_m(:,1))
    grid on
end

if strcmp(ThisFigure,'ErrorPlot_alfaR_M')
    indata = load(file);
    figure()
    for j = 1:indata.iterations1
        plot(indata.delta_m(:,j),indata.error(:,j)/2,'.-','MarkerSize',25,'LineWidth',4-0.5*j)
        hold on
    end
    set(gca,'FontSize',15);
    xlabel('α_r','FontSize',axisfontsize)
    ylabel('e_{rel}','FontSize',axisfontsize,'Rotation',90)
    set(gca,'YScale','log');
    %ylim([1E-19 1E-1])
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1;
    set(yh,'position',p);
    xticks(indata.delta_m(:,1))
    lgd = legend('M=8','M=7','M=6','M=5','M=4','Location','southeast');
    lgd.FontSize = legendfontsize;
    grid on
end

if strcmp(ThisFigure,'ErrorPlot_realradii')
    indata = load(file);
    figure()
    for j = 1:indata.iterations1
        plot(indata.delta_m(:,j),indata.error(:,j)/2,'.-','MarkerSize',25,'LineWidth',4-0.5*j)
        hold on
    end
    set(gca,'FontSize',15);
    xlabel('M','FontSize',axisfontsize)
    ylabel('e_{rel}','FontSize',axisfontsize,'Rotation',90)
    set(gca,'YScale','log');
    %ylim([1E-19 1E-1])
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1;
    set(yh,'position',p);
    xticks(indata.delta_m(:,1))
    lgd = legend('K^D = 1000','K^D = 100','K^D = 10','K^D = 1','K^D = 0.1','Location','southeast');
    lgd.FontSize = legendfontsize;
    grid on
end

if strcmp(ThisFigure,'ErrorPlot_RRT')
    indata = load(file);
    figure()
    for j = 1:size(indata.error,2)
        plot(indata.x(:,j),indata.error(:,j),'.-','MarkerSize',25,'LineWidth',4-0.5*j)
        hold on
    end
    set(gca,'FontSize',15);
    xlabel('Terminal nodes','FontSize',axisfontsize)
    ylabel('e_{rel}','FontSize',axisfontsize,'Rotation',90)
    set(gca,'YScale','log');
    %ylim([1E-19 1E-1])
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1;
    set(yh,'position',p);
    lgd = legend('K^D = 1000','K^D = 100','K^D = 10','K^D = 1','K^D = 0.1','Location','southeast');
    lgd.FontSize = legendfontsize;
    grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Analytically K^T vs. computationally calculated K^T value %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure, 'TrueVScodeKT')
    figure()
    indata = load(file);
    map = turbo(size(indata.kT,2));
    for iter = 1:size(indata.kT,2)
        plot(indata.r(indata.NotZeroVals(:,iter)==1,iter),indata.kT(indata.NotZeroVals(:,iter)==1,iter),'.','MarkerSize',20,'Color',map(iter,:));
        yline(indata.TrueKT(iter),'LineWidth',2.5,'Color',map(iter,:),'LineStyle','--')
        hold on
    end
    ylim([min(max(indata.kT))-0.1*min(max(indata.kT)) max(max(indata.kT))+0.3*max(max(indata.kT))])
    yyaxis left
    %yticks(sort(indata.TrueKT));
    ylabel('K^{T} [(s\cdot kPa)^{-1}]','FontSize',axisfontsize)
    yyaxis right
    ylim([min(max(indata.kT))-0.1*min(max(indata.kT)) max(max(indata.kT))+0.3*max(max(indata.kT))])
    ylabel('L_{M}-L_{M''}','FontSize',axisfontsize)
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1.055;
    set(yh,'position',p);
    yticks(sort(indata.TrueKT));
    yticklabels({'8','','','','4','3','2','1'}');
    lgd = legend('Computationally calculated value','Analytically calculated value');
    lgd.FontSize = axisfontsize;
    xlabel('R [mm]','FontSize',axisfontsize)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% K^T varying with radius rate %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure, 'KTradiusrateFIG')
    figure
    indata = load(file);
    map=turbo(size(indata.TrueKT,1));
    for j = 1:size(indata.TrueKT,1)
    plot(indata.rr(j,:),indata.TrueKT(j,:),'.-','MarkerSize',20,'LineWidth',1.5,'Color',map(j,:))
    hold on
    end
    xlabel('α_r','FontSize',axisfontsize)
    ylabel('K^T [(s\cdot kPa)^{-1}]','FontSize',axisfontsize)
    lgd = legend('L_{M}-L_{M''} = 1','L_{M}-L_{M''} = 2','L_{M}-L_{M''} = 3','L_{M}-L_{M''} = 4','L_{M}-L_{M''} = 5','L_{M}-L_{M''} = 6','L_{M}-L_{M''} = 7','L_{M}-L_{M''} = 8','Location','southeast');
    lgd.FontSize = axisfontsize-2;
    set(gca,'Yscale','log')
    grid on
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*0.98;
    p(2) =p(2)*1.3;
    set(yh,'position',p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% KT vs. R on DLA tree %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'KTvsRDLAfig')
    figure()
    maxKT = 0;
    minKT = 10000;
    for i = 1:3
        if i == 3
            indata = load('./Files/KTvsR_DLA2000');
        elseif i == 2
            indata = load('./Files/KTvsR_DLA4000');
        elseif i == 1
            indata = load('./Files/KTvsR_DLA8000');
        end
        map = flip(turbo(7),1);
        % med map = flip(turbo(7)): map(6) = blå, map(4) = grønn, map(2) = rød
        ind_area = zeros(size(indata.kT,1),1);
        for iter = 1:size(indata.kT,1)
            ind_area(iter) = indata.cell_area(iter)/(max(indata.cell_area)-min(indata.cell_area));
        end
        (log2(ind_area)+0.1+abs(min(log2(ind_area))))*10
        scatter(indata.r,indata.kT,(log2(ind_area)+0.1+abs(min(log2(ind_area))))*15,'MarkerFaceColor',map(8-2*i,:),'MarkerEdgeColor',map(8-2*i,:));
        hold on
        if max(max(indata.kT)) > maxKT
            maxKT = max(max(indata.kT));
        end
        if min(min(indata.kT)) < minKT
            minKT = min(min(indata.kT));
        end
    end
    lgd = legend('8000 particles','4000 particles','2000 particles');
    lgd.FontSize= legendfontsize;
    set(gca,'YScale','log');
    ylim([minKT-0.1*minKT maxKT+4*maxKT])
    xlim([0 0.4])
    %xlim([0 max(max(indata.r))])
    ylabel('K^T [(s \cdot kPa)^{-1}]','FontSize',axisfontsize)
    xlabel('R [mm]','FontSize',axisfontsize)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Coefficient vs. particles DLA  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'aVSparticlesDLA')
map = flip(turbo(7),1);
    for i = 1:5
        if i == 1
            indata = load('./Files/KTvsR_DLA1000')
            data1 = indata;
        elseif i == 2
            indata = load('./Files/KTvsR_DLA2000');
            data2 = indata;
        elseif i == 3
            indata = load('./Files/KTvsR_DLA4000');
            data3 = indata;
        elseif i == 4
            indata = load('./Files/KTvsR_DLA8000');
            data4 = indata;
        elseif i == 5
            indata = load('./Files/KTvsR_DLA16000');
        end
        r_points = [];
        kt_points = [];
        for j = 1:size(indata.kT,2)
            r_points = [r_points;indata.r(:,j)];
            kt_points = [kt_points;indata.kT(:,j)];
        end
        linreg = fitlm(r_points,log10(kt_points))
        a_vect(i)=linreg.Coefficients.Estimate(2);
        b_vect(i)=linreg.Coefficients.Estimate(1);
        p_vect1(i)=linreg.Coefficients.pValue(2);
    end
    save('./Files/coefficients','a_vect','b_vect');
    % Plot slope vs. nodes
    rr = [1000 2000 4000 8000 16000];
    plot(rr,a_vect,'.-','MarkerSize',20,'LineWidth',2.5,'Color',[0 0 0])
    hold on
    plot(rr(1),a_vect(1),'.','MarkerSize',35,'Color',map(2,:));
    hold on
    plot(rr(2),a_vect(2),'.','MarkerSize',35,'Color',map(4,:));
    hold on
    plot(rr(3),a_vect(3),'.','MarkerSize',35,'Color',map(6,:));
    hold on
    plot(rr(4),a_vect(4),'.','MarkerSize',35,'Color',map(3,:));
    hold on 
    xlabel('Number of particles','FontSize',axisfontsize)
    ylabel('Slope coefficient','FontSize',axisfontsize)
    grid on
    set(gca,'FontSize',15)

    figure()

    % Plot data points and linear regression line

    ind_area4 = data4.cell_area./(max(data4.cell_area)-min(data4.cell_area));
    ind_area3 = data3.cell_area./(max(data3.cell_area)-min(data3.cell_area));
    ind_area2 = data2.cell_area./(max(data2.cell_area)-min(data2.cell_area));
    ind_area1 = data1.cell_area./(max(data1.cell_area)-min(data1.cell_area));

    map1 = turbo(300);
    scatter(data4.r,log10(data4.kT),(log2(ind_area4)+0.1+abs(min(log2(ind_area4))))*10,'MarkerFaceColor',map(3,:),'MarkerEdgeColor',map(3,:));
    %plot(data3.r(:,end), log10(data3.kT(:,end)),'.','MarkerSize',12,'Color',map(6,:));
    hold on

    X4 = [0:0.001:max(data4.r(:,end))];
    y4 = @(x) a_vect(4)*X4+b_vect(4);

    scatter(data3.r,log10(data3.kT),(log2(ind_area3)+0.1+abs(min(log2(ind_area3))))*10,'MarkerFaceColor',map(6,:),'MarkerEdgeColor',map(6,:));
    %plot(data3.r(:,end), log10(data3.kT(:,end)),'.','MarkerSize',12,'Color',map(6,:));
    hold on

    X3 = [0:0.001:max(data3.r(:,end))];
    y3 = @(x) a_vect(3)*X3+b_vect(3);

    scatter(data2.r,log10(data2.kT),(log2(ind_area2)+0.1+abs(min(log2(ind_area2))))*10,'MarkerFaceColor',map(4,:),'MarkerEdgeColor',map(4,:));
    %plot(data2.r(:,end),log10(data2.kT(:,end)),'.','MarkerSize',16,'Color',map(4,:));
    hold on
    
    X2 = [0:0.001:max(data2.r(:,end))];
    y2 = @(x) a_vect(2)*X2 + b_vect(2);

    scatter(data1.r,log10(data1.kT),(log2(ind_area1)+0.1+abs(min(log2(ind_area1))))*10,'MarkerFaceColor',map(2,:),'MarkerEdgeColor',map(2,:));
    %plot(data1.r(:,end),log10(data1.kT(:,end)),'.','MarkerSize',24,'Color',map(2,:));
    hold on
    
    X1 = [0:0.001:max(data1.r(:,end))];
    y1 = @(x) a_vect(1)*X1 +b_vect(1);

    plot(X1,y1(X1),'-','LineWidth',4,'Color',map1(254,:));
    hold on
    plot(X2,y2(X2),'-','LineWidth',4,'Color',map1(160,:));
    hold on
    plot(X3,y3(X3),'-','LineWidth',4,'Color',map1(25,:));
    hold on
    plot(X4,y4(X4),'-','LineWidth',4,'Color',map1(200,:));
    
    %set(gca,'Yscale','log')
    ylabel('log_{10}(K^T [(s \cdot kPA)^{-1}])','FontSize',axisfontsize)
    xlabel('R [mm]','FontSize',axisfontsize)
    lgd = legend('8000 particles','','','','','','','','','','4000 particles','','','','','','','','','',...
        '2000 particles','','','','','','','','','','1000 particles');
    %lgd = legend('','','','1600 particles','2000 particles','500 particles');
    lgd.FontSize = legendfontsize;
    set(gca,'FontSize',15)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Pressure plot of network and Darcy domain %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'PressurePlotNetworkDarcy')
    indata = load(file);
    figure('Name','Pressure plot');
    nodes = indata.Tree.nodes;
    edges = indata.Tree.edges;
    p_network = indata.p_networkEx;
    TNinfo = indata.TNinfo;
    D = indata.D;
    for i = 1:size(edges,1)
        x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
        y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
        pressure = [p_network(edges(i,2)) p_network(edges(i,3))];
        plot3(x,y,pressure,'.-','Color',[0 0 0],'LineWidth',edges(i,4)*100)
        hold on
        grid on
    end
    p_map = autumn(size(nodes,1)+size(TNinfo,1)*3);
    network_map = p_map(size(TNinfo,1)*2:end,:);
    for i = 1:size(nodes,1)
        x = nodes(i,1);
        y = nodes(i,2);
        pressure = p_network(i);
        p_sorted = sort(p_network);
        ind = find(p_sorted==pressure);
        ind = ind(1);
        plot3(x,y,pressure,'.','Color',network_map(ind,:),'MarkerSize',20);
        hold on
    end
    zlabel('Node pressure','FontSize',axisfontsize)
    xlabel('x','FontSize',axisfontsize)
    ylabel('y','FontSize',axisfontsize)

    darcy_map = p_map(1:3*size(TNinfo,1),:);
    
    [xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
    p_points = [TNinfo(:,1:2);indata.boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
    p_values = [indata.p_darcyEx;indata.Bv_darcy*ones(size(indata.boundary_cells,1)+4,1)];
    vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
    surf(xq,yq,vq)
    colormap(darcy_map)
    axis tight
end

if strcmp(ThisFigure,'PressurePlotNetworkDarcyCoarse')
    indata = load(file);
    figure('Name','Pressure plot');
    nodes = indata.DETtree.nodes;
    edges = indata.DETtree.edges
    p_network = indata.p_network;
    TNinfo = indata.TNinfo;
    D = indata.D;
    p_map = autumn(size(nodes,1)+size(TNinfo,1)*3);
    network_map = p_map(size(TNinfo,1)*2:end,:);
    for i = 2:size(nodes,1)
        x = nodes(i,1);
        y = nodes(i,2);
        pressure = p_network(i);
        p_sorted = sort(p_network);
        ind = find(p_sorted==pressure);
        ind = ind(1);
        plot3(x,y,pressure,'.','Color',p_map(end-3,:),'MarkerSize',20);
        hold on
    end
    zlabel('Node pressure')
    
    darcy_map = p_map(1:3*size(TNinfo,1),:);
    
    [xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
    p_points = [TNinfo(:,1:2);indata.boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
    p_values = [indata.p_darcyCoarse;indata.Bv_darcy*ones(size(indata.boundary_cells,1)+4,1)];
    vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
    surf(xq,yq,vq)
    colormap(darcy_map)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Illustration of exact model %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'ExactNetworkModel')
    indata = load(file);
    nodes = indata.Tree.nodes;
    edges = indata.Tree.edges;
    p_network = indata.p_networkEx;
    TNinfo = indata.TNinfo;
    D = indata.D;
    for i = 1:size(edges,1)
        x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
        y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
        pressure = [p_network(edges(i,2)) p_network(edges(i,3))];
        if nodes(edges(i,3),3) <= 3
            %plot3(x,y,pressure,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',5)
            plot3(x,y,pressure,'.-','Color',[0 0.4470 0.7410],'LineWidth',5)
        else
            %jenny = 1;
            plot3(x,y,pressure,'.-','Color',[0 0.4470 0.7410],'LineWidth',2)
        end
        if nodes(edges(i,3),3) == 3
            %plot3(nodes(edges(i,3),1),nodes(edges(i,3),2),pressure(2),'.','MarkerSize',30,'Color',[0.8500 0.3250 0.0980])
            jenny = 2;
        elseif nodes(edges(i,3),3)==5
            plot3(nodes(edges(i,3),1),nodes(edges(i,3),2),pressure(2),'.','MarkerSize',25,'Color',[0 0.4470 0.7410])
        end
        hold on
        grid on
    end
    
    p_map = gray(size(nodes,1)+size(TNinfo,1)*3);
    network_map = p_map(size(TNinfo,1)*2:end,:);
    
    darcy_map = gray(300); %p_map(1:3*size(TNinfo,1),:);
    
    [xq,yq]=meshgrid(D(1):0.005:D(2), D(3):0.005:D(4));
    p_points = [TNinfo(:,1:2);indata.boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
    p_values = [indata.p_darcyEx;indata.Bv_darcy*ones(size(indata.boundary_cells,1)+4,1)];
    vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
    surf(xq,yq,vq)
    colormap(darcy_map)
    axis off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error plot different values, DLA %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'ErrorPlotDifferentValuesDLA')
    data1 = load('./Files/DLAerror1');
    data2 = load('./Files/DLAerror2');
    data3 = load('./Files/DLAerror3');
    data4 = load('./Files/DLAerror4');
    data5 = load('./Files/DLAerror5');
    data6 = load('./Files/DLAerror6');
    data7 = load('./Files/DLAerror7');
    data8 = load('./Files/DLAerror8');
    data9 = load('./Files/DLAerror9');
    data10 = load('./Files/DLAerror10');

    x = data1.x;
    e1 = data1.e_rel;
    e2 = data2.e_rel;
    e3 = data3.e_rel;
    e4 = data4.e_rel;
    e5 = data5.e_rel;
    e6 = data6.e_rel;
    e7 = data7.e_rel;
    e8 = data8.e_rel;
    e9 = data9.e_rel;
    e10 = data10.e_rel;


    data_mean = zeros(size(x));
    data_std = zeros(size(x));
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            e_vect = [e1(i,j);e2(i,j);e3(i,j);e4(i,j);e5(i,j);e6(i,j);e7(i,j);e8(i,j);e9(i,j);e10(i,j)];
            data_mean(i,j)=mean(e_vect);
            data_std(i,j)=std(e_vect);
        end
    end
    
    rgb=zeros(size(x,2),3);
    rgb(1,:)=[0 0.4470 0.7410];
    rgb(2,:)=[0.8500 0.3250 0.0980];
    rgb(3,:)=[0.9290 0.6940 0.1250];
    rgb(4,:)=[0.4940 0.1840 0.5560];
    rgb(5,:)=[0.4660 0.6740 0.1880];

    for k = 1:size(x,2)
        plot(x(:,k),data_mean(:,k),'.-','LineWidth',2,'MarkerSize',15);
        hold on
        patch([x(:,k); flipud(x(:,k))],[data_mean(:,k)-data_std(:,k); flipud(data_mean(:,k)+data_std(:,k))],rgb(k,:),'FaceAlpha',0.3,'EdgeColor','none');
        hold on
    end
    
    set(gca,'FontSize',15);
    xlabel('Number of particles','FontSize',axisfontsize)
    ylabel('e_{rel}','FontSize',axisfontsize,'Rotation',90)
    set(gca,'YScale','log');
    yh = get(gca,'ylabel');
    p = get(yh,'position');
    p(1) =p(1)*1;
    set(yh,'position',p);
    ylim([1E-6 1E0]);
    %xticks(indata.x(:,1))
    lgd = legend('K^D = 1000','','K^D = 100','','K^D = 10','','K^D = 1','','K^D = 0.1','','Location','northeast');
    lgd.FontSize = legendfontsize;
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error plot different values, RRT %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'ErrorPlotDifferentValuesRRT')
    data1 = load('./Files/RRTerror1');
    data2 = load('./Files/RRTerror2');
    data3 = load('./Files/RRTerror3');
    data4 = load('./Files/RRTerror4');
    data5 = load('./Files/RRTerror5');
    data6 = load('./Files/RRTerror6');
    data7 = load('./Files/RRTerror7');
    data8 = load('./Files/RRTerror8');
    data9 = load('./Files/RRTerror9');
    data10 = load('./Files/RRTerror10');
    x = data1.x;
    e1 = data1.e_rel;
    e2 = data2.e_rel;
    e3 = data3.e_rel;
    e4 = data4.e_rel;
    e5 = data5.e_rel;
    e6 = data6.e_rel;
    e7 = data7.e_rel;
    e8 = data8.e_rel;
    e9 = data9.e_rel;
    e10 = data10.e_rel;
    data_mean = zeros(size(x));
    data_std = zeros(size(x));
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            data_mean(i,j)=mean([e1(i,j);e2(i,j);e3(i,j);e4(i,j);e5(i,j);e6(i,j);e7(i,j);e8(i,j);e9(i,j);e10(i,j)]);
            data_std(i,j)=std([e1(i,j);e2(i,j);e3(i,j);e4(i,j);e5(i,j);e6(i,j);e7(i,j);e8(i,j);e9(i,j);e10(i,j)]);
%             if data_mean(i,j)<data_std(i,j)
%                 data_std(i,j)=data_mean(i,j)-0.0000000001;
%             end

        end
    end
    rgb=zeros(size(x,2),3);
    rgb(1,:)=[0 0.4470 0.7410];
    rgb(2,:)=[0.8500 0.3250 0.0980];
    rgb(3,:)=[0.9290 0.6940 0.1250];
    rgb(4,:)=[0.4940 0.1840 0.5560];
    rgb(5,:)=[0.4660 0.6740 0.1880];
    
    for k = 1:size(x,2)
        plot(x(:,k),data_mean(:,k),'.-','LineWidth',2,'MarkerSize',15);
        hold on
        if k == 2
            patch([x(:,k); flipud(x(:,k))],[data_mean(:,k); flipud(data_mean(:,k)+data_std(:,k))],rgb(k,:),'FaceAlpha',0.3,'EdgeColor','none');
            hold on
        else
            patch([x(:,k); flipud(x(:,k))],[data_mean(:,k)-data_std(:,k); flipud(data_mean(:,k)+data_std(:,k))],rgb(k,:),'FaceAlpha',0.3,'EdgeColor','none');
            hold on
        end
    end
    
    set(gca,'FontSize',15);
    xlabel('Number of terminal nodes','FontSize',axisfontsize)
    ylabel('e_{rel}','FontSize',axisfontsize,'Rotation',90)
    ylim([1E-6 1E0])
    %ylim([-0.001 0.12]);
    xlim([100 1600])
    set(gca,'YScale','log');
    %xticks(indata.x(:,1))
    lgd = legend('K^D = 1000','','K^D = 100','','K^D = 10','','K^D = 1','','K^D = 0.1','','Location','northwest');
    lgd.FontSize = legendfontsize;
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Coeff vs. nodes, RRT %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ThisFigure,'aVSparticlesRRT')
map = flip(turbo(7),1);
figure()
    for i = 1:9
        if i == 1
            indata = load('./Files/KTvsR_RRT100');
            data1 = indata;
        elseif i == 2 
            indata = load('./Files/KTvsR_RRT200');
        elseif i == 3
            indata = load('./Files/KTvsR_RRT400');
            data2 = indata;
        elseif i == 4
            indata = load('./Files/KTvsR_RRT600');
        elseif i == 5
            indata = load('./Files/KTvsR_RRT800');
        elseif i == 6 
            indata = load('./Files/KTvsR_RRT1000');
        elseif i == 7
            indata = load('./Files/KTvsR_RRT1200');
        elseif i == 8 
            indata = load('./Files/KTvsR_RRT1400');
        elseif i == 9
            indata = load('./Files/KTvsR_RRT1600');
            data3 = indata;
        end
        r_points = [];
        kt_points = [];
        for j = 1:size(indata.kT,2)
            r_points = [r_points;indata.r(:,j)];
            kt_points = [kt_points;indata.kT(:,j)];
        end
        linreg = fitlm(r_points,log10(kt_points))
        a_vect(i)=linreg.Coefficients.Estimate(2);
        b_vect(i)=linreg.Coefficients.Estimate(1);
        p_vect1(i)=linreg.Coefficients.pValue(2);
    end
    map1 = turbo(300);
    save('./Files/Coefficients_RRT','a_vect','b_vect');
    % Plot slope vs. nodes
    rr = [100 200 400 600 800 1000 1200 1400 1600];
    plot(rr,a_vect,'.-','MarkerSize',20,'LineWidth',2.5,'Color',map(3,:))
    hold on
    plot(rr(1),a_vect(1),'.','MarkerSize',30,'Color',map(2,:));
    hold on
    plot(rr(3),a_vect(3),'.','MarkerSize',30,'Color',map(4,:));
    hold on
    plot(rr(9),a_vect(9),'.','MarkerSize',30,'Color',map(6,:));
    hold on 
    xlabel('Terminal nodes','FontSize',axisfontsize)
    ylabel('Slope coefficient','FontSize',axisfontsize)
    grid on
    set(gca,'FontSize',15)

    figure()

    % Plot data points and linear regression line
    plot(data3.r(:,end), log10(data3.kT(:,end)),'.','MarkerSize',12,'Color',map(6,:));
    hold on

    X = 0:0.001:max(data3.r(:,end));

    y3 = @(x) a_vect(9)*X+b_vect(9);

    plot(data2.r(:,end),log10(data2.kT(:,end)),'.','MarkerSize',16,'Color',map(4,:));
    hold on

    y2 = @(x) a_vect(3)*X + b_vect(3);

    plot(data1.r(:,end),log10(data1.kT(:,end)),'.','MarkerSize',24,'Color',map(2,:));
    hold on

    y1 = @(x) a_vect(1)*X +b_vect(1);
% 
    plot(X,y1(X),'-','LineWidth',4,'Color',map1(254,:));
    hold on
    plot(X,y2(X),'-','LineWidth',4,'Color',map1(160,:));
    hold on
    plot(X,y3(X),'-','LineWidth',4,'Color',map1(25,:));
    hold on
    
    %set(gca,'Yscale','log')
    xlim([0 0.8])
    ylim([-8 8])
    ylabel('log_{10}(K^T [(s \cdot kPA)^{-1}])','FontSize',axisfontsize)
    xlabel('R [mm]','FontSize',axisfontsize)
    lgd = legend('1600 terminal nodes','400 terminal nodes','1400 terminal nodes','','','');
    lgd.FontSize = legendfontsize;
    set(gca,'FontSize',15)
end
set(gca,'FontSize',15)
end