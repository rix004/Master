%clc;
close all
clear;
TN = 1;
iterations = 10;
iterations1 = 10;

% Deterministic tree data
RootNode = [2.5 0];
DT.Levels = 4;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.05;    % mm
DT.RadiusRate = 0.7;
DT.TrunkLength = 5/sqrt(2)/2; %mm
DT.LengthRate = 1/sqrt(2);

% Random tree data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 0.7;

% Domain
D = [0 5 0 5/sqrt(2)];
D_area = (D(2)-D(1))*(D(4)-D(3));

% Make tree
for iter1 = 1:iterations1
    Ncells = 100;
    for iter = 1:iterations
    CutLevel = DT.Levels;
    trees = {'Deterministic','Random','Combinated','Half Deterministic'};
    flag.case = trees{3};
    Tree = ChooseTree(flag.case,RandomTree,DT,D,Ncells);
    nodes = Tree.nodes; edges = Tree.edges;
    DETnodes = nodes(1:2^(CutLevel-1),:);
    DETedges = edges(1:2^(CutLevel-1)-1,:);
    
    % Find terminal nodes
    [TNinfo,TNlogic]=FindTerminals(nodes,edges);
    MicroTermIndexes = find(TNlogic==1);
    MacroTermIndexes = find(nodes(:,3)==CutLevel);
    
    % Fix edge radiee
    rel = edges(:,1)./edges(:,4);
    fix = find(rel<20);
    edges(fix,4)=edges(fix,1)/100;
    
    %%% Voronoi diagram %%%
    [cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');
    
    %%% Parameters %%%
    mu = 3E-6;                                        % viscosity (Ns/mm^2)
%     if iter1 == 1
%         k = 3E-12;                                     % permeability (mm^2)
%     elseif iter1 == 2
%         k = 3E-8;
%     elseif iter1 == 3
%         k = 3E-10;
%     elseif iter1 == 4
%         k = 3E-12;
%     end
    k = 3E-6;
    K_N = (pi*edges(:,4).^4./(8*mu*edges(:,1)));      % Conductance
    K_D = k/mu;                                       % Hydraulic conductivity [mm^4/Ns)
    f = @(x,y) 0;                                     % External sources or sinks
    
    %%% Set boundary values %%%
    BCs = {'Dirichlet','Neumann'};
    flag.case = BCs{1};
    Neu_network = 5;
    Dir_network = 1;
    Bv_darcy = -1;
    
    %%% Solve coupled system with stocastic tree %%%
    [Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells,bv] = TPFA(cells,vertices,f,K_D,1,Bv_darcy,edges(TNinfo(end,3),4));
    [p_darcyEx,q_T,q_network,p_network]=SolveSystemEx(nodes,edges,TNinfo,TNlogic,Dir_network,Neu_network,mu,k,K_N,LHS,RHS,cell_area,flag.case);
    
    % Test system
    %SystemTestEx(nodes,edges,TNlogic,TNinfo,mu,k,p_darcy,q_darcy,p_network,q_network,Grad_D,boundary_cells,cell_area,D_bvs,K_N,MacroTermIndexes,MicroTermIndexes)
    
    %% Find K_T %%%
    [q_network,p_network]=SolveNetwork(nodes,edges,K_N,TNlogic,-1);
    [K_T,connections] = findKT1(edges,cell_area,MicroTermIndexes,MacroTermIndexes,TNinfo,q_network,p_network,mu);

    % Solve system with exact K^T value
    [TNinfoDT,TNlogicDT]=FindTerminals(DETnodes,DETedges);
    [p_darcyCoarse,q_T,q_network,p_network]=SolveSystemUpS1(DETnodes,DETedges,TNlogicDT,TNinfoDT,Dir_network,Neu_network,mu,k,K_N,K_T,connections,LHS,RHS,cell_area,flag.case);

    %%% Use linear regression %%%
    newK_T = zeros(size(K_T));
    for iter2 = 1:length(MacroTermIndexes)
        r = sqrt((TNinfo(:,1)-nodes(MacroTermIndexes(iter2),1)).^2+(TNinfo(:,2)-nodes(MacroTermIndexes(iter2),2)).^2);
        kT=K_T(:,iter2);
        [a,b] = LinReg(kT,r,1,MicroTermIndexes,MacroTermIndexes,0);
        a_val(iter2)=a;
        con = connections';
        newK_T(find(con(:,iter2)==1),iter2) = exp(a*r(find(con(:,iter2)))+b);
    end
    aValue(iter,iter1)=mean(a_val);
    %StdDev(iter,iter1)=std(a_val);
    % Solve system with kT function
    [TNinfoDT,TNlogicDT]=FindTerminals(DETnodes,DETedges);
    [p_darcyCoarse1,q_T,q_network,p_network]=SolveSystemUpS1(DETnodes,DETedges,TNlogicDT,TNinfoDT,Dir_network,Neu_network,mu,k,K_N,newK_T,connections,LHS,RHS,cell_area,flag.case);
    
    %Test system
    %SystemTestCoarse(DETnodes,DETedges,TNlogicDT,TNinfoDT,mu,k,p_darcy1,q_T,p_network,q_network,Grad_D,boundary_cells,cell_area,D_bvs,K_N,MacroTermIndexes,connections)
    
    %%% Root mean square error %%%
%     e = 0;
%     e1 = 0;
%     for i = 1:length(p_darcyEx)
%         e = e + (p_darcyEx(i)-p_darcyCoarse(i))^2*cell_area(i);
%         e1 = e1 + (p_darcyEx(i)-p_darcyCoarse1(i))^2*cell_area(i);
%     end
%     error(iter,iter1) = sqrt(e/D_area);
%     error1(iter,iter1) = sqrt(e1/D_area);
%     l(iter,iter1)=Ncells;
    
%     if iter == iterations
%         figure(1)
%     %     IntensityMap(cells,vertices,kT(:,iter))
%     %     hold on
%     %     axis(D)
%         DrawTree(Tree,150,'b',D);
%         DETtree.nodes = DETnodes;
%         DETtree.edges = DETedges;
%         DrawTree(DETtree,150,[0.8500, 0.3250, 0.0980],D);
%         plot(nodes(MicroTermIndexes,1),nodes(MicroTermIndexes,2),'b.','MarkerSize',10)
%         plot(nodes(MacroTermIndexes(TN),1),nodes(MacroTermIndexes(TN),2),'.','MarkerSize',30,'Color',[0.8500, 0.3250, 0.0980])
%         hold on
%         axis off
%     end
    l(iter,iter1)=Ncells;
    Ncells = Ncells + 100;
    end
    disp(iter1);
end

% %%%% PLOTS %%%%
% figure(4)
% subplot(1,2,1);
% for j = 1:iterations1
%     plot(l(:,j),error(:,j),'.','MarkerSize',20)
%     hold on
% end
% title('Exact K^T value')
% xlabel('Number of cells','FontSize',15)
% ylabel('RMS','FontSize',15,'Rotation',0)
% set(gca,'YScale','log');
% ylim([1E-9 1E2])
% %lgd = legend('K^D = 1','K^D = 1E-2','K^D = 1E-4','K^D = 1E-6','Location','northeast');
% %lgd.FontSize = 12;
% yh = get(gca,'ylabel');
% p = get(yh,'position');
% p(1) = -1;
% set(yh,'position',p);
% xticks(l(:,1))
% subplot(1,2,2);
% for j = 1:iterations1
%     plot(l(:,j),error1(:,j),'.','MarkerSize',20)
% end
% title('K^T value with linear regression')
% xlabel('Number of cells','FontSize',15)
% ylabel('RMS','FontSize',15,'Rotation',0)
% set(gca,'YScale','log');
% ylim([1E-9 1E2])
% %lgd = legend('K^D = 1','K^D = 1E-2','K^D = 1E-4','K^D = 1E-6','Location','northeast');
% %lgd.FontSize = 12;
% yh = get(gca,'ylabel');
% p = get(yh,'position');
% p(1) = -1;
% set(yh,'position',p);
% xticks(l(:,1))
% 
% % Pressure plot
% h1 = figure('Name','Pressure in Darcy domain');
% subplot(1,3,1)
% [xq,yq]=meshgrid(D(1):0.1:D(2), D(3):0.1:D(4));
% p_points = [cell_center;boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
% p_values = [p_darcyEx;Bv_darcy*ones(size(boundary_cells,1)+4,1)];
% vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
% surf(xq,yq,vq)
% title('Exact model','FontSize',15)
% xlabel('x','FontSize',15)
% ylabel('y','FontSize',15)
% zlabel('P^{D}_h','FontSize',13,'Rotation',0)
% zh = get(gca,'zlabel'); % handle to the label object
% p = get(zh,'position'); % get the current position property
% p(1) = 0.5*p(1) ;        % double the distance, 
%                        % negative values put the label below the axis
% set(zh,'position',p)   % set the new position
% 
% h2 = figure();
% subplot(1,3,1)
% imagesc(vq)
% colorbar
% 
% figure(h1)
% subplot(1,3,2)
% p_points = [cell_center;boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
% p_values = [100*(p_darcyCoarse-p_darcyEx)./p_darcyEx;Bv_darcy*zeros(size(boundary_cells,1)+4,1)];
% vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
% surf(xq,yq,vq)
% title('Deviation (%)','FontSize',15)
% xlabel('x','FontSize',15)
% ylabel('y','FontSize',15)
% colormap jet 
% 
% figure(h2)
% subplot(1,3,2)
% imagesc(vq)
% colorbar
% 
% figure(h1)
% subplot(1,3,3)
% p_points = [cell_center;boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
% p_values = [100*(p_darcyCoarse1-p_darcyEx)./p_darcyEx;Bv_darcy*zeros(size(boundary_cells,1)+4,1)];
% vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
% surf(xq,yq,vq)
% title('Deviation (%) (Linear regression)','FontSize',15)
% xlabel('x','FontSize',15)
% ylabel('y','FontSize',15)
% colormap jet 
% 
% figure(h2)
% subplot(1,3,3)
% imagesc(vq)
% colorbar

figure(5)
errorbar(l(:,1),mean(aValue,2),std(aValue,0,2),'.-','LineWidth',2.5)
xlabel('Number of cells','FontSize',15)
ylabel('a','FontSize',15,'Rotation',0)

figure(6)
for h = 1:iter1
plot(l(:,h),aValue(h,:),'.','MarkerSize',10)
end
xlabel('Number of cells','FontSize',15)
ylabel('a','FontSize',15,'Rotation',0)
