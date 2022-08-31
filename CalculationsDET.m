%clc;
%close all
clear;
TN = 1;
Ncells = 10;
iterations = 1;    % How many times to change radiusrate
iterations1 = 1;    % How many times to change K^D

% Deterministic tree data
RootNode = [2.5 0];
DT.Levels = 9;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.05;    % mm
DT.RadiusRate = 0.4;
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
    for iter = 1:iterations
    CutLevel = DT.Levels-5;
    trees = {'Deterministic','Random','Combinated','Half Deterministic'};
    flag.case = trees{1};
    Tree = ChooseTree(flag.case,RandomTree,DT,D,Ncells);
    nodes = Tree.nodes; edges = Tree.edges;
    DETnodes = nodes(1:2^(CutLevel-1),:);
    DETedges = edges(1:2^(CutLevel-1)-1,:);
    DETtree.nodes = DETnodes;
    DETtree.edges = DETedges;
    DETtree.RootNodeIdx = 1;
    % Find terminal nodes
    [TNinfo,TNlogic]=FindTerminals(Tree);
    MicroTermIndexes = find(TNlogic==1);
    MacroTermIndexes = find(nodes(:,3)==CutLevel);
    
    % Fix edge radiee
%     rel = edges(:,1)./edges(:,4);
%     fix = find(rel<20);
%     edges(fix,4)=edges(fix,1)/100;
    
    %%% Voronoi diagram %%%
    [cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');
    
    %%% Parameters %%%
    mu = 3E-6;                                        % viscosity (Ns/mm^2)
%     if iter1 == 1
%         k = 3E-3;                                     % permeability (mm^2)
%     elseif iter1 == 2
%         k = 3E-4;
%     elseif iter1 == 3
        k = 3E-5;
%     elseif iter1 == 4
%         k = 3E-6;
%     elseif iter1 == 5
%         k = 3E-7;
%     end
    K_N = (pi*edges(:,4).^4./(8*mu*edges(:,1)));      % Conductance
    K_D = k/mu                                       % Hydraulic conductivity [mm^4/Ns)
    f = @(x,y) 0;                                     % External sources or sinks
    
    %%% Set boundary values %%%
    BCs = {'Dirichlet','Neumann'};
    flag.case = BCs{1};
    Neu_network = 5;
    Dir_network = 1;
    Bv_darcy = -1;
    
    %%% Solve coupled system with stocastic tree %%%
    [Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells,bv] = TPFA(cells,vertices,f,K_D,1,Bv_darcy,edges(TNinfo(end,3),4));
    [p_darcyEx,q_network,p_network]=SolveSystemEx(Tree,TNinfo,TNlogic,Dir_network,Neu_network,mu,k,K_N,LHS,RHS,cell_area,flag.case);
    
    %% Find K_T %%%
    [q_network,p_network]=SolveNetwork(Tree,K_N,TNlogic,-1);
    [K_T,connections] = findKT(edges,cell_area,MicroTermIndexes,MacroTermIndexes,TNinfo,q_network,p_network,mu);

    % Solve system with exact K^T value
    [TNinfoDT,TNlogicDT]=FindTerminals(DETtree);
    [p_darcyCoarse,q_T,q_network,p_network]=SolveSystemUpS(DETtree,TNlogicDT,TNinfoDT,Dir_network,Neu_network,mu,k,K_N,K_T,connections,LHS,RHS,cell_area,flag.case);
    
    %%% Root mean square error %%%
    e = 0;
    e1 = 0;
    for i = 1:length(p_darcyEx)
        e = e + (p_darcyEx(i)-p_darcyCoarse(i))^2*cell_area(i);
    end
    error(iter,iter1) = sqrt(e/D_area);
    delta_m(iter,iter1)=iter;
    p_mean(iter,iter1)=mean(p_darcyEx);
    
     if iter == iterations
        KT = K_T(:,TN);
        save('DeterministicIF','cells','vertices','KT','D','Tree','DETtree','MicroTermIndexes','MacroTermIndexes','TN','nodes');
        Plots('ImpactFieldDeterministic');
    end
    DT.RadiusRate = DT.RadiusRate + 0.1;
    end
end

%%%% PLOTS %%%%
%save('CalculationsDetErrorPlot','iterations1','delta_m','error');
%Plots('ErrorPlotDifferentKD')

save ('PressurePlotData','D','cell_center','boundary_cells','p_darcyEx','Bv_darcy','p_darcyCoarse');
Plots('PressureAndDeviation');




