%clc;
close all
clear;
TN = 1;

% Domain
D = [0 5 0 5/sqrt(2)];
D_area = (D(2)-D(1))*(D(4)-D(3));

% How many times to change delta m
iterations = 8;

% How many times to change K^D
iterations1 = 5;    

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

% RRT data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 0.7;
RandomTree.Ncells = 1600;
RandomTree.TerminalRadius = 0.004;

% DLA tree data
DLA.Nparticles = 100;
DLA.RootRadius = 0.5;
DLA.TerminalRadius = 0.004;

% Make tree
for iter1 = 1:iterations1
    % For testing
    %DT.RadiusRate = 0.4;
    DT.Levels = 4;
    M = DT.Levels-1;
    for iter = 1:iterations
    M = DT.Levels-1;
    CutLevel = 1;
    trees = {'Deterministic','Random','Combinated','Half Deterministic'};
    ChosenTree = trees{1};
    Tree = ChooseTree(ChosenTree,RandomTree,DT,DLA,D);
    nodes = Tree.nodes; edges = Tree.edges;
    DETnodes = nodes(1:2^(CutLevel),:);
    DETedges = edges(1:2^(CutLevel)-1,:);
    DETtree.nodes = DETnodes;
    DETtree.edges = DETedges;
    DETtree.RootNodeIdx = 1;

    % Find terminal nodes
    [TNinfo,TNlogic]=FindTerminals(Tree);
    MicroTermIndexes = find(TNlogic==1);
    MacroTermIndexes = find(nodes(:,3)==CutLevel);
    
    if strcmp(ChosenTree,'Deterministic')
        RootRadius = DT.TrunkRadius;
        TerminalRadius = RandomTree.TerminalRadius;
        MicroTermIndexes = find(TNlogic==1);
        [Graph, ~, ~]=makegraph(Tree.edges(:,2:3),size(Tree.edges,1),size(Tree.nodes,1));
        NodeLevels = zeros(size(Tree.nodes,1),1);
        for i = 1:length(MicroTermIndexes)
            SP = shortestpath(Graph,MicroTermIndexes(i),1);
            for j = 1:size(Tree.nodes,1)
                if any(ismember(SP,j)) && NodeLevels(j) < find(ismember(SP,j))-1
                    NodeLevels(j) = find(ismember(SP,j))-1;
                end
            end
        end
        
        
        power = log(RootRadius/TerminalRadius)/(max(Tree.nodes(:,3))-1);
        RadiusRate = exp(power);
        for i = 2:length(NodeLevels)
            edge = find(Tree.edges(:,3)==i);
            edgerad = TerminalRadius*RadiusRate^(NodeLevels(i));
            Tree.edges(edge,4)=edgerad;
        end

        edges = Tree.edges;
        nodes = Tree.nodes;
     
    end 
    
    %%% Voronoi diagram %%%
    [cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');
    
    %%% Parameters %%%

    % viscosity (Ns/mm^2)
    mu = 3E-6;   

    % permeability (mm^2)
    %k = 3E-5;
    if iter1 == 1
        k = 3E-3;                                     
    elseif iter1 == 2
        k = 3E-4;
    elseif iter1 == 3
        k = 3E-5;
    elseif iter1 == 4
        k = 3E-6;
    elseif iter1 == 5
        k = 3E-7;
    end

    % Edge conductance
    K_N = (pi*edges(:,4).^4./(8*mu*edges(:,1))); 

    % Hydraulic conductivity [mm^4/Ns)
    K_D = k/mu; 

    % External sources or sinks
    f = @(x,y) 0;                                     
    
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

    % Change when testing (x-axis in RMS-plot)
    delta_m(iter,iter1)=M;
    p_mean(iter,iter1)=mean(p_darcyEx);
    
     if iter == iterations
        KT = K_T(:,TN);
        save('DeterministicIF','cells','vertices','KT','D','Tree','DETtree','MicroTermIndexes','MacroTermIndexes','TN','nodes');
     end
    % For testing
     %DT.RadiusRate = DT.RadiusRate + 0.1;
     DT.Levels = DT.Levels+1;
    end
    % For testing
    disp(iter1)
end

%%%% PLOTS %%%%
%save('CalculationsDetErrorPlot','iterations1','delta_m','error');
%Plots('ErrorPlotDifferentKD')
%Plots('ImpactFieldDeterministic');

%save ('PressurePlotData','D','cell_center','boundary_cells','p_darcyEx','Bv_darcy','p_darcyCoarse');
%Plots('PressureAndDeviation');

%save('VaryDeltam_and_M','iterations1','delta_m','error')
%save('VaryDeltam_and_KD','iterations1','delta_m','error')
%save('VaryAlphaR_and_KD','iterations1','delta_m','error')
save('RealRadiiSetup','iterations1','delta_m','error')
Plots('ErrorPlotDifferentValues','RealRadiiSetup')


