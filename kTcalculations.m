%clc;
%close all
clear;
iterations = 1;
iterations1 = 1;
TN = 1;

% Domain
%D = [0 5 0 5/sqrt(2)];
D = [0 1 0 1];
D_area = (D(2)-D(1))*(D(4)-D(3));

% Deterministic tree data
RootNode = [0.5*(D(2)-D(1)) 0];
DT.Levels = 2;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.05;    % mm
DT.RadiusRate = 0.7;
DT.TrunkLength = 1/sqrt(2)/2; %mm
DT.LengthRate = 1/sqrt(2);
DT.RootNodeIdx = 1;

% RRT data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 1;
RandomTree.Ncells = 100;
RandomTree.TerminalRadius = 0.004;

% DLA tree data
DLA.Nparticles = 16000;
DLA.RootRadius = 0.05;
DLA.TerminalRadius = 0.004;
DLA.version = 1;

for iter1 = 1:iterations1
    for iter = 1:iterations
    % Where to cut
    M = DT.Levels -1;
    CutLevel = 1;
    trees = {'Deterministic','Random','Combinated','DLA'};
    ChosenTree = trees{3};
    Tree = ChooseTree(ChosenTree,RandomTree,DT,DLA,D);
    nodes = Tree.nodes; edges = Tree.edges;
    if strcmp(ChosenTree,'Deterministic') || strcmp(ChosenTree,'Combinated')
        % Fix edge radii
        VisibleNodes = nodes(1:2^(CutLevel),:);
        VisibleEdges = edges(1:2^(CutLevel)-1,:);
        VisibleTree.nodes = VisibleNodes;
        VisibleTree.edges = VisibleEdges;
        VisibleTree.RootNodeIdx = 1;
        edges = Tree.edges;
    elseif ChosenTree == 'DLA'
        VisibleTree.nodes = nodes(1:2,:);
        VisibleTree.edges = edges(1,:);
        InvisibleTree.nodes = nodes(2:end,:);
        InvisibleEdges = edges(2:end,:);
        InvisibleEdges(:,2)=InvisibleEdges(:,2)-1;
        InvisibleEdges(:,3)=InvisibleEdges(:,3)-1;
        InvisibleTree.edges = InvisibleEdges;
    end

    % Find terminal nodes
    [TNinfo,TNlogic]=FindTerminals(Tree);
    MicroTermIndexes = find(TNlogic==1);
    MacroTermIndexes = find(nodes(:,3)==CutLevel);

    % Fix edge radii
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
        1/RadiusRate
        for i = 2:length(NodeLevels)
            edge = find(Tree.edges(:,3)==i);
            edgerad = TerminalRadius*RadiusRate^(NodeLevels(i));
            Tree.edges(edge,4)=edgerad;
        end

        edges = Tree.edges;
        nodes = Tree.nodes;
     
    end 

    %DrawTree(Tree,10,[0.8500, 0.3250, 0.0980],D);

    %%% Voronoi diagram %%%
    [cells, vertices,orig_cells] = VoronoiDiagram(TNinfo(:,1:2),[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');
%     
    %%% Parameters %%%
    mu = 3E-6;                                        % viscosity (Ns/mm^2)
    k = 3E-5;                                         % permeability (mm^2)
    K_N = (pi*edges(:,4).^4./(8*mu*edges(:,1)));      % Conductance
    K_D = k/mu;                                       % Hydraulic conductivity [mm^4/Ns)
    f = @(x,y) 0;                                     % External sources or sinks
    
    %%% Set boundary values %%%
    BCs = {'Dirichlet','Neumann'};
    flag.case = BCs{1};
    Neu_network = 5;
    Dir_network = 1;
    Bv_darcy = -1;

    %%% Solve for pressure and flux in system %%%
    [Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells,bv] = TPFA(cells,vertices,f,K_D,1,Bv_darcy,edges(TNinfo(end,3),4));
    [q_network,p_network]=SolveNetwork(Tree,K_N,TNlogic,-1);

    %save('PressurePlotDataDLANetworkOnly','p_networkEx','Tree','TNinfo','D');
    %Plots('PressurePlotNetworkDarcy')

    %%% For validation %%%%
    if strcmp(ChosenTree,'Deterministic')
        summ = 0;
        for j = CutLevel:M-1
            summ = summ + DT.LengthRate^(j)/(2^(j)*DT.RadiusRate^(4*j));
        end
        TrueKT(iter,iter1) = pi*DT.TrunkRadius^4/(8*mu*DT.TrunkLength*D_area*summ);
    else
        TrueKT = 0;
    end

    %%% Find K^T %%%
    [K_T,connections] = findKT(edges,cell_area,MicroTermIndexes,MacroTermIndexes,TNinfo,q_network,p_network,mu);
    connections=connections';
    kT(:,iter)=K_T(:,TN);
    %kTnotZero = kT(connections(:,TN)==1,iter);
    kT_DET(iter)=K_T(1);
    Numb_levels(iter) = M;

    r(:,iter) = sqrt((TNinfo(:,1)-nodes(MacroTermIndexes(TN),1)).^2+(TNinfo(:,2)-nodes(MacroTermIndexes(TN),2)).^2);
    rVSk = [r(:,iter) kT(:,iter)];
    rVSkSorted = sortrows(rVSk,1);
    r(:,iter)=rVSkSorted(:,1);
    kT(:,iter)=rVSkSorted(:,2);
%     disp(iter)
%     rr(iter,iter1)=DT.RadiusRate;
%     levelcut(iter,iter1)=iter;
    %DT.Levels = DT.Levels-1;
    DLA.Version = DLA.version+1;
    end
    %DT.RadiusRate = DT.RadiusRate + 0.1;

end
NotZeroVals = zeros(size(kT));
NotZeroVals(kT>0)=1;

if strcmp(ChosenTree,'Deterministic')
    %save('TestDETtree','cells','vertices','K_T','Tree','VisibleTree','MicroTermIndexes','MacroTermIndexes','D','TN')
    % save('KTvsDeltaRData','rr','TrueKT');
    % save('AnalyticVScompKT1','NotZeroVals','kT','r','TrueKT')
    save('KT_Test_TrueRadii','kT_DET','Numb_levels');
    %Plots('KTvsRDET_trueradii');
elseif strcmp(ChosenTree,'DLA')
    %%% Impact field plot %%%
save('TestDLAtree','cells','vertices','K_T','Tree','VisibleTree','MicroTermIndexes','MacroTermIndexes','D','TN')
    if DLA.Nparticles == 2000
        save('KTvsR_DLA2000','kT','NotZeroVals','r','cell_area');
    elseif DLA.Nparticles == 4000
        save('KTvsR_DLA4000','kT','NotZeroVals','r','cell_area');
    elseif DLA.Nparticles == 8000
        save('KTvsR_DLA8000','kT','NotZeroVals','r','cell_area');
    elseif DLA.Nparticles == 100
        save('KTvsR_DLA100','kT','NotZeroVals','r','cell_area');
    elseif DLA.Nparticles == 200
        save('KTvsR_DLA200','kT','NotZeroVals','r','cell_area');
    elseif DLA.Nparticles == 500
        save('KTvsR_DLA500','kT','NotZeroVals','r','cell_area');
    elseif DLA.Nparticles == 1000
        save('KTvsR_DLA1000','kT','NotZeroVals','r','cell_area');
    elseif DLA.Nparticles == 250
        save('KTvsR_DLA250','kT','NotZeroVals','r','cell_area');
    elseif DLA.Nparticles == 16000
        save('KTvsR_DLA16000','kT','NotZeroVals','r','cell_area');
    end
elseif strcmp(ChosenTree,'Combinated')
    save('TestRRTtree','cells','vertices','K_T','Tree','VisibleTree','MicroTermIndexes','MacroTermIndexes','D','TN')
    Plots('ImpactFieldAndUnstructuredTree','TestRRTtree')
    if RandomTree.Ncells == 100
        save('KTvsR_RRT100','kT','NotZeroVals','r');
    elseif RandomTree.Ncells == 200
        save('KTvsR_RRT200','kT','NotZeroVals','r');
    elseif RandomTree.Ncells == 400
        save('KTvsR_RRT400','kT','NotZeroVals','r');
    elseif RandomTree.Ncells == 600
        save('KTvsR_RRT600','kT','NotZeroVals','r');
    elseif RandomTree.Ncells == 800
        save('KTvsR_RRT800','kT','NotZeroVals','r');
    elseif RandomTree.Ncells == 1000
        save('KTvsR_RRT1000','kT','NotZeroVals','r');
    elseif RandomTree.Ncells == 1200
        save('KTvsR_RRT1200','kT','NotZeroVals','r');
    elseif RandomTree.Ncells == 1400
        save('KTvsR_RRT1400','kT','NotZeroVals','r');
    elseif RandomTree.Ncells == 1600
        save('KTvsR_RRT1600','kT','NotZeroVals','r');
    end

    % Mean and std.dev. plot
    if RandomTree.Ncells == 16
        save('RRT_16cells','K_T');
    elseif RandomTree.Ncells == 64
        save('RRT_64cells','K_T');
    elseif RandomTree.Ncells == 256
        save('RRT_256cells','K_T')
    elseif RandomTree.Ncells == 1024
        save('RRT_1024cells','K_T')
    end
end




