%clc;
close all
clear;
iterations = 5;
iterations1 = 5;

% Deterministic tree data
RootNode = [0.5 0];
DT.Levels = 2;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.5;    % mm
DT.RadiusRate = 0.8;
DT.TrunkLength = 1/sqrt(2)/2; %mm
DT.LengthRate = 1/sqrt(2);

% Random tree data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 0.8;
RandomTree.Ncells = 100;
RandomTree.TerminalRadius = 0.004;

% DLA tree data
DLA.Nparticles = 250;
DLA.RootRadius = 0.5;
DLA.TerminalRadius = 0.004;
DLA.version = 3;

% Domain
%D = [0 5 0 5/sqrt(2)];
D = [0 1 0 1];
D_area = (D(2)-D(1))*(D(4)-D(3));


% Make tree
for iter1 = 1:iterations1
    RandomTree.Ncells = 100;
    for iter = 1:iterations
    CutLevel = 1;
    trees = {'Deterministic','Random','Combinated','DLA'};
    ChosenTree = trees{3};
    Tree = ChooseTree(ChosenTree,RandomTree,DT,DLA,D);
    if strcmp(ChosenTree,'Deterministic') || strcmp(ChosenTree,'Combinated')
        VisibleNodes = Tree.nodes(1:2^(CutLevel),:);
        VisibleEdges = Tree.edges(1:2^(CutLevel)-1,:);
        VisibleTree.nodes = VisibleNodes;
        VisibleTree.edges = VisibleEdges;
        VisibleTree.RootNodeIdx = 1;
    elseif strcmp(ChosenTree,'DLA')
        VisibleTree.nodes = Tree.nodes(1:2,:);
        VisibleTree.edges = Tree.edges(1,:);
        VisibleTree.RootNodeIdx = 1;
        InvisibleTree.nodes = Tree.nodes(2:end,:);
        InvisibleEdges = Tree.edges(2:end,:);
        InvisibleEdges(:,2)=InvisibleEdges(:,2)-1;
        InvisibleEdges(:,3)=InvisibleEdges(:,3)-1;
        InvisibleTree.edges = InvisibleEdges;
        InvisibleTree.RootNodeIdx = 1;
    end
    
    nodes = Tree.nodes; edges = Tree.edges;

    % Draw tree
%     if iter1 == iterations1 && iter == iterations
%         DrawTree(Tree,20,[0.8500, 0.3250, 0.0980],D);
%     end

    % Find terminal nodes
    [TNinfo,TNlogic]=FindTerminals(Tree);
    MicroTermIndexes = find(TNlogic==1);
    MacroTermIndexes = find(nodes(:,3)==CutLevel);
    
    %%% Voronoi diagram %%%
    [cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');

    %%% Parameters %%%
    mu = 3E-6;                                        % viscosity (Ns/mm^2)
    if iter1 == 1
        k = 3E-4;                                     % permeability (mm^2)
    elseif iter1 == 2
        k = 3E-5;
    elseif iter1 == 3
        k = 3E-6;
    elseif iter1 == 4
        k = 3E-7;
    elseif iter1 == 5
        k = 3E-8;
    end
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
    [p_darcyEx,q_networkEx,p_networkEx]=SolveSystemEx(Tree,TNinfo,TNlogic,Dir_network,Neu_network,mu,k,K_N,LHS,RHS,cell_area,flag.case);
    
    %save('PressurePlotDataRRT','p_darcyEx','p_networkEx','Tree','TNinfo','D','cell_center','boundary_cells','Bv_darcy');

    %SystemTestEx(Tree,TNlogic,TNinfo,mu,k,p_darcyEx,p_networkEx,q_networkEx,Grad_D,boundary_cells,cell_area,D_bvs,K_N,MacroTermIndexes,MicroTermIndexes)
    
    %% Find K_T %%%
    [q_network,p_network]=SolveNetwork(Tree,K_N,TNlogic,-1);
    [K_T,connections] = findKT(edges,cell_area,MicroTermIndexes,MacroTermIndexes,TNinfo,q_network,p_network,mu);

    % Solve system with exact K^T value
    %[TNinfoDT,TNlogicDT]=FindTerminals(VisibleTree);
    %[p_darcyCoarse,q_T,q_networkCoarse,p_networkCoarse]=SolveSystemUpS(VisibleTree,TNlogicDT,TNinfoDT,Dir_network,Neu_network,mu,k,K_N,K_T,connections,LHS,RHS,cell_area,flag.case);

    %%% Use linear regression %%%
    newK_T = zeros(size(K_T));
    R = sqrt((TNinfo(:,1)-nodes(MacroTermIndexes,1)).^2+(TNinfo(:,2)-nodes(MacroTermIndexes,2)).^2);
    coeff = load('./Filer/Coefficients_RRT');
    newK_T = 10.^(coeff.a_vect(iter)*R+coeff.b_vect(iter));

    % Solve system with kT function
    [TNinfoDT,TNlogicDT]=FindTerminals(VisibleTree);

    [p_darcyCoarse1,q_T,q_network,p_network]=SolveSystemUpS(VisibleTree,TNlogicDT,TNinfoDT,Dir_network,Neu_network,mu,k,K_N,newK_T,connections,LHS,RHS,cell_area,flag.case);
    
    %%% Root mean square error %%%
    %e = 0;
    e1 = 0;
    for i = 1:length(p_darcyEx)
        %e = e + (p_darcyEx(i)-p_darcyCoarse(i))^2*cell_area(i);
        e1 = e1 + (p_darcyEx(i)-p_darcyCoarse1(i))^2*cell_area(i);
    end
    %error(iter,iter1) = sqrt(e/D_area);
    error1(iter,iter1) = sqrt(e1/D_area);
    x(iter,iter1)=RandomTree.Ncells;
    RandomTree.Ncells = RandomTree.Ncells*2;
    end
    disp(iter1);
end

%%% Error plot DLA tree %%%
e_rel = error1/(Dir_network-Bv_darcy);

%save('./Filer/DLAerror_KD_and_Nparticles_3','e_rel','x');

save('./Filer/RRTerror_KD_ncells_7','e_rel','x');
%save('ToDrawDLAtree','Tree','VisibleTree','cells','vertices','D','MicroTermIndexes','MacroTermIndexes','TN','K_T')
%save ('PressurePlotDataDLA','D','cell_center','boundary_cells','p_darcyEx','Bv_darcy','p_darcyCoarse');
%Plots('PressureAndDeviation');
%Plots('PressurePlotNetworkDarcy');
 
 
