%clc;
close all
clear;
iterations = 8;
iterations1 = 6;
TN = 1;
Ncells =200; 

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
DT.RootNodeIdx = 1;

% Random tree data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 0.7;

% Domain
D = [0 5 0 5/sqrt(2)];
%D = [0.2 0.8 0.2 0.8];
D_area = (D(2)-D(1))*(D(4)-D(3));

for iter1 = 1:iterations1
    for iter = 1:iterations
    % Where to cut
    CutLevel = DT.Levels-iter;
    trees = {'Deterministic','Random','Combinated','DLA'};
    ChosenTree = trees{1};
    Tree = ChooseTree(ChosenTree,RandomTree,DT,D,Ncells);
    nodes = Tree.nodes; edges = Tree.edges;
    if strcmp(ChosenTree,'Deterministic') || strcmp(ChosenTree,'Combinated')
        VisibleNodes = nodes(1:2^(CutLevel-1),:);
        VisibleEdges = edges(1:2^(CutLevel-1)-1,:);
        VisibleTree.nodes = VisibleNodes;
        VisibleTree.edges = VisibleEdges;
        VisibleTree.RootNodeIdx = 1;
    elseif ChosenTree == 'DLA'
        VisibleNodes = nodes(Tree.RootNodeIdx);
    end

    % Find terminal nodes
    [TNinfo,TNlogic]=FindTerminals(Tree);
    MicroTermIndexes = find(TNlogic==1);
    MacroTermIndexes = find(nodes(:,3)==CutLevel);
    
    % Fix edge radii
%     rel = edges(:,1)./edges(:,4);
%     fix = find(rel<20);
%     edges(fix,4)=edges(fix,1)/100;

    %%% Voronoi diagram %%%
    [cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');
    
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
    
    %%% For validation %%%%
    if strcmp(ChosenTree,'Deterministic')
        summ = 0;
        for j = CutLevel:DT.Levels-1
            summ = summ + DT.LengthRate^(j-1)/(2^(j-1)*DT.RadiusRate^(4*(j-1)));
        end
        TrueKT(iter,iter1) = pi*DT.TrunkRadius^4/(8*mu*DT.TrunkLength*D_area*summ);
    else
        TrueKT = 0;
    end

    %%% Find K^T %%%
    [K_T,connections] = findKT(edges,cell_area,MicroTermIndexes,MacroTermIndexes,TNinfo,q_network,p_network,mu);
    connections=connections';
    kT(:,iter)=K_T(:,TN);
    kTnotZero = kT(connections(:,TN)==1,iter);
    p_termnodes = p_network(MicroTermIndexes(connections(:,TN)==1));
    for it = 1:sum(connections(:,TN)==1)
        qT(it)=kTnotZero(it)*(p_network(MacroTermIndexes(TN))-p_termnodes(it));
    end

    r(:,iter) = sqrt((TNinfo(:,1)-nodes(MacroTermIndexes(TN),1)).^2+(TNinfo(:,2)-nodes(MacroTermIndexes(TN),2)).^2);
    rVSk = [r(:,iter) kT(:,iter)];
    rVSkSorted = sortrows(rVSk,1);
    r(:,iter)=rVSkSorted(:,1);
    kT(:,iter)=rVSkSorted(:,2);
    disp(iter)
    rr(iter,iter1)=DT.RadiusRate;
    levelcut(iter,iter1)=iter;
    end
    DT.RadiusRate = DT.RadiusRate + 0.1;
end
NotZeroVals = zeros(size(kT));
NotZeroVals(kT>0)=1;

%save('KTvsDeltaRData','rr','TrueKT');
%save('AnalyticVScompKT1','NotZeroVals','kT','r','TrueKT')
% Plots('TrueVScodeKT');
% Plots('KTradiusrateFIG');


%[a,b] = LinReg(kT,r,iter,MicroTermIndexes,MacroTermIndexes,TrueKT);
