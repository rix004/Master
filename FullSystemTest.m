% Neumann BC i "stamme"
% Riktig randbetingelser på domenet (må ordnes i TPFA)
% Konduktivitet lik 0 på alle enteterminaler utenom den i midten
% Fjern delen i ChooseTree som gjør at det produserest eksakt Ncells celler
% (for å få den siste terminalnoden akkurat der man vil, i origo)

clc;
clear;
close all
iterations1 = 10;
iterations2 = 6;
L2_error = zeros(iterations1,iterations2);
h = zeros(iterations1,iterations2);

% Deterministic tree data
RootNode = [0 -1];
DT.Levels = 4;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 60;
DT.TrunkRadius = 0.05;    % mm
DT.RadiusRate = 0.9;
DT.TrunkLength = 1/sqrt(2)/2; %mm
DT.LengthRate = 1/sqrt(2);

% Random tree data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 0.7;

% Domain
D = [-1 1 -1 1];
D_area = (D(2)-D(1))*(D(4)-D(3));

for iter1 = 1:iterations1
    Ncells = 50;
    for iter2 = 1:iterations2
    flag.case = 'Combinated';
    Tree = ChooseTree(flag.case,RandomTree,DT,D,Ncells);
    nodes = Tree.nodes; edges = Tree.edges;
%     figure()
%     DrawTree(Tree,150,'b',D);
    
    % Find terminal nodes
    [TNinfo,TNlogic]=FindTerminals(nodes,edges);
    
    % Fix edge radiee
    rel = edges(:,1)./edges(:,4);
    fix = find(rel<20);
    edges(fix,4)=edges(fix,1)/100;
    fix = find(edges(:,4)<1E-6);
    edges(fix,4)=1E-6;
    
    %%% Voronoi diagram %%%
    [cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');
    
    %%% Parameters %%%
    k = 3E-6;
    mu = 3E-6;
    K_N = (pi*edges(:,4).^4./(8*mu*edges(:,1)));      % Conductance
    K_N(TNinfo(1:end-1,3))=0;
    K_D = k/mu;                                       % Hydraulic conductivity [mm^4/Ns)
    f = @(x,y) 0;                                     % External sources or sinks
    
    %%% Set boundary values %%%
    BC = 'Neumann';
    Neu_network = 5;
    Dir_network = 1;
    Bv_darcy = -1;
    
    %%% Solve coupled system with stocastic tree %%%
    [Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells,bv] = TPFA(cells,vertices,f,K_D,1,Neu_network,edges(TNinfo(end,3),4));

   [p_darcy,~,~,~]=SolveSystemEx(nodes,edges,TNinfo,TNlogic,Dir_network,Neu_network,mu,k,K_N,LHS,RHS,cell_area,BC);

    %%%% Test with fundamental solution %%%%%%
    p_exact = @(x,y) -Neu_network*mu/(2*pi*k)*log(sqrt(x.^2+y.^2)/edges(TNinfo(end,3),4));
    epsilon = 0.1;
    distToOrigin = sqrt(cell_center(:,1).^2+(cell_center(:,2).^2));
    ok_cells = find(distToOrigin>epsilon);
    notok_cell = find(distToOrigin<epsilon);
    cells_x = cell_center(ok_cells,1);
    cells_y = cell_center(ok_cells,2);
    error = p_darcy(ok_cells)-p_exact(cells_x,cells_y);
    l2 = 0;
    for i =1:length(error)
        l2 = l2 + error(i)^2*cell_area(i);
    end
    l2 = sqrt(l2);
    L2_error(iter1,iter2) = l2;
    h(iter1,iter2) = sqrt(D_area/size(cells,1));
    Ncells = Ncells*2;
    end
    disp(iter1)
end

save('L2errorDataSystemTest','L2_error','h');
Plots('TotalSystemTest')