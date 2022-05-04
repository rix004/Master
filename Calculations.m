%clc;
clear;
close all
RootNode = [0 -1];

% Deterministic tree data
DT.Levels = 5;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 42.25;
DT.TrunkRadius = 0.1;    % mm
DT.RadiusRate = 0.7;
DT.TrunkLength = 0.7; %mm
DT.LengthRate = 0.65;

% Random tree data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 0.5;
Np = 200;
D = [-1 1 -1 1];

% Make tree
trees = {'Deterministic','Random','Combinated','Half Deterministic'};
flag.case = trees{3};
Tree = ChooseTree(flag.case,RandomTree,DT,D,Np);
nodes = Tree.nodes; edges = Tree.edges;
nodes = FindLevels(nodes,edges);

% Find terminal nodes
[MicroTN,MicroTe,TermNodes]=FindTerminals(nodes,edges);
MicroTermIndexes = find(TermNodes==1);
MacroTermIndexes = find(nodes(:,3)==DT.Levels);

%%% Voronoi diagram %%%
[cells, vertices] = VoronoiDiagram(MicroTN,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');

Ne = length(edges);
Nn = length(nodes);
Nin = sum(1-TermNodes)-1;                        % Interior nodes
Npn = sum(1-TermNodes);                          % Parent nodes
Ntn = size(MicroTN,1);                            % Terminal nodes (micro)
mu = 3E-6;                                        % viscosity (Ns/mm^2)
k = 3E-6;                                         % permeability (mm^2)
K_N(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));   % Conductance
%K_N(Te(1:end-1))=0;
f = @(x,y) 0;
K_D = @(x,y) k/mu;

%%% Set boundary values %%%
BCs = {'Dirichlet','Neumann'};
Neu_network = 5;
Dir_network = 1;
Bv_darcy = -1;



%%% Solve coupled system %%%

% Construct matrices
[Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells,bv] = TPFA(cells,vertices,f,K_D,1,Bv_darcy,edges(MicroTN(end,3),4));
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]);

I = spdiags(ones(Ntn,1),0,Ntn,Ntn);

flag.case = BCs{1};
if strcmp(flag.case, 'Dirichlet')
    Bc_nodes = [1;zeros(Nn-1,1)];
    Term_edges = zeros(Ne,1);
    Term_edges(MicroTe)=1;
    connections_trimmed = connections(:,(TermNodes+Bc_nodes)==0);

    P = zeros(Ntn,Ne);
    for i = 1:size(P,1)
        P(i,MicroTN(i,3))=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(MicroTN(i,3),4));
    end
    P = sparse(P);

    PM_vect = zeros(Ntn,1);
    for i = 1:size(PM_vect,1)
        PM_vect(i)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(MicroTN(i,3),4));
    end
    PM = spdiags(PM_vect,0,Ntn,Ntn);

    C = zeros(Ntn,Ne);
    for i = 1:size(C,1)
        C(i,MicroTN(i,3))=1;
    end
    C = sparse(C);

    I2 = zeros(Ne,Ntn);
    for i = 1:size(I2,2)
        I2(MicroTe(i),i)=1;
    end

    A = [LHS, sparse(Ntn,Ne), sparse(Ntn,Nin), sparse(Ntn,Ntn),-I;
    sparse(Ne,Ntn), spdiags(K_N.^-1,0,Ne,Ne), connections_trimmed, I2, sparse(Ne,Ntn);
    sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Nin), sparse(Nin,Ntn), sparse(Nin,Ntn);
    sparse(Ntn,Ntn), C, sparse(Ntn,Nin), sparse(Ntn,Ntn), -I;
    I, sparse(Ntn,Ne), sparse(Ntn,Nin), -I, PM];
    
    rhs = [RHS; -connections*Dir_network*Bc_nodes; zeros(Nin,1); zeros(Ntn,1); zeros(Ntn,1)];

    SOL = A\rhs;

    p_darcy = SOL(1:Ntn);
    q_network = SOL(Ntn+1:Ntn+Ne);
    p_int = SOL(Ntn+Ne+1:Ntn+Ne+Nin);
    p_tn = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);
    q_darcy = SOL(Ntn+Ne+Nin+Ntn+1:end);

    int_indexes = find(Bc_nodes+TermNodes==0);
    p_network = Dir_network*Bc_nodes;
    p_network(int_indexes) = p_int;
    p_network(MicroTermIndexes) = p_tn;

elseif strcmp(flag.case, 'Neumann')
    Bc_edges = [1;zeros(Ne-1,1)];
    connections_trimmed = connections(Bc_edges==0,(TermNodes+[1;zeros(Nn-1,1)])==0);
    connections_trimmed1 = connections(:,TermNodes==0);
    
    I_term = zeros(Ntn,Nn);
    for i = 1:size(I_term,1)
        I_term(i,MicroTermIndexes(i))=1;
    end

    P = zeros(Ntn,Ne-1);
    for i = 1:size(P,1)
        P(i,MicroTN(i,3)-1)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(MicroTN(i,3),4));
    end
    P = sparse(P);

    PM_vect = zeros(Ntn,1);
    for i = 1:size(PM_vect,1)
        PM_vect(i)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(MicroTN(i,3),4));
    end
    PM = spdiags(PM_vect,0,Ntn,Ntn);
    
    C= zeros(Ntn,Ne-1);
    for i = 1:size(C,1)
        C(i,MicroTN(i,3)-1)=1;
    end
    C = sparse(C);

    I_q = zeros(Ne,Ne-1);
    I_q(2:end,1:end)=spdiags(ones(Ne-1,1),0,Ne-1,Ne-1);

    A = [LHS, sparse(Ntn,Ne-1), sparse(Ntn,Npn), sparse(Ntn,Ntn), -I;
        sparse(Ne,Ntn), I_q, K_N.*connections, sparse(Ne,Ntn);
        sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Npn), sparse(Nin,Ntn), sparse(Nin,Ntn);
        sparse(Ntn,Ntn), C, sparse(Ntn,Npn), sparse(Ntn,Ntn), -I;
        I, sparse(Ntn,Ne-1), -I_term, PM];
    rhs = [RHS; -Neu_network*Bc_edges; [-Neu_network;zeros(Nin-1,1)]; zeros(Ntn,1); zeros(Ntn,1)];


    SOL = A\rhs;
    
    p_darcy = SOL(1:Ntn);
    q_network = [Neu_network;SOL(Ntn+1:Ntn+Ne-1)];
    p_pn = SOL(Ntn+Ne:Ntn+Ne+Nin);
    p_tn = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);
    q_darcy = SOL(Ntn+Ne+Nin+Ntn+1:end);

    p_network = [p_pn;p_tn];
end

K_T = findKT(edges,cells,MicroTermIndexes,MacroTermIndexes,q_network,p_network,p_darcy);

figure
for i = 1:length(MacroTermIndexes)
a = num2str(i);
b = 'Macro terminal ';
subplot(2,length(MacroTermIndexes)/2,i)
IntensityMap(cells,vertices,K_T(:,i),[b a])
plot(nodes(MacroTermIndexes(i),1),nodes(MacroTermIndexes(i),2),'r*')
hold on
title([b a])
axis(D)
end

rMax = 0;
KTmax = 0;
figure
for i = 1:length(MacroTermIndexes)
a = num2str(i);
b = 'Macro terminal ';
subplot(2,length(MacroTermIndexes)/2,i)
ChildCells=find(K_T(:,i)>0);
r = sqrt((cell_center(ChildCells,1)-nodes(MacroTermIndexes(i),1)).^2+(cell_center(ChildCells,2)-nodes(MacroTermIndexes(i),2)).^2);
loglog(r,K_T(K_T(:,i)>0,i),'.','MarkerSize',15)
title([b a])
xlabel('|x-x_{I}|')
ylabel('K^{T}_{I}')
hold on
if max(r)>rMax
    rMax = max(r);
end
if max(max(K_T))>KTmax;
    KTmax = max(max(K_T));
end
end
