clc;clear;
close all
RootNode = [0 -1];

% Deterministic tree data
DT.Levels = 4;
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
Np = 1;
D = [-1 1 -1 1];

trees = {'Deterministic','Random','Combinated','Half Deterministic'};
flag.case = trees{3};
Tree = ChooseTree(flag.case,RandomTree,DT,D,Np);
nodes = Tree.nodes; edges = Tree.edges;

% Find terminal nodes
Ne = length(edges);
Nn = length(nodes);
StartNodes = edges(:,2);
EndNodes = edges(:,3);
find_tn = ismember(EndNodes,StartNodes);
Tn = [];
Te = [];
Bc_nodes = [1;zeros(Nn-1,1)];
for i = 1:length(EndNodes)
    if find_tn(i) == 0
        Tn = [Tn;nodes(EndNodes(i),1:2),i];
        Te = [Te;i];
        Bc_nodes(EndNodes(i))=1; 
    end
end

%%% Set boundary values %%%
Bv_network = [1;zeros(Nn-1,1)]; % Only BV on root node
Bv_darcy = -1;

%%% Voronoi diagram %%%
[cells, vertices] = VoronoiDiagram(Tn,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');

%%% Solve coupled system %%%
Nin = sum(1-Bc_nodes);                            % Interior nodes
Ntn = size(Tn,1);                                 % Terminal nodes
mu = 3E-6;                                        % viscosity (Ns/mm^2)
k = 3E-6;                                         % permeability (mm^2)
K_N(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));   % Conductance
%K_N(Te(1:end-1))=0;
f = @(x,y) 0;
K_D = @(x,y) k/mu;

% Construct matrices
[Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells] = TPFA(cells,vertices,f,K_D,1,Bv_darcy);
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]);
connections_trimmed = connections(:,Bc_nodes==0);

I = spdiags(ones(Ntn,1),0,Ntn,Ntn);
q_term = zeros(Ntn,1);
q_term(Te)=1;
Q_term = spdiags(q_term,0,Ntn,Ne);

P = zeros(Ntn,Ne);
for i = 1:size(P,1)
    P(i,Tn(i,3))=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
end
P = sparse(P);

C = zeros(Ntn,Ne);
for i = 1:size(C,1)
    C(i,Tn(i,3))=1;
end
C = sparse(C);

A = [LHS, sparse(Ntn,Ne), sparse(Ntn,Nin), -I, sparse(Ntn,Ntn);
    sparse(Ne,Ntn), spdiags(K_N.^-1,0,Ne,Ne), connections_trimmed, sparse(Ne,Ntn), sparse(Ne,Ntn);
    sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Nin), sparse(Nin,Ntn), sparse(Nin,Ntn);
    sparse(Ntn,Ntn), C, sparse(Ntn,Nin) -I, sparse(Ntn,Ntn);
    I, P, sparse(Ntn,Nin), sparse(Ntn,Ntn), -I];
rhs = [RHS; -connections*Bv_network; zeros(Nin,1); zeros(Ntn,1); zeros(Ntn,1)];

SOL = A\rhs;

p_darcy = SOL(1:Ntn);
q_network = SOL(Ntn+1:Ntn+Ne);
p_int = SOL(Ntn+Ne+1:Ntn+Ne+Nin);
darcy_source = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);
p_tn = SOL(Ntn+Ne+Nin+Ntn+1:end);

int_indexes = find(Bc_nodes==0);
term_indexes = find([0;Bc_nodes(2:end)]==1);
p_network = Bv_network;
p_network(int_indexes) = p_int;
p_network(term_indexes) = p_tn;


                
