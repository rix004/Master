%clc;
clear;
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
Np = 50;
D = [-1 1 -1 1];

trees = {'Deterministic','Random','Combinated','Half Deterministic'};
flag.case = trees{3};
Tree = ChooseTree(flag.case,RandomTree,DT,D,Np);
nodes = Tree.nodes; edges = Tree.edges;
Ne = length(edges);
Nn = length(nodes);

% Find terminal nodes
[Tn,Te,Term_nodes]=FindTerminals(nodes,edges);

%%% Set boundary values %%%
Dir_network = [1;zeros(Nn-1,1)]; % Only BV on root node
Neu_network = [1;zeros(Ne-1,1)];
Bv_darcy = -1;

%%% Voronoi diagram %%%
[cells, vertices] = VoronoiDiagram(Tn,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');

%%% Solve coupled system %%%
Nin = sum(1-Term_nodes-Dir_network);                           % Interior nodes
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
connections_trimmed = connections(:,(Term_nodes+Dir_network)==0);

I = spdiags(ones(Ntn,1),0,Ntn,Ntn);
q_term = zeros(Ntn,1);
q_term(Te)=1;
Q_term = spdiags(q_term,0,Ntn,Ne);

P = zeros(Ntn,Ne);
for i = 1:size(P,1)
    P(i,Tn(i,3))=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
end
P = sparse(P);

PM_vect = zeros(Ntn,1);
    for i = 1:size(PM_vect,1)
        PM_vect(i)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
    end
PM = spdiags(PM_vect,0,Ntn,Ntn);

C = zeros(Ntn,Ne);
for i = 1:size(C,1)
    C(i,Tn(i,3))=1;
end
C = sparse(C);

A = [LHS, sparse(Ntn,Ne), sparse(Ntn,Nin), sparse(Ntn,Ntn),-I;
    sparse(Ne,Ntn), spdiags(K_N.^-1,0,Ne,Ne), connections_trimmed, sparse(Ne,Ntn), sparse(Ne,Ntn);
    sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Nin), sparse(Nin,Ntn), sparse(Nin,Ntn);
    sparse(Ntn,Ntn), C, sparse(Ntn,Nin), sparse(Ntn,Ntn), -I;
    I, sparse(Ntn,Ne), sparse(Ntn,Nin), -I, PM];
rhs = [RHS; -connections*Dir_network; zeros(Nin,1); zeros(Ntn,1); zeros(Ntn,1)];

SOL = A\rhs;

p_darcy = SOL(1:Ntn);
q_network = SOL(Ntn+1:Ntn+Ne);
p_int = SOL(Ntn+Ne+1:Ntn+Ne+Nin);
p_tn = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);
darcy_source = SOL(Ntn+Ne+Nin+Ntn+1:end);


int_indexes = find((Term_nodes+Dir_network)==0);
term_indexes = find(Term_nodes==1);
p_network = Dir_network;
p_network(int_indexes) = p_int;
p_network(term_indexes) = p_tn;


figure('Name','Pressure plot')
for i = 1:Ne
    x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
    y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
    pressure = [p_network(edges(i,2)) p_network(edges(i,3))];
    plot3(x,y,pressure,'.-','Color',[0 0 0],'LineWidth',edges(i,4)*10)
    hold on
    grid on
end
p_map = autumn(Nn+Ntn*3);
network_map = p_map(Ntn*2:end,:);
for i = 1:Nn
    x = nodes(i,1);
    y = nodes(i,2);
    pressure = p_network(i);
    p_sorted = sort(p_network);
    ind = find(p_sorted==pressure);
    ind = ind(1);
    plot3(x,y,pressure,'.','Color',network_map(ind,:),'MarkerSize',20);
    hold on
end
zlabel('Node pressure')

darcy_map = p_map(1:3*Ntn,:);

[xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
p_points = [Tn(:,1:2);boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
p_values = [p_darcy;Bv_darcy*ones(size(boundary_cells,1)+4,1)];
vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
surf(xq,yq,vq)
colormap(darcy_map)

% Check that the flux out of the domain (over dOmega) is the same as the
% flux coming into the domain and the flux into the network
sum(Grad_D(boundary_cells(:,3),:)*p_darcy-D_bvs(boundary_cells(:,3)))
sum(darcy_source)
q_network(1)

%% Check that Poiseuilles law is also valid in terminal edges
disp('Poiseuilles law in terminal edge (should be zero)')
for i = 1:length(Te)
    sum_PL(i) = K_N(Te(i))*(p_network(edges(Te(i),2))-p_network(edges(Te(i),3)))-q_network(Te(i));
end
sum(sum_PL)
                
