%clc;
clear;
close all
function[l2_error,h]=main(Np)
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
Ne = length(edges);
Nn = length(nodes);

% Find terminal nodes
[Tn,Te,Term_nodes]=FindTerminals(nodes,edges);
term_indexes = find(Term_nodes==1);

%%% Set boundary values %%%
BCs = {'Dirichlet','Neumann'};
Bv_darcy = -1;
Neu_network = 5;
Dir_network = 1;

%%% Voronoi diagram %%%
[cells, vertices] = VoronoiDiagram(Tn,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');

%%% Solve coupled system %%%
Nin = sum(1-Term_nodes)-1;                        % Interior nodes
Npn = sum(1-Term_nodes);                          % Parent nodes
Ntn = size(Tn,1);                                 % Terminal nodes
mu = 3E-6;                                        % viscosity (Ns/mm^2)
k = 3E-6;                                         % permeability (mm^2)
K_N(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));   % Conductance
K_N(Te(1:end-1))=0;
f = @(x,y) 0;
K_D = @(x,y) k/mu;

% Construct matrices
[Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells,bv] = TPFA(cells,vertices,f,K_D,1,Neu_network*mu/k,edges(Tn(end,3),4));
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]);

I = spdiags(ones(Ntn,1),0,Ntn,Ntn);

flag.case = BCs{2};
if strcmp(flag.case, 'Dirichlet')
    Bc_nodes = [1;zeros(Nn-1,1)];
    Term_edges = zeros(Ne,1);
    Term_edges(Te)=1;
    connections_trimmed = connections(:,(Term_nodes+Bc_nodes)==0);

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

    I2 = zeros(Ne,Ntn);
    for i = 1:size(I2,2)
        I2(Te(i),i)=1;
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

    int_indexes = find(Bc_nodes+Term_nodes==0);
    p_network = Dir_network*Bc_nodes;
    p_network(int_indexes) = p_int;
    p_network(term_indexes) = p_tn;

elseif strcmp(flag.case, 'Neumann')
    Bc_edges = [1;zeros(Ne-1,1)];
    connections_trimmed = connections(Bc_edges==0,(Term_nodes+[1;zeros(Nn-1,1)])==0);
    connections_trimmed1 = connections(:,Term_nodes==0);
    
    I_term = zeros(Ntn,Nn);
    for i = 1:size(I_term,1)
        I_term(i,term_indexes(i))=1;
    end

    P = zeros(Ntn,Ne-1);
    for i = 1:size(P,1)
        P(i,Tn(i,3)-1)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
    end
    P = sparse(P);

    PM_vect = zeros(Ntn,1);
    for i = 1:size(PM_vect,1)
        PM_vect(i)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
    end
    PM = spdiags(PM_vect,0,Ntn,Ntn);
    
    C= zeros(Ntn,Ne-1);
    for i = 1:size(C,1)
        C(i,Tn(i,3)-1)=1;
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


figure('Name','Pressure plot')


[xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
p_points = [Tn(:,1:2);boundary_cells(:,1:2)];
p_values = [p_darcy;bv];
vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
surf(xq,yq,vq)

p_exact = @(x,y) -Neu_network*mu/(2*pi*k)*log(sqrt(x.^2+y.^2)/edges(Tn(i,3),4));
error = p_darcy-p_exact(cell_center(:,1),cell_center(:,2));
l2_error = 0;
for i =1:length(error)
    l2_error = l2_error + error(i)^2*cell_area(i);
end

l2_error
h=sqrt(4/size(cells,1))

figure
loglog([0.25 0.25/2 0.25/4 0.25/8],[0.1 0.05 0.025 0.025/2],'-','LineWidth',2.5,'Color','r')
hold on
loglog([0.5547 0.3592 0.2673 0.2 0.1414 0.0979 0.0695],[0.1997 0.0317 0.0125 0.0150 0.0268 0.0036 0.0013],'-.','LineWidth',2,'Color','b')
hold on

end








