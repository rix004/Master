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
Np = 10;
D = [-1 1 -1 1];

trees = {'Deterministic','Random','Combinated','Half Deterministic'};
flag.case = trees{3};
Tree = ChooseTree(flag.case,RandomTree,DT,D,Np);
nodes = Tree.nodes; edges = Tree.edges;
nodes = FindLevels(nodes,edges);


% Find terminal nodes
[Tn,Te,Term_nodes]=FindTerminals(nodes,edges);


%%% Voronoi diagram %%%
[cells, vertices] = VoronoiDiagram(Tn,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');

%%% Set boundary values %%%
BCs = {'Dirichlet','Neumann'};
Bv_darcy = -1;
Neu_network = 5;
Dir_network = 1;

% %%% Solve coupled system %%%
Ne = length(edges);
Nn = length(nodes);
Nin = sum(1-Term_nodes)-1;                        % Interior nodes
Npn = sum(1-Term_nodes);                          % Parent nodes
Ntn = size(Tn,1);                                 % Terminal nodes
Nc = size(cells,1);
mu = 3E-6;                                        % viscosity (Ns/mm^2)
k = 3E-6;                                         % permeability (mm^2)
K_N(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));   % Conductance
f = @(x,y) 0;
K_D = @(x,y) k/mu;

% Construct matrices
[darcy_vertices,~,triangles]= GridGeneration(10,10,D);
[Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells,bv] = TPFA(triangles,darcy_vertices,f,K_D,1,Bv_darcy);

O = Overlap(cells,vertices,triangles,darcy_vertices);
sumArea = zeros(size(O,1),1);
for i = 1:size(O,1)
    sumArea(i) = sum(O(i,:));
end
% connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]);
% 
% I1 = spdiags(ones(Ntn,1),0,Ntn,Ntn);
% 
% flag.case = BCs{1};
% if strcmp(flag.case, 'Dirichlet')
%     Bc_nodes = [1;zeros(Nn-1,1)];
%     Term_edges = zeros(Ne,1);
%     Term_edges(Te)=1;
%     connections_trimmed = connections(:,(Term_nodes+Bc_nodes)==0);
% 
%     Areas = spdiags(sumArea,0,Nc,Nc);
% 
%     I2 = zeros(Ne,Ntn);
%     for i = 1:size(I2,2)
%         I2(Te(i),i)=1;
%     end
% 
%     I3 = spdiags(ones(Nc,1),0,Nc,Nc);
% 
%     C = zeros(Ntn,Ne);
%     for i = 1:size(C,1)
%         C(i,Tn(i,3))=1;
%     end
%     C = sparse(C);
%     
%     TermToCell = zeros(size(O));
%     TermToCell(O>0)=1;
% 
% 
%     A = [LHS, sparse(Nc,Ne), sparse(Nc,Nin), sparse(Nc,Ntn),-I3;
%     sparse(Ne,Nc), spdiags(K_N.^-1,0,Ne,Ne), connections_trimmed, I2, sparse(Ne,Nc);
%     sparse(Nin,Nc), connections_trimmed', sparse(Nin,Nin), sparse(Nin,Ntn), sparse(Nin,Nc);
%     -Areas, zeros(Nc,Ne), sparse(Nc,Nin), O, -I3;
%     sparse(Ntn,Nc), C, sparse(Ntn,Nin), sparse(Ntn,Ntn), -TermToCell'];
%     
%     rhs = [RHS; -connections*Dir_network*Bc_nodes; zeros(Nin,1); zeros(Ntn,1); zeros(Nc,1)];
% 
%     SOL = A\rhs;
% 
%     p_darcy = SOL(1:Nc);
%     q_network = SOL(Nc+1:Nc+Ne);
%     p_int = SOL(Nc+Ne+1:Nc+Ne+Nin);
%     p_tn = SOL(Nc+Ne+Nin+1:Nc+Ne+Nin+Ntn);
%     q_darcy = SOL(Nc+Ne+Nin+Ntn+1:end);
% 
%     int_indexes = find(Bc_nodes+Term_nodes==0);
%     p_network = Dir_network*Bc_nodes;
%     p_network(int_indexes) = p_int;
%     p_network(Term_nodes==1) = p_tn;
% 
% elseif strcmp(flag.case, 'Neumann')
%     Bc_edges = [1;zeros(Ne-1,1)];
%     connections_trimmed = connections(Bc_edges==0,(Term_nodes+[1;zeros(Nn-1,1)])==0);
%     connections_trimmed1 = connections(:,Term_nodes==0);
%     
%     I_term = zeros(Ntn,Nn);
%     for i = 1:size(I_term,1)
%         I_term(i,term_indexes(i))=1;
%     end
% 
%     P = zeros(Ntn,Ne-1);
%     for i = 1:size(P,1)
%         P(i,Tn(i,3)-1)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
%     end
%     P = sparse(P);
% 
%     PM_vect = zeros(Ntn,1);
%     for i = 1:size(PM_vect,1)
%         PM_vect(i)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
%     end
%     PM = spdiags(PM_vect,0,Ntn,Ntn);
%     
%     C= zeros(Ntn,Ne-1);
%     for i = 1:size(C,1)
%         C(i,Tn(i,3)-1)=1;
%     end
%     C = sparse(C);
% 
%     I_q = zeros(Ne,Ne-1);
%     I_q(2:end,1:end)=spdiags(ones(Ne-1,1),0,Ne-1,Ne-1);
% 
%     A = [LHS, sparse(Ntn,Ne-1), sparse(Ntn,Npn), sparse(Ntn,Ntn), -I;
%         sparse(Ne,Ntn), I_q, K_N.*connections, sparse(Ne,Ntn);
%         sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Npn), sparse(Nin,Ntn), sparse(Nin,Ntn);
%         sparse(Ntn,Ntn), C, sparse(Ntn,Npn), sparse(Ntn,Ntn), -I;
%         I, sparse(Ntn,Ne-1), -I_term, PM];
%     rhs = [RHS; -Neu_network*Bc_edges; [-Neu_network;zeros(Nin-1,1)]; zeros(Ntn,1); zeros(Ntn,1)];
% 
% 
%     SOL = A\rhs;
%     
%     p_darcy = SOL(1:Ntn);
%     q_network = [Neu_network;SOL(Ntn+1:Ntn+Ne-1)];
%     p_pn = SOL(Ntn+Ne:Ntn+Ne+Nin);
%     p_tn = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);
%     q_darcy = SOL(Ntn+Ne+Nin+Ntn+1:end);
% 
%     p_network = [p_pn;p_tn];
%  end

% figure('Name','Pressure plot')
% for i = 1:Ne
%     x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
%     y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
%     pressure = [p_network(edges(i,2)) p_network(edges(i,3))];
%     plot3(x,y,pressure,'.-','Color',[0 0 0],'LineWidth',edges(i,4)*10)
%     hold on
%     grid on
% end
% % p_map = autumn(Nn+Nc*3);
% % network_map = p_map(Ntn*2:end,:);
% for i = 1:Nn
%     x = nodes(i,1);
%     y = nodes(i,2);
%     pressure = p_network(i);
% %     p_sorted = sort(p_network);
% %     ind = find(p_sorted==pressure);
% %     ind = ind(1);
%     plot3(x,y,pressure,'.','Color',[0 0 0],'MarkerSize',20);
%     hold on
% end
% zlabel('Node pressure')
% % darcy_map = p_map(1:3*Ntn,:);
% 
% [xq,yq]=meshgrid(D(1):0.1:D(2), D(3):0.1:D(4));
% p_points = [cell_center(:,1:2);boundary_cells(:,1:2)];
% p_values = [p_darcy;bv];
% vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
% surf(xq,yq,vq)
% %colormap(darcy_map)
% 
% % Check that the flux out of the domain (over dOmega) is the same as the
% % flux coming into the domain and the flux into the network
% disp('Flux out of domain, into domain and into network (should be the same):')
% sum(Grad_D(boundary_cells(:,3),:)*p_darcy-D_bvs(boundary_cells(:,3)))
% sum(q_darcy)
% q_network(1)
% 
% % Check that the sum of fluxes in and out of an interior node equals 0.
% int_nodes = find((Term_nodes+[1;zeros(Nn-1,1)])==0);
% total_flux = zeros(Nin,1);
% for i =1:length(int_nodes)
%     edge_in=find(edges(:,3)==int_nodes(i));
%     edge_out=find(edges(:,2)==int_nodes(i));
%     flux_in = q_network(edge_in);
%     flux_out = q_network(edge_out);
%     total_flux(i) = sum(flux_in)-sum(flux_out);
% end
% disp('Sum flux in and out of interior nodes (should be zero):')
% sum(total_flux)





