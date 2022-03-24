close all
clc; clear;
RootNode = [0.5 0];

% Deterministic tree data
DT.Levels = 4;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.05;    
DT.RadiusRate = 0.5;
DT.TrunkLength = 0.5;
DT.LengthRate = 0.65;

% Random tree data
RandomTree.TrunkRadius = 0.1;
RandomTree.RadiusRate = 0.5;
Np = 300;
Domain = [0 1 0 1];

trees = {'Deterministic','Random','Combinated','Half Deterministic'};
flag.case = trees{3};
Tree = ChooseTree(flag.case,RandomTree,DT,Domain,Np);
DrawTree(Tree);
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

%%% Find distance from nodes to start node %%%
dist = sqrt((nodes(:,1)-RootNode(1)).^2+(nodes(:,2)-RootNode(2)).^2);
dist_tn = dist(EndNodes(find_tn==0));

%%% Set boundary values %%%
Bc_vals = [1;zeros(Nn-1,1)];
for i = 2:length(Bc_vals)
    Bc_vals(i)=Bc_nodes(i)*0;
end

%%% Voronoi diagram %%%
[cells, vertices] = VoronoiDiagram(Tn,[0 0 1 1 0;0 1 1 0 0]');

%%% Solve coupled system %%%
Nin = sum(1-Bc_nodes);  % Interior nodes
Ntn = size(Tn,1);       % Terminal nodes
mu = 5.5;               % viscosity (mPa*s)
g(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));   % Conductance
f = @(x,y) 1;
k = @(x,y) 1;

% Construct matrices
[LHS,~,cell_center,cell_edges,cell_area] = TPFA(cells,vertices,f,k,1,0);
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]);
connections_trimmed = connections(:,Bc_nodes==0);

I = spdiags(ones(Ntn,1),0,Ntn,Ntn);
q_term = zeros(Ntn,1);
q_term(Te)=1;
Q_term = spdiags(q_term,0,Ntn,Ne);

B = zeros(Ntn,Ne);
for i = 1:size(B,1)
    B(i,Tn(i,3))=1/(2*pi*k(cell_center(i,1),cell_center(i,2)))*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
end
B = sparse(B);

C = zeros(Ntn,Ne);
for i = 1:size(C,1)
    C(i,Tn(i,3))=1;
end
C = sparse(C);

A = [LHS, sparse(Ntn,Ne), sparse(Ntn,Nin), -I, sparse(Ntn,Ntn);
    sparse(Ne,Ntn), spdiags(g.^-1,0,Ne,Ne), connections_trimmed, sparse(Ne,Ntn), sparse(Ne,Ntn);
    sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Nin), sparse(Nin,Ntn), sparse(Nin,Ntn);
    sparse(Ntn,Ntn), C, sparse(Ntn,Nin) -I, sparse(Ntn,Ntn);
    I, -B, sparse(Ntn,Nin), sparse(Ntn,Ntn), -I];
rhs = [zeros(size(LHS,1),1); -connections*Bc_vals; zeros(Nin,1); zeros(Ntn,1); zeros(Ntn,1)];

SOL = A\rhs;

p_darcy = SOL(1:Ntn);
q_network = SOL(Ntn+1:Ntn+Ne);
p_int = SOL(Ntn+Ne+1:Ntn+Ne+Nin);
darcy_source = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);
p_tn = SOL(Ntn+Ne+Nin+Ntn+1:end);

indexes = find(Bc_nodes==0);
p_all = Bc_vals;
p_all(indexes) = p_int;


%%% PLOTS %%%%
IntensityMap(cells,vertices,q_network(Te),'Flux in terminal network edges')
plot(Tn(:,1),Tn(:,2),'b.')
hold on
plot(RootNode(1),RootNode(2),'r*');
hold on
axis([0 1 0 1])

figure('Name','Pressure plot for Darcy domain')
plot3(cell_center(:,1),cell_center(:,2),p_darcy,'.','MarkerSize',20);
grid on
IntensityMap(cells,vertices,p_darcy,'Pressure in Darcy domain');
axis([0 1 0 1])
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
zlabel('p','FontSize',15)

% figure('Name','Pressure plot')
% for i = 1:Ne
%     x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
%     y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
%     pressure = [p_all(edges(i,2)) p_all(edges(i,3))];
%     plot3(x,y,pressure,'.-','Color',[0 0 0],'LineWidth',edges(i,4)*10)
%     hold on
%     grid on
% end
% p_map = autumn(Nn);
% for i = 1:Nn
%     x = nodes(i,1);
%     y = nodes(i,2);
%     pressure = p_all(i);
%     p_sorted = sort(p_all,'descend');
%     ind = find(p_sorted==pressure);
%     ind = ind(1);
%     plot3(x,y,pressure,'.','Color',p_map(ind,:),'MarkerSize',20);
%     hold on
% end
% zlabel('Node pressure')
% 
% figure('Name','Flux in terminal edges vs. distance to root node');
% loglog(dist_tn,q(Te),'*')
% xlabel('Distance (L)')
% ylabel('Flux (L^3T^{-1}L^{-2})')



                
                
