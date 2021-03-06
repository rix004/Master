%clc;
close all
clear;
TN = 1;
Ncells = 100; 
iterations = 6;

% Deterministic tree data
RootNode = [2.5 0.2];
DT.Levels = 4;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.RotationRate = 1.05;
DT.TrunkRadius = 0.05;    % mm
DT.RadiusRate = 0.7;
DT.TrunkLength = 1.8; %mm
DT.LengthRate = 0.7;

% Random tree data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 0.7;

% Domain
D = [0 5 0 5];
D_area = (D(2)-D(1))*(D(4)-D(3));

% Make tree
% h3 = figure('Name','Error plot');
% error = zeros(iterations,1);
% l = zeros(iterations,1);
%for iter = 1:iterations
CutLevel = DT.Levels;
trees = {'Deterministic','Random','Combinated','Half Deterministic'};
flag.case = trees{3};
Tree = ChooseTree(flag.case,RandomTree,DT,D,Ncells);
nodes = Tree.nodes; edges = Tree.edges;
if strcmp(flag.case, 'Combinated')
    DETnodes = nodes(1:2^(CutLevel-1),:);
    DETedges = edges(1:2^(CutLevel-1)-1,:);
elseif strcmp(flag.case, 'Deterministic')
    DETnodes = nodes(1:2^(CutLevel-1),:);
    DETedges = edges(1:2^(CutLevel-1)-1,:);
end

% Find terminal nodes
[TNinfo,TNlogic]=FindTerminals(nodes,edges);
MicroTermIndexes = find(TNlogic==1);
MacroTermIndexes = find(nodes(:,3)==CutLevel);

% Fix edge radiee
rel = edges(:,1)./edges(:,4);
fix = find(rel<20);
edges(fix,4)=edges(fix,1)/100;

%%% Voronoi diagram %%%
[cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');

%%% Parameters %%%
mu = 3E-6;                                        % viscosity (Ns/mm^2)
k = 3E-6;                                        % permeability (mm^2)
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
[p_darcy,q_T,q_network,p_network]=SolveSystemEx(nodes,edges,TNinfo,TNlogic,Dir_network,Neu_network,mu,k,K_N,LHS,RHS,cell_area,flag.case);

% Test system
%SystemTestEx(nodes,edges,TNlogic,TNinfo,mu,k,p_darcy,q_darcy,p_network,q_network,Grad_D,boundary_cells,cell_area,D_bvs,K_N,MacroTermIndexes,MicroTermIndexes)

%% Find K_T %%%
[q_network,p_network]=SolveNetwork(nodes,edges,K_N,TNlogic,-1);
[K_T,connections] = findKT1(edges,cell_area,MicroTermIndexes,MacroTermIndexes,TNinfo,q_network,p_network,mu);

% Solve system with kT function
[TNinfoDT,TNlogicDT]=FindTerminals(DETnodes,DETedges);
[p_darcy1,q_T,q_network,p_network]=SolveSystemUpS1(DETnodes,DETedges,TNlogicDT,TNinfoDT,Dir_network,Neu_network,mu,k,K_N,K_T,connections,LHS,RHS,cell_area,flag.case);

%Test system
%SystemTestCoarse(DETnodes,DETedges,TNlogicDT,TNinfoDT,mu,k,p_darcy1,q_T,p_network,q_network,Grad_D,boundary_cells,cell_area,D_bvs,K_N,MacroTermIndexes,connections)

%%% Root mean square error %%%
% e = 0;
% for i = 1:length(p_darcy)
%     e = e + (p_darcy(i)-p_darcy1(i))^2*cell_area(i);
% end
% error(iter) = sqrt(e/D_area);
% l(iter)=DT.Levels-iter;
%end

%%%% PLOTS %%%%
% figure(h3)
% plot(ones(iterations,1)*DT.Levels-l,error,'.','MarkerSize',20)
% xlabel('Number of cutted levels','FontSize',15)
% ylabel('RMS','FontSize',15,'Rotation',0)
% xticks(sort(ones(iterations,1)*DT.Levels-flip(l)))


% Pressure plot
h1 = figure('Name','Pressure in Darcy domain')
subplot(1,3,1)
[xq,yq]=meshgrid(D(1):0.1:D(2), D(3):0.1:D(4));
p_points = [cell_center;boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
p_values = [p_darcy;Bv_darcy*ones(size(boundary_cells,1)+4,1)];
vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
surf(xq,yq,vq)
title('Exact model','FontSize',18)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
zlabel('P^{D}_h','FontSize',13,'Rotation',0)

h2 = figure()
subplot(1,3,1)
imagesc(vq)
colormap gray
colorbar

figure(h1)
subplot(1,3,2)
p_points = [cell_center;boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
p_values = [p_darcy1;Bv_darcy*ones(size(boundary_cells,1)+4,1)];
vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
surf(xq,yq,vq)
title('Coarse model','FontSize',18)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
zlabel('P^{D}_h','FontSize',13,'Rotation',0)
colormap jet

figure(h2)
subplot(1,3,2)
imagesc(vq)
colormap gray
colorbar

figure(h1)
subplot(1,3,3)
p_points = [cell_center;boundary_cells(:,1:2);[D(1) D(3)];[D(1) D(4)];[D(2) D(3)]; [D(2) D(4)]];
p_values = [p_darcy1-p_darcy;Bv_darcy*zeros(size(boundary_cells,1)+4,1)];
vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
surf(xq,yq,vq)
title('Difference','FontSize',20)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
colormap jet 

figure(h2)
subplot(1,3,3)
imagesc(vq)
colormap gray
colorbar

a = 'Impact field of one terminal node';
figure
IntensityMap(cells,vertices,K_T(:,1),a)
plot(nodes(MacroTermIndexes(TN),1),nodes(MacroTermIndexes(TN),2),'r*','MarkerSize',20)
hold on
axis(D)
DrawTree(Tree,20,'b',D);
DETtree.nodes = DETnodes;
DETtree.edges = DETedges;
DrawTree(DETtree,50,[0.8500, 0.3250, 0.0980],D);
axis off

