close all
clc; clear;
RootNode = [0 0];

% Deterministic tree data
DT.Levels = 5;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.1;    
DT.RadiusRate = 0.5;
DT.TrunkLength = 1;
DT.LengthRate = 0.65;

% Random tree data
RandomTree.TrunkRadius = 0.1;
RandomTree.RadiusRate = 0.5;
Np = 50;
Domain = [-2 2 0 3];

trees = {'Deterministic','Random','Combinated','Half Deterministic'};
flag.case = trees{1};

if strcmp(flag.case, 'Deterministic')
    Tree= GetTree(DT);
elseif strcmp(flag.case, 'Random')
    RandomTree.nodes = [RootNode 1];
    RandomTree.edges = [];
    RandomTree = RRT_Tree(RandomTree,Domain,Np,0);
    Tree = RandomTree;
elseif strcmp(flag.case, 'Combinated')
    DetTree= GetTree(DT);
    RandomTree.nodes = DetTree.nodes;
    RandomTree.edges = DetTree.edges;
    Tree = RRT_Tree(RandomTree,Domain,Np,0);
elseif strcmp(flag.case, 'Half Deterministic')
    DetTree= GetTree(DT);
    pos_nodes = find(DetTree.nodes(:,1)>=0);
    neg_nodes = find(DetTree.nodes(:,1)<0);
    half_nodes = DetTree.nodes(pos_nodes,:);
    pos_edges = ismember(DetTree.edges(:,3),pos_nodes);
    half_edges = DetTree.edges(pos_edges==1,:);
    for i = 1:length(half_edges)
        ii = find(pos_nodes==half_edges(i,3));
        jj = find(pos_nodes==half_edges(i,2));
        half_edges(i,3)=ii;
        half_edges(i,2)=jj;
    end
    Tree.nodes = half_nodes;
    Tree.edges = half_edges;
end

DrawTree(Tree);
nodes = Tree.nodes; edges = Tree.edges;


%%% CALCULATE PRESSURE AND FLUX IN NETWORK %%%
Ne = length(edges); Nn = length(nodes);
mu = 5.5; % viscosity (mPa*s)

% Conductance
g(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));

% Find terminal nodes
StartNodes = edges(:,2);
EndNodes = edges(:,3);
find_tn = ismember(EndNodes,StartNodes);
Tn = [];
Te = [];
Bc_nodes = [1;zeros(Nn-1,1)];
for i = 1:length(EndNodes)
    if find_tn(i) == 0
        Tn = [Tn;nodes(EndNodes(i),:)];
        Te = [Te;i];
        Bc_nodes(EndNodes(i))=1;
    end
end

%%% Find distance from nodes to start node %%%
dist = sqrt((nodes(:,1)-RootNode(1)).^2+(nodes(:,2)-RootNode(2)).^2);
dist_tn = dist(EndNodes(find_tn==0));

%%% Set boundary values %%%
Bc_vals = [3*max(dist);zeros(Nn-1,1)];
for i = 2:length(Bc_vals)
    Bc_vals(i)=Bc_nodes(i)*dist(i);
end

%%% Solve linear system %%%
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]); % Forutsetter at kolonne 2 og 3 i «edges» angir node-nummer
connections_trimmed = connections(:,Bc_nodes==0);
A = [spdiags(g.^-1,0,Ne,Ne), connections_trimmed; connections_trimmed', sparse(size(connections_trimmed,2),size(connections_trimmed,2))];
rhs = [-connections*Bc_vals; zeros(sum(1-Bc_nodes),1)];

sol = A\rhs;
p = sol(Ne+1:end);
q = sol(1:Ne);
indexes = find(Bc_nodes==0);
p_all = Bc_vals;
p_all(indexes) = p;

% %%% PLOTS %%%%
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



%%% VORONOI DIAGRAM %%%
Ntn = length(Tn); % Antall terminalnoder
bs_ext = 2*[-2 2 -2 2];
[v,c,xy]=VoronoiLimit(Tn(:,1),Tn(:,2),'bs_ext',[-.8 .5 1.80 -.8;-.05 1.7 -.05 -.05]');
for i = 1:size(c,1)
    cell_index = find(xy==Tn(i,1:2));
    cells(i)=c(cell_index(1));
end
cells = cells';

%find(v(:,1)>bs_ext(2))
%find(v(:,2)>bs_ext(4))
% v(find(v(:,1)>bs_ext(2)),1)=bs_ext(2);
% v(find(v(:,1)<bs_ext(1)),1)=bs_ext(1);
% v(find(v(:,2)>bs_ext(4)),2)=bs_ext(4);
% v(find(v(:,2)<bs_ext(2)),2)=bs_ext(3);

figure('Name','Flux intensity plot')
map = gray(Ne);
q_sorted = sort(q,'descend');
for i = 1:Ntn
    coords = v(cells{i},:);
%     bigyval = find(coords(:,2) > bs_ext(4));
%     smallyval = find(coords(:,2) < bs_ext(3));
%     bigxval = find(coords(:,1) > bs_ext(2));
%     smallxval = find(coords(:,1) < bs_ext(1));
%     if any(bigyval)
%         coords(bigyval,2)=bs_ext(4);
%         v(c{cell_index(1)},:) = coords;
%     end
%     if any(smallyval)
%         coords(smallyval,2)=bs_ext(3);
%         v(c{cell_index(1)},:) = coords;
%     end
%     if any(bigxval)
%         coords(bigxval,1)=bs_ext(2);
%         v(c{cell_index(1)},:) = coords;
%     end
%     if any(smallxval)
%         coords(smallxval,1)=bs_ext(1);
%         v(c{cell_index(1)},:) = coords;
%     end
    pgon = polyshape(coords(:,1),coords(:,2));
    pg = plot(pgon);
    q_here = q(Te(i));
    ind = find(q_sorted==q_here);
    ind = ind(1);
    pg.FaceColor = map(ind,:);
    hold on
end
plot(Tn(:,1),Tn(:,2),'b.')
hold on
plot(RootNode(1),RootNode(2),'r*');
hold on
axis(gca,'equal')
%axis([-5 5 -5 5])


%%% Darcy flow %%%
p_cell = TPFA(cells,v)
IntensityMap(cells,v,p_cell)




                
                
