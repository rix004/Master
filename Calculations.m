clc; clear;
RootNode = [0 0];

% Deterministic tree data
DT.Levels = 4;
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
Np = 500;
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
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]); % Forutsetter at kolonne 2 og 3 i �edges� angir node-nummer
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
bs_ext = [-5 5 5 -5;-5 -5 5 5]';
[v,c,xy]=VoronoiLimit(Tn(:,1),Tn(:,2));

figure('Name','Flux intensity plot')
map = gray(Ne);
for i = 1:Ntn
    cell_index = find(xy==Tn(i,1:2));
    coords = v(c{cell_index(1)},:);
    pgon = polyshape(coords(:,1),coords(:,2));
    pg = plot(pgon);
    q_here = q(Te(i));
    q_sorted = sort(q,'descend');
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
axis(Domain)



%%% Darcy flow %%%
cells = zeros(Ntn,2);
k = 1; % Permeability
m = 1;
for i = 1:size(c,1)
    cell_index = find(xy==Tn(i,1:2));
    coords = v(c{cell_index(1)},:);
    cell_coords = v(c{cell_index(1)},:);
    cells(i,1:2)=[mean(cell_coords(:,1)) mean(cell_coords(:,2))];
    num_vertices = length(c{cell_index(1),:});
    for j = 1:num_vertices
        if j < num_vertices
            pedgesloc=c{cell_index(1)}(j:j+1);
        elseif j == length(c{cell_index(1),:})
            pedgesloc=[c{cell_index(1)}(j) c{cell_index(1)}(1)];
        end
        if i > 1
            existing_edge=ismember(pedges, pedgesloc, "rows")+ismember(flip(pedges,2), pedgesloc, "rows");
        else 
            existing_edge = 0;
        end
        if any(existing_edge)
            edge_number=find(existing_edge==1);
            p_con(i,edge_number)=-1;
        else
            pedges(m,:)=pedgesloc;
            edge_number = m;
            m = m+1;
            p_con(i,edge_number)=1;
        end
    end
                
end

sp_con = sparse(p_con);
p_cell = p_all(Te);
bc_edges = zeros(size(p_con,2),1);
x1 = v(pedges(:,1),1);
x2 = v(pedges(:,2),1);
y1 = v(pedges(:,1),2);
y2 = v(pedges(:,2),2);
l_edge = sqrt((x1-x2).^2+(y1-y2).^2);
d = zeros(size(p_con,2),1);
for i = 1:size(p_con,2)
    n_cells = find(p_con(:,i)~=0);    % neighbouring cells
    if length(n_cells)==2
        d(i) = sqrt((cells(n_cells(1),1)-cells(n_cells(2),1))^2+(cells(n_cells(1),2)-cells(n_cells(2),2))^2);
    elseif length(n_cells) == 1
            X1 = [v(pedges(i,1),1) v(pedges(i,1),2)];
            X2 = [v(pedges(i,2),1) v(pedges(i,2),2)];
            p = [cells(n_cells,1) cells(n_cells,2)];
            D1 = norm(p-X1);
            D2 = norm(X2-X1);
            D3 = norm(p-X2);
            theta = acosd((D1^2+D2^2-D3^2)/(2*D1*D2));
            d_half = D1*sind(theta);
            d(i) = d_half*2;
            bc_edges(i)=1;
    end
end

% LHS
T = k*l_edge./d;
con_trans = sp_con';
B = zeros(size(con_trans));
for i = 1:size(con_trans,1)
    if bc_edges(i)==0
        B(i,:)=con_trans(i,:)*T(i);
    elseif bc_edges(i)==1
    end
end







                
                