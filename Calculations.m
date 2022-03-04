clc; clear;
addpath(genpath('C:///Users/jennyhognestad/Downloads/github_repo/'));
savepath;
DT.Levels = 7;
DT.StartPos = [0 0];
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.1;    
DT.RadiusRate = 0.5;
DT.TrunkLength = 1;
DT.LengthRate = 0.7;

DetTree= GetTree(DT);
%DrawTree(DetTree)


RandomTree.nodes = DetTree.nodes;
RandomTree.edges = DetTree.edges;
RandomTree.TrunkRadius = 2;
RandomTree.RadiusRate = 0.5;

RandomTree = RRT_Tree(RandomTree,-2,2,-2,2,50);
DrawTree(RandomTree)

nodes = RandomTree.nodes; edges = RandomTree.edges;
%nodes = DetTree.nodes; edges = DetTree.edges;

% CALCULATE PRESSURE AND FLUX
Ne = length(edges); Nn = length(nodes);
mu = 5.5; % viscosity (mPa*s)

% Conductance
g(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));

StartNodes = edges(:,2);
EndNodes = edges(:,3);
% Finner terminalnoder
Tn = [];
Te = [];
Bc_nodes = [1;zeros(Nn-1,1)];
for i = 1:length(EndNodes)
    v = ismember(EndNodes,StartNodes);
    if v(i) == 0
        Tn = [Tn;nodes(EndNodes(i),:)];
        Te = [Te;i];
        Bc_nodes(i+1)=1;
    end
end
Bc_vals = [1;zeros(Nn-1,1)];


%Bc_nodes = [1; zeros(Nn/2-1,1); ones(Nn/2,1)];  % Noder som du gir randverdier i
%Bc_vals = [1; zeros(Nn/2-1,1); zeros(Nn/2,1)];  % Randverdier


connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]); % Forutsetter at kolonne 2 og 3 i �edges� angir node-nummer
connections_trimmed = connections(:,Bc_nodes==0);
A = [spdiags(g.^-1,0,Ne,Ne), connections_trimmed; connections_trimmed', sparse(size(connections_trimmed,2),size(connections_trimmed,2))];
rhs = [-connections*Bc_vals; zeros(sum(1-Bc_nodes),1)];

sol = A\rhs;
p = sol(Ne+1:end);
q = sol(1:Ne);
p_all = Bc_vals;
p_all(2:(length(p)+1)) = p;

%%% PRESSURE PLOT %%%%
figure
for i = 1:Ne
    x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
    y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
    pressure = [p_all(edges(i,2)) p_all(edges(i,3))];
    plot3(x,y,pressure,'.-','Color','Red','LineWidth',edges(i,4)*20)
    hold on
    grid on
end


%%% VORONOI-DIAGRAM %%%
Ntn = length(Tn); % Antall terminalnoder
Tn_x = Tn(:,1);
Tn_y = Tn(:,2);
bs_ext = [-5 5 5 -5;-5 -5 5 5]';
%p = Tn(:,1:2);
[v,c,xy]=VoronoiLimit(Tn_x,Tn_y);
[v1,c1]=voronoin([Tn_x Tn_y]);
figure
map = winter(Ntn);

for i = 1:Ntn
    coords = v(c{i},:);
    pgon = polyshape(coords(:,1),coords(:,2));
    pg = plot(pgon);
    q_here = q(Te(i));
    q_max = max(q(Te(:)));
    q_min = min(q(Te(:)));
    ind = ceil(Ntn * q_here/q_max);
    if ind == 0 || isnan(ind)
        pg.FaceColor = map(1,:);
        hold on
    else
        pg.FaceColor = map(ind,:);
        hold on
    end
end
plot(Tn_x,Tn_y,'r.')
axis(gca,'equal')
%axis([min(Tn(:,1))-abs(Tn(end,1)-Tn(end-1,1)) max(Tn(:,1))+abs(Tn(end,1)-Tn(end-1,1)) min(Tn(:,2))-abs(Tn(end,1)-Tn(end-1,1)) max(Tn(:,2))+abs(Tn(end,1)-Tn(end-1,1))]);
