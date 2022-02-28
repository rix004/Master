clc; clear;
DT.Levels = 6;
DT.StartPos = [0 0];
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.1;    % Radius of trunk (mm)
DT.RadiusRate = 0.5;
DT.TrunkLength = 1;      % Length of trunk (mm)
DT.LengthRate = 0.7;

DetTree= GetTree(DT);
% RandomTree.nodes = DetTree.nodes;
% RandomTree.edges = DetTree.edges;
% RandomTree.TrunkRadius = 2;
% RandomTree.RadiusRate = 0.5;
% RandomTree = RRT_Tree(RandomTree,-2,2,-2,2,100);
% nodes = RandomTree.nodes; edges = RandomTree.edges;

% CALCULATE PRESSURE AND FLUX
Ne = length(edges); Nn = length(nodes);
mu = 5.5; % viscosity (mPa*s)

% Conductance
g(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));

Bc_nodes = [1; zeros(Nn/2-1,1); ones(Nn/2,1)];  % Noder som du gir randverdier i
Bc_vals = [1; zeros(Nn/2-1,1); zeros(Nn/2,1)];  % Randverdier

connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]); % Forutsetter at kolonne 2 og 3 i �edges� angir node-nummer
connections_trimmed = connections(:,Bc_nodes==0);
A = [spdiags(g.^-1,0,Ne,Ne), connections_trimmed; connections_trimmed', sparse(Nn/2-1,Nn/2-1)];
rhs = [-connections*Bc_vals; zeros(sum(1-Bc_nodes),1)];

sol = A\rhs;
p = sol(Ne+1:end);
q = sol(1:Ne);
p_all = Bc_vals;
p_all(2:(length(p)+1)) = p;

%%% PRESSURE PLOT %%%%
figure(2)
for i = 1:Ne
    x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
    y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
    pressure = [p_all(edges(i,2)) p_all(edges(i,3))];
    plot3(x,y,pressure,'.-','Color','Red','LineWidth',edges(i,4)*20)
    hold on
    grid on
end


%%% VORONOI-DIAGRAM %%%
StartNodes = edges(:,2);
EndNodes = edges(:,3);
Tn = []; Te = [];          % Terminalnoder og terminalkanter
for i = 1:length(EndNodes)
    v = ismember(EndNodes,StartNodes);
    if v(i) == 0
        Tn = [Tn;nodes(EndNodes(i),:)];
        Te = [Te;i];
    end
end
Ntn = length(Tn); % Antall terminalnoder
Nte = length(Te);
Tn_x = Tn(:,1);
Tn_y = Tn(:,2);
%p = Tn(:,1:2);
[v,c,xy]=VoronoiLimit(Tn_x,Tn_y,'bs_ext',[-2 2 -2 2;-2 -2 2 2]);

figure(3)
map = winter(Ntn);

for i = 1:Ntn
    coords = v(c{i},:);
    pgon = polyshape(coords(:,1),coords(:,2));
    pg = plot(pgon);
    q_here = q(Nte+i-1);
    q_max = max(q(Nte+i-1:end));
    q_min = min(q(Nte+i-1:end));
    ind = ceil(Ntn * q_here/q_max);
    pg.FaceColor = map(ind,:);
    hold on
end
plot(Tn_x,Tn_y,'r.')
axis(gca,'equal')
%axis([min(Tn(:,1))-abs(Tn(end,1)-Tn(end-1,1)) max(Tn(:,1))+abs(Tn(end,1)-Tn(end-1,1)) min(Tn(:,2))-abs(Tn(end,1)-Tn(end-1,1)) max(Tn(:,2))+abs(Tn(end,1)-Tn(end-1,1))]);
