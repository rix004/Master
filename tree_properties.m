clc; clear;
mu = 5.5;           % viscosity (mPa*s)
levels = 9;
R_trunk = 0.1;    % Radius of trunk (mm)
L_trunk = 1;      % Length of trunk (mm)
x0 = 0;
y0 = 0;
theta = 90;
R_rate = 0.5;
L_rate = 0.7;
rotation_angle = 90;
[nodes, edges]= GetMatrices(levels,R_trunk,L_trunk,x0,y0,theta,R_rate,L_rate,rotation_angle);
Ne = length(edges); Nn = length(nodes);


% Conductance
g(:,1) = pi*edges(:,4).^4./(8*mu*edges(:,1));

Bc_nodes = [1; zeros(Nn/2-1,1); ones(Nn/2,1)];  % Noder som du gir randverdier i
Bc_vals = [1; zeros(Nn/2-1,1); zeros(Nn/2,1)];  % Randverdier

connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]); % Forutsetter at kolonne 2 og 3 i «edges» angir node-nummer
connections_trimmed = connections(:,Bc_nodes==0);
A = [spdiags(g.^-1,0,Ne,Ne), connections_trimmed; connections_trimmed', sparse(Nn/2-1,Nn/2-1)];
rhs = [-connections*Bc_vals; zeros(sum(1-Bc_nodes),1)];

sol = A\rhs;
p = sol(Ne+1:end);
q = sol(1:Ne);
p_all = Bc_vals;
p_all(2:(length(p)+1)) = p;

figure(2)
for i = 1:Ne
    x = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
    y = [nodes(edges(i,2),3) nodes(edges(i,3),3)];
    pressure = [p_all(edges(i,2)) p_all(edges(i,3))];
    plot3(x,y,pressure,'.-','Color','Red','LineWidth',edges(i,4)*20)
    hold on
    grid on
end


%%% VORONOI-DIAGRAM %%%%
Ntn = 2^(levels-2); % Antall terminalnoder
terminal_nodes = nodes(Nn-Ntn+1:end,:);
tn_x = terminal_nodes(:,2);
tn_y = terminal_nodes(:,3);
p = terminal_nodes(:,2:3);
[v,c,xy]=VoronoiLimit(tn_x,tn_y,'bs_ext',[-2 2 -2 2;-2 -2 2 2]);

figure(3)
for i = 1:Ntn
    coords = v(c{i},:);
    pgon = polyshape(coords(:,1),coords(:,2));
    pg = plot(pgon);
    pg.FaceColor = i/200*[0.4 0.6 0.7];
    hold on
end
plot(tn_x,tn_y,'.')
axis([-2 2 -0.5 2]);
