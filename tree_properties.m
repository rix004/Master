clc; clear;
mu = 5.5;           % viscosity (mPa*s)
levels = 9;
R_trunk = 0.1;    % Radius of trunk (mm)
L_trunk = 1;      % Length of trunk (mm)
x0 = 0;
y0 = 0;
z0 = 0;
start_angle = 90;
R_rate = 0.5;
L_rate = 0.7;
rotation_angle = 90;
[nodes, edges]= tree(levels,R_trunk,L_trunk,x0,y0,z0,start_angle,R_rate,L_rate,rotation_angle);
Ne = length(edges); Nn = length(nodes);

% Conductance
g = zeros(length(edges),1);
for i = 1:length(edges)
    g(i) = pi*edges(i,4)^4/(8*mu*edges(i,1));
end


Bc_nodes = [1; zeros(Nn/2-1,1); ones(Nn/2,1)];  % Noder som du gir randverdier i
Bc_vals = [1; zeros(Nn/2-1,1); zeros(Nn/2,1)];  % Randverdier

connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]); % Forutsetter at kolonne 2 og 3 i «edges» angir node-nummer
connections_trimmed = connections(:,Bc_nodes==0);
A = [spdiags(g.^-1,0,Ne,Ne), connections_trimmed; connections_trimmed', sparse(Nn/2-1,Nn/2-1)];
rhs = [-connections*Bc_vals; zeros(sum(1-Bc_nodes),1)];

sol = A\rhs;
p = sol(1:Nn);
% q = sol(1:Ne)

for i = 1:Ne
    q(i) = g(i)*(p(i+1)-p(i));
end

figure(2)
for i = 1:Ne-1
    x = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
    y = [nodes(edges(i,2),3) nodes(edges(i,3),3)];
    pressure = [p(edges(i,2)) p(edges(i,3))];
    plot3(x,y,pressure,'.-','Color','Red')
    hold on
    grid on
end
