clc; clear;
mu = 5.5;           % viscosity (mPa*s)
levels = 9;
R_trunk = 0.1;    % Radius of trunk (mm)
L_trunk = 1;      % Length of trunk (mm)
x0 = 0;
y0 = 0;
z0 = 0;
theta = 90;
phi = 0;
R_rate = 0.5;
L_rate = 0.7;
delta_theta = 90;
delta_phi = 30;

[nodes, edges]= tree3D(levels,R_trunk,L_trunk,x0,y0,z0,theta,phi,R_rate,L_rate,delta_theta,delta_phi);

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
p_inner = sol(Ne+1:end);
q = sol(1:Ne);
p = Bc_vals;
p(2:2^(levels-2)) = p_inner;

figure(2)
for i = 1:Ne
    x = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
    y = [nodes(edges(i,2),3) nodes(edges(i,3),3)];
    pressure = [p(edges(i,2)) p(edges(i,3))];
    plot3(x,y,pressure,'.-','Color','Red')
    hold on
    grid on
end
