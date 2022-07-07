
%%% Testing %%%

Ne = length(edges);
Nn = length(nodes);
Nin = sum(1-TermNodes)-1;                        % Interior nodes
Npn = sum(1-TermNodes);                          % Parent nodes
Ntn = size(MicroTN,1);

% Check that the flux out of the domain (over dOmega) is the same as the
% flux coming into the domain and the flux into the network
disp('Flux out of domain, into domain and into network (should be the same):')
sum(Grad_D(boundary_cells(:,3),:)*p_darcy-D_bvs(boundary_cells(:,3)))
sum(q_darcy)
q_network(1)

% Check that the sum of fluxes in and out of an interior node equals 0.
int_nodes = find((TermNodes+[1;zeros(Nn-1,1)])==0);
total_flux = zeros(Nin,1);
for i =1:length(int_nodes)
    edge_in=find(edges(:,3)==int_nodes(i));
    edge_out=find(edges(:,2)==int_nodes(i));
    flux_in = q_network(edge_in);
    flux_out = q_network(edge_out);
    total_flux(i) = sum(flux_in)-sum(flux_out);
end
disp('Sum flux in and out of interior nodes (should be zero):')
sum(total_flux)
% % Check that terminal node corresponds to darcy cell
% term_nodes = find(TermNodes==1);
% for i = 1:Ntn
%     coords = vertices(cells{i},:);
%     PointInside(i) = inpolygon(nodes(TermNodes(i),1),nodes(TermNodes(i),2),coords(:,1),coords(:,2));
% end
% 
% % Check that the flux into terminal nodes equals the flux into the Darcy domain
flux_sum = zeros(Ntn,1);
for i = 1:Ntn
    edge_in=find(edges(:,3)==MicroTermIndexes(i));
    flux_in=q_network(edge_in);
    darcy_flux = q_darcy(i);
    flux_sum(i) = flux_in-darcy_flux;
end
disp('Flux into terminal nodes-flux into Darcy domain (should be zero):')
sum(flux_sum)

% % Chech that the pressure in terminal nodes equals the pressure in the
% % Darcy-cell + Peaceman
p_sum = zeros(Ntn,1);
for i = 1:Ntn
    p_term = p_network(MicroTermIndexes(i));
    edge_in=find(edges(:,3)==MicroTermIndexes(i));
    flux_in=q_network(edge_in);
    p_peaceman = flux_in*mu/(2*pi*k)*log(0.2*sqrt(cell_area(i))/edges(edge_in,4));
    p_sum(i) = p_darcy(i)-p_term+p_peaceman;
end
disp('p_darcy-p_term+p_peaceman (should be zero):')
sum(p_sum)
                
%% Check that Poiseuilles law is also valid in terminal edges
disp('Poiseuilles law in terminal edge (should be zero)')
for i = 1:size(MicroTN,1)
    sum_PL(i) = K_N(MicroTN(i,3))*(p_network(edges(MicroTN(i,3),2))-p_network(edges(MicroTN(i,3),3)))-q_network(MicroTN(i,3));
end
sum(sum_PL)

%%% Test med fundamentalløsningen %%%
p_exact = @(x,y) -Neu_network*mu/(2*pi*k)*log(sqrt(x.^2+y.^2)/edges(MicroTN(i,3),4));
epsilon = 0.1;
distToOrigin = sqrt(cell_center(:,1).^2+(cell_center(:,2).^2));
ok_cells = find(distToOrigin>epsilon);
notok_cell = find(distToOrigin<epsilon);
error = p_darcy(ok_cells)-p_exact(cell_center(ok_cells,1),cell_center(ok_cells,2));
l2_error = 0;
for i =1:length(error)
    l2_error = l2_error + error(i)^2*cell_area(i);
end
sqrt(l2_error);
h=sqrt(4/size(cells,1));


figure('Name','Pressure plot')
for i = 1:Ne
    x = [nodes(edges(i,2),1) nodes(edges(i,3),1)];
    y = [nodes(edges(i,2),2) nodes(edges(i,3),2)];
    pressure = [p_network(edges(i,2)) p_network(edges(i,3))];
    plot3(x,y,pressure,'.-','Color',[0 0 0],'LineWidth',edges(i,4)*10)
    hold on
    grid on
end
p_map = autumn(Nn+Ntn*3);
network_map = p_map(Ntn*2:end,:);
for i = 1:Nn
    x = nodes(i,1);
    y = nodes(i,2);
    pressure = p_network(i);
    p_sorted = sort(p_network);
    ind = find(p_sorted==pressure);
    ind = ind(1);
    plot3(x,y,pressure,'.','Color',network_map(ind,:),'MarkerSize',20);
    hold on
end
zlabel('Node pressure')
darcy_map = p_map(1:3*Ntn,:);

[xq,yq]=meshgrid(D(1):0.05:D(2), D(3):0.05:D(4));
p_points = [cell_center(:,1:2);boundary_cells(:,1:2);[-1 -1;-1 1;1 1;1 -1]];
p_values = [p_darcy;bv;Bv_darcy*ones(4,1)];
pEx_values = [p_exact(cell_center(ok_cells,1),cell_center(ok_cells,2));bv];
vq = griddata(p_points(:,1),p_points(:,2),p_values,xq,yq);
surf(xq,yq,vq)
colormap(darcy_map)