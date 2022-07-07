function[]=SystemTestEx(nodes,edges,TNlogic,TNinfo,mu,k,p_darcy,q_darcy,p_network,q_network,Grad_D,boundary_cells,cell_area,D_bvs,K_N,MacroTermIndexes,MicroTermIndexes)
% Testing
Ne = length(edges);
Nn = length(nodes);
Nin = sum(1-TNlogic)-1;                        % Interior nodes
Npn = sum(1-TNlogic);                          % Parent nodes
Ntn = size(TNinfo,1);

% Check that the flux out of the domain (over dOmega) is the same as the
% flux coming into the domain and the flux into the network
disp('Flux out of domain, into domain and into network (should be the same):')
sum(Grad_D(boundary_cells(:,3),:)*p_darcy-D_bvs(boundary_cells(:,3)))
sum(q_darcy)
q_network(1)
sum(Grad_D(boundary_cells(:,3),:)*p_darcy-D_bvs(boundary_cells(:,3)))-q_network(1)
q_network(1)-sum(q_darcy)

% Check that the sum of fluxes in and out of an interior node equals 0.
int_nodes = find((TNlogic+[1;zeros(Nn-1,1)])==0);
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

flux_sum = zeros(Ntn,1);
for i = 1:Ntn
    edge_in=find(edges(:,3)==MicroTermIndexes(i));
    flux_in=q_network(edge_in);
    darcy_flux = q_darcy(i);
    flux_sum(i) = flux_in-darcy_flux;
end
disp('Flux into terminal nodes-flux into Darcy domain (should be zero):')
sum(flux_sum)

% Chech that the pressure in terminal nodes equals the pressure in the
% Darcy-cell + Peaceman
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
for i = 1:size(TNinfo,1)
    sum_PL(i) = K_N(TNinfo(i,3))*(p_network(edges(TNinfo(i,3),2))-p_network(edges(TNinfo(i,3),3)))-q_network(TNinfo(i,3));
end
sum(sum_PL)
end