
%%% Testing %%%

% Check that the flux out of the domain (over dOmega) is the same as the
% flux coming into the domain and the flux into the network
disp('Flux out of domain, into domain and into network (should be the same):')
sum(Grad_D(boundary_cells(:,3),:)*p_darcy-D_bvs(boundary_cells(:,3)))
sum(q_darcy)
q_network(1)

% Check that the sum of fluxes in and out of an interior node equals 0.
int_nodes = find((Term_nodes+[1;zeros(Nn-1,1)])==0);
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
term_nodes = find(Term_nodes==1);
for i = 1:Ntn
    coords = vertices(cells{i},:);
    PointInside(i) = inpolygon(nodes(term_nodes(i),1),nodes(term_nodes(i),2),coords(:,1),coords(:,2));
end
% 
% % Check that the flux into terminal nodes equals the flux into the Darcy domain
flux_sum = zeros(Ntn,1);
for i = 1:Ntn
    edge_in=find(edges(:,3)==term_nodes(i));
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
    p_term = p_network(term_nodes(i));
    edge_in=find(edges(:,3)==term_nodes(i));
    flux_in=q_network(edge_in);
    p_peaceman = flux_in*mu/(2*pi*k)*log(0.2*sqrt(cell_area(i))/edges(edge_in,4));
    p_sum(i) = p_darcy(i)-p_term+p_peaceman;
end
disp('p_darcy-p_term+p_peaceman (should be zero):')
sum(p_sum)
                
%% Check that Poiseuilles law is also valid in terminal edges
disp('Poiseuilles law in terminal edge (should be zero)')
for i = 1:length(Te)
    sum_PL(i) = K_N(Te(i))*(p_network(edges(Te(i),2))-p_network(edges(Te(i),3)))-q_network(Te(i));
end
sum(sum_PL)

                