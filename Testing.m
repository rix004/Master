
%%% Testing %%%

% Check that the sum of fluxes in and out of an interior node equals 0.
% int_nodes = find(Bc_nodes==0);
% total_flux = zeros(Nin,1);
% for i =1:length(int_nodes)
%     edge_in=find(edges(:,3)==int_nodes(i));
%     edge_out=find(edges(:,2)==int_nodes(i));
%     flux_in = q_network(edge_in);
%     flux_out = q_network(edge_out);
%     total_flux(i) = sum(flux_in)-sum(flux_out);
% end
% 
% % Check that terminal node corresponds to darcy cell
% term_nodes = find([0;Bc_nodes(2:end)]==1);
% for i = 1:Ntn
%     coords = vertices(cells{i},:);
%     PointInside(i) = inpolygon(nodes(term_nodes(i),1),nodes(term_nodes(i),2),coords(:,1),coords(:,2));
% end
% 
% % Check that the flux into terminal nodes equals the flux into the Darcy domain
% sum = zeros(Ntn,1);
% for i = 1:Ntn
%     edge_in=find(edges(:,3)==term_nodes(i));
%     flux_in=q_network(edge_in);
%     darcy_flux = darcy_source(i);
%     sum(i) = flux_in-darcy_flux;
% end
% 
% 
% % Chech that the pressure in terminal nodes equals the pressure in the
% % Darcy-cell + Peaceman
% p_sum = zeros(Ntn,1);
% for i = 1:Ntn
%     p_term = p_all(term_nodes(i));
%     edge_in=find(edges(:,3)==term_nodes(i));
%     flux_in=q_network(edge_in);
%     p_peaceman = flux_in*mu/(2*pi*k)*log(0.2*sqrt(cell_area(i))/edges(edge_in,4));
%     p_sum(i) = p_darcy(i)-p_term+p_peaceman;
% end


                