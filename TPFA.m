function[p] = TPFA(cells,vertices)
cell_center = zeros(size(cells,1),2);
k = 1; % Permeability
m = 1;
for i = 1:size(cells,1)
    coords = vertices(cells{i},:);
    cell_center(i,1:2)=[mean(coords(:,1)) mean(coords(:,2))];
    num_vertices = length(cells{i,:});
    for j = 1:num_vertices
        if j < num_vertices
            pedgesloc=cells{i}(j:j+1);
        elseif j == length(cells{i,:})
            pedgesloc=[cells{i}(j) cells{i}(1)];
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
num_edges = size(p_con,2);

div = sparse(p_con);
bc_edges = zeros(num_edges,1);
x1 = vertices(pedges(:,1),1);
x2 = vertices(pedges(:,2),1);
y1 = vertices(pedges(:,1),2);
y2 = vertices(pedges(:,2),2);
l_edge = sqrt((x1-x2).^2+(y1-y2).^2);
d = zeros(num_edges,1);
for i = 1:num_edges
    n_cells = find(p_con(:,i)~=0);    % neighbouring cells
    if length(n_cells)==2
        d(i) = sqrt((cell_center(n_cells(1),1)-cell_center(n_cells(2),1))^2+(cell_center(n_cells(1),2)-cell_center(n_cells(2),2))^2);
    elseif length(n_cells) == 1
            X1 = [vertices(pedges(i,1),1) vertices(pedges(i,1),2)];
            X2 = [vertices(pedges(i,2),1) vertices(pedges(i,2),2)];
            p = [cell_center(n_cells,1) cell_center(n_cells,2)];
            d_half = DistanceToEdge(p,X1,X2);
            d(i) = d_half*2;
            bc_edges(i)=1;
    end
end
pedges
bc_edges
% Construct matrices
T = k*l_edge./d;
flux_boundary = zeros(num_edges);
for i = 1:num_edges
    B(:,i)=div(:,i)*T(i);
    if bc_edges(i) == 1
        flux_boundary(i,i)=T(i);
    end
end
B = B';

% Boundary values
bv = zeros(num_edges,1);
bv(4)=1;
bv(8)=1;
LHS = div*B
RHS = div*flux_boundary*bv;
p = LHS\RHS;
end