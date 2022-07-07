% Two point flux approximation.
% Dirichlet BC: BC_type = 1.
% Neumann BC: BC_type = 2.
function[A,LHS,bound_vals,RHS,cell_center,pedges,cell_area,boundary_points,bv_out] = TPFA(cells,vertices,f,K_D,BC_type,bc,r_0)
cell_center = zeros(size(cells,1),2);
m = 1;

% Find all edges
for i = 1:size(cells,1)
    coords = vertices(cells{i},:);
    cell_area(i) = polyarea(coords(:,1),coords(:,2));
    cell_center(i,1:2)=[mean(coords(:,1)) mean(coords(:,2))];
    num_vertices = length(cells{i,:});
    for j = 1:num_vertices
        if j < num_vertices
            pedgesloc=[cells{i}(j) cells{i}(j+1)];
        elseif j == length(cells{i,:})
            pedgesloc=[cells{i}(j) cells{i}(1)];
        end
        if i > 1
            existing_edge=ismember(pedges, pedgesloc, "rows")+ismember(flip(pedges,2), pedgesloc, "rows");
        else 
            existing_edge = 0;
        end
        if any(existing_edge)
            edge_number=existing_edge==1;
            p_con(i,edge_number)=-1;
        else
            pedges(m,:)=pedgesloc;
            edge_number = m;
            m = m+1;
            p_con(i,edge_number)=1;
        end
    end
                
end

% Find length of edges
num_edges = size(p_con,2);
bc_edges = zeros(num_edges,1);
x1 = vertices(pedges(:,1),1);
x2 = vertices(pedges(:,2),1);
y1 = vertices(pedges(:,1),2);
y2 = vertices(pedges(:,2),2);
l_edge = sqrt((x1-x2).^2+(y1-y2).^2);
d = zeros(num_edges,1);

% Make connectivity matrix
Div = sparse(p_con);
Grad = Div';
bound_count = 1;
for i = 1:num_edges
    n_cells = find(p_con(:,i)~=0);    % neighbouring cells
    if length(n_cells)==2
        d(i) = sqrt((cell_center(n_cells(1),1)-cell_center(n_cells(2),1))^2+(cell_center(n_cells(1),2)-cell_center(n_cells(2),2))^2);
    elseif length(n_cells) == 1
            X1 = [vertices(pedges(i,1),1) vertices(pedges(i,1),2)];
            X2 = [vertices(pedges(i,2),1) vertices(pedges(i,2),2)];
            p = [cell_center(n_cells,1) cell_center(n_cells,2)];
            d(i) = DistanceToEdge(p,X1,X2);
            bc_edges(i)=1;
            boundary_points(bound_count,1:3)=[mean([X1(1),X2(1)]) mean([X1(2),X2(2)]) i];
            bound_count = bound_count+1;
    end
end

bc_index = find(bc_edges==1);
boundary_cells = zeros(length(bc_index),1);
for i = 1:length(bc_index)
    boundary_cells(i)=find(Div(:,bc_index(i))==1);
end
boundary_cells=unique(boundary_cells);

% Construct matrices
B = l_edge./d;
flux_boundary = zeros(num_edges);
for i = 1:num_edges
    A(i,:)=Grad(i,:)*B(i);
    if bc_edges(i) == 1 && BC_type == 2
        flux_boundary(i,i)=1;
    elseif bc_edges(i) == 1 && BC_type == 1
        flux_boundary(i,i)=B(i)*K_D; % Må endres dersom K_D skal avhenge av posisjon
    end
end

for i = 1:size(cells,1)
    A(:,i)=K_D*A(:,i);
end
LHS = Div*A;

% RHS
% Source term
for i = 1:size(cell_center,1)
    b(i) = f(cell_center(i,1),cell_center(i,2))*cell_area(i);
end
b = b';

% Boundary values
bv = zeros(num_edges,1);
bv(bc_edges==1)=bc;
% Test with Peaceman correction
% r = sqrt((0-boundary_points(:,1)).^2+(0-boundary_points(:,2)).^2);
% bv(bc_edges==1)=-bc/(2*pi)*log(r./r_0);
bound_vals = flux_boundary*bv;
RHS = Div*flux_boundary*bv + b;
bv_out = bv(bc_edges==1);
end

