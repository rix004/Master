%cells are sorted by starting in the lower left corner (origo) counting
%first horizontally, then one layer at a time upwards.

%Edges are are sorted by first counting all horizontally aligned edges
%starting from the bottom left counting horizontally, then one layer at a time
%upwards. Then vertically alined edges are counted starting from the
%leftmost bottom edge counting upwards, then one layer at a time towards
%right. 
function [A, b, G,D,K, cells,edges] = assembleMatrices(num_cells_x,num_cells_y, f, k)

%Define mesh
dx = 1/(num_cells_x-1);
dy = 1/(num_cells_y-1);
num_cells = num_cells_x*num_cells_y;
cells = zeros(num_cells,6);

%Relevant edges in the dual grid
num_edges_vertical = (num_cells_x-1)*num_cells_y;
num_edges_horizontal = (num_cells_y-1)*num_cells_x;

num_edges = num_edges_vertical + num_edges_horizontal;
edges = zeros(num_edges,4);

%Set node values
for j = 1:num_cells_y
    for i = 1:num_cells_x
        cells(num_cells_x*(j-1)+i,1:2) = [dx*(i-1),dy*(j-1)];
    end
end

%Set edges connected to interior cells
for j = 2:num_cells_y-1
    for i = 2:num_cells_x-1
        cells(num_cells_x*(j-1)+i,3:6) = [num_edges_vertical+(i-1)*(num_cells_y-1)+(j-1), (j-1)*(num_cells_x-1)+i, num_edges_vertical + (i-1)*(num_cells_y-1)+ j, (j-1)*(num_cells_x-1)+i-1];
    end
end

%Set edges connected to lower/top cells
for i = 2:num_cells_x-1
    cells(i,4:6) = [i,num_edges_vertical+(i-1)*(num_cells_y-1)+1,i-1]; 
    cells(i+(num_cells_y-1)*num_cells_x,[3,4,6]) = [num_edges_vertical+i*(num_cells_y-1), (num_cells_y-1)*(num_cells_x-1)+i, (num_cells_y-1)*(num_cells_x-1)+i-1];
end

%Set edges connected to left/right cells
for i = 2:num_cells_y-1
    cells((i-1)*num_cells_x+1,3:5) = [num_edges_vertical+i-1,(i-1)*(num_cells_x-1)+1 , num_edges_vertical+i]; 
    cells(i*num_cells_x,[3,5,6]) = [num_edges_vertical+(num_cells_x-1)*(num_cells_y-1)+i-1,num_edges_vertical+(num_cells_x-1)*(num_cells_y-1)+i, i*(num_cells_x-1)];
end

%Set edges connected to boundary cells
cells(1,4:5) = [1,num_edges_vertical+1];
cells(num_cells_x, 5:6) = [num_edges_vertical+(num_cells_x-1)*(num_cells_y-1)+1, num_cells_x-1];
cells((num_cells_y-1)*num_cells_x+1,3:4) = [num_edges_vertical+num_cells_y-1, (num_cells_y-1)*(num_cells_x-1)+1];
cells(num_cells_x*num_cells_y,[3,6]) = [num_edges, num_edges_vertical];

%Set edge values + Connectivity with cells
for j = 1:num_cells_y
    for i = 1:(num_cells_x-1)
        edges((num_cells_x-1)*(j-1)+i,:) = [0.5*dx+(i-1)*dx,(j-1)*dy,i+(num_cells_x)*(j-1), i+(num_cells_x)*(j-1)+1];
    end
end

for j = 1:num_cells_x
    for i = 1:(num_cells_y-1)
        edges(num_edges_vertical+(num_cells_y-1)*(j-1)+i,:) = [(j-1)*dx,0.5*dy+(i-1)*dy,(j-1)+1 + num_cells_x*(i-1),(j-1)+1 + num_cells_x*i];
    end
end


%Grad operator:
grad_occupation_i = zeros(1,2*num_edges);
grad_occupation_j = grad_occupation_i;
grad_values = grad_occupation_i;
diffusion_values = grad_occupation_i;
for i = (num_cells_x):num_edges_vertical-(num_cells_x-1)
    grad_occupation_i((2*i-1):2*i) = i;
    grad_occupation_j((2*i-1):2*i) = edges(i,3:4);
    grad_values((2*i-1):2*i) = [-1,1]/dx; 
    diffusion_values((2*i-1):2*i) = 2*[1,1]/(1/k(cells(edges(i,3),1),cells(edges(i,3),2))+1/k(cells(edges(i,4),1),cells(edges(i,4),2))); 
end
%Remove impact from cells on the boundary LEFT/RIGHT
for i = 1:num_cells_y-2
    grad_occupation_i(2*(i*(num_cells_x-1)+1)-1) = 0;
    grad_occupation_j(2*(i*(num_cells_x-1)+1)-1) = 0;
    grad_values(2*(i*(num_cells_x-1)+1)-1) = 0;
    diffusion_values(2*(i*(num_cells_x-1)+1)-1:2*(i*(num_cells_x-1)+1)) = [0,1]*k(cells(edges(i*(num_cells_x-1)+1,4),1),cells(edges(i*(num_cells_x-1)+1,4),2));

    grad_occupation_i(2*((i+1)*(num_cells_x-1))) = 0;
    grad_occupation_j(2*((i+1)*(num_cells_x-1))) = 0;
    grad_values(2*((i+1)*(num_cells_x-1))) = 0;
    diffusion_values(2*((i+1)*(num_cells_x-1))-1:2*((i+1)*(num_cells_x-1))) = [1,0]*k(cells(edges((i+1)*(num_cells_x-1),3),1),cells(edges((i+1)*(num_cells_x-1),3),2));
end
for i = (num_cells_y):(num_edges_horizontal-(num_cells_y-1))
    grad_occupation_i(2*num_edges_vertical+(2*i-1):2*num_edges_vertical+2*i) = num_edges_vertical+i;
    grad_occupation_j(2*num_edges_vertical+(2*i-1):2*num_edges_vertical+2*i) = edges(num_edges_vertical+i,3:4);
    grad_values(2*num_edges_vertical+(2*i-1):2*num_edges_vertical+2*i) = [-1,1]/dy; 
    diffusion_values(2*num_edges_vertical+(2*i-1):2*num_edges_vertical+2*i) = 2*[1,1]/(1/k(cells(edges(num_edges_vertical+i,3),1),cells(edges(num_edges_vertical+i,3),2))+1/k(cells(edges(num_edges_vertical+i,4),1),cells(edges(num_edges_vertical+i,4),2))); 
end
%Remove impact from cells on the bottom/top boundary
for i = 1:num_cells_x-2
    grad_occupation_i(2*(num_edges_vertical+i*(num_cells_y-1)+1)-1) = 0;
    grad_occupation_j(2*(num_edges_vertical+i*(num_cells_y-1)+1)-1) = 0;
    grad_values(2*(num_edges_vertical+i*(num_cells_y-1)+1)-1) = 0;
    diffusion_values(2*(num_edges_vertical+i*(num_cells_y-1)+1)-1:2*(num_edges_vertical+i*(num_cells_y-1)+1)) = [0,1]*k(cells(edges(num_edges_vertical+i*(num_cells_y-1)+1,4),1),cells(edges(num_edges_vertical+i*(num_cells_y-1)+1,4),2));

    grad_occupation_i(2*(num_edges_vertical+(i+1)*(num_cells_y-1))) = 0;
    grad_occupation_j(2*(num_edges_vertical+(i+1)*(num_cells_y-1))) = 0;
    grad_values(2*(num_edges_vertical+(i+1)*(num_cells_y-1))) = 0;
    diffusion_values(2*(num_edges_vertical+(i+1)*(num_cells_y-1))-1:2*(num_edges_vertical+(i+1)*(num_cells_y-1))) = [1,0]*k(cells(edges(num_edges_vertical+(i+1)*(num_cells_y-1),3),1),cells(edges(num_edges_vertical+(i+1)*(num_cells_y-1),3),2));
end
grad_occupation_i(grad_occupation_i ==0)= [];
grad_occupation_j(grad_occupation_j ==0)= [];
grad_values(grad_values ==0)= [];
diffusion_values(diffusion_values ==0)= [];
G = sparse(grad_occupation_i,grad_occupation_j,grad_values,num_edges,num_cells);
K = sparse(grad_occupation_i,grad_occupation_j,diffusion_values,num_edges,num_cells);


%Div operator
div_occupation_i = zeros(1,4*num_cells);
div_occupation_j = div_occupation_i;
div_values = div_occupation_i;

%interior cells
for j = 2:num_cells_y-1
    for i = 2:num_cells_x-1
        index = num_cells_x*(j-1)+i;
        div_occupation_i(4*index-3:4*index) = index;
        div_occupation_j(4*index-3:4*index) = cells(index,3:6);
        div_values(4*index-3:4*index) = [-1/dy, 1/dx, 1/dy, -1/dx];
    end
end


div_occupation_i(div_occupation_i ==0)= [];
div_occupation_j(div_occupation_j ==0)= [];
div_values(div_values ==0)= [];


D = sparse(div_occupation_i, div_occupation_j,div_values,num_cells, num_edges);

A = -D*(K.*G);

%Create RHS
b = ones(num_cells,1);
for i=1:num_cells
    b(i) = f(cells(i,1),cells(i,2));
end

%Hack boundary cells
boundary_cells = zeros(1,2*num_cells_x+2*(num_cells_y-2));
for i=1:num_cells_x
    boundary_cells(i)= i;
    boundary_cells(num_cells_x+i) = (num_cells_y-1)*num_cells_x+i;
end
for i=2:num_cells_y-1
    boundary_cells(2*num_cells_x +(i-1)) = (i-1)*num_cells_x+1;
    boundary_cells(2*num_cells_x +num_cells_y-2+(i-1)) = i*num_cells_x; 
end
for i=boundary_cells
    A(i,i) = 1;
    %b(i) = 0;
end