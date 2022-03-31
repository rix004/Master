function [nodes,cells] = GridGeneration(num_nodes_x,num_nodes_y,Domain)
 num_nodes = num_nodes_x*num_nodes_y;
 num_elements = (num_nodes_x-1)*(num_nodes_y-1);
 
 % Node matrix - position of nodes 
 dx = (Domain(2)-Domain(1))/(num_nodes_x - 1);
 dy = (Domain(4)-Domain(3))/(num_nodes_y - 1);
 nodes = zeros(num_nodes,2);
 for i = 1:num_nodes_y
     for j = 1:num_nodes_x
        nodes(num_nodes_x*(i-1)+j,1:2)=[Domain(1)+dx*(j-1),Domain(3)+dy*(i-1)];
     end
 end
 
 % Element matrix - connections of nodes into triangles
 k = 0;
 cells = cell(num_elements,1);
 %elements = zeros(num_elements,4);
 for i = 1:num_nodes_y-1
     for j = 1:num_nodes_x-1
        node_p = num_nodes_x*(i-1)+j;
        cells{k+1}=[node_p,node_p+1,num_nodes_x+node_p+1,num_nodes_x+node_p];
        k=k+1;
     end
 end
 end