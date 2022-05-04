function [nodes,cells,triangles] = GridGeneration(num_nodes_x,num_nodes_y,Domain)
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
 
 % Cell matrix - connections of nodes into rectangles
 % Triangle matrix - connections of nodes into triangles
 k1 = 0;
 k2 = 0;
 cells = cell(num_elements,1);
 triangles = cell(num_elements,1);
 for i = 1:num_nodes_y-1
     for j = 1:num_nodes_x-1
        node_p = num_nodes_x*(i-1)+j;
        cells{k1+1}=[node_p,node_p+1,num_nodes_x+node_p+1,num_nodes_x+node_p];
        k1=k1+1;
        triangles{k2+1}=[node_p,num_nodes_x+node_p+1,num_nodes_x+node_p];
        triangles{k2+2}=[node_p,node_p+1,num_nodes_x+node_p+1];
        k2=k2+2;
     end
 end
 end