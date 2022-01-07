clc;clear;
format short
axis(gca,'equal')
distance = 0.7;
levels = 5;
num_nodes = 2^(levels-1);
num_edges = num_nodes-1;
x1 = 0;
y1 = 0;
z1 = 0;
degrees = 90;
root_radius = 5;
x_val = [x1];
y_val = [y1];
z_val = [z1];

allnodes = newbranch(x_val,y_val,[levels],x1,y1,degrees,distance,levels,root_radius)
%allnodes3D = newbranch3D(x_val,y_val,z_val,x1,y1,z1,degrees,distance,levels,root_radius);
%allnodes3D = allnodes3D(length(allnodes3D)+1-num_nodes:end,:)

nodes = zeros(length(allnodes),3);
p = levels;
k = 1;
for j = 1:levels
    for i = 1:num_nodes
        if allnodes(i,3) == p
            nodes(k,1:2) = allnodes(i,1:2);
            nodes(k,3) = j;
            k = k+1;
        end
    end
    p = p-1;
end

edges = zeros(num_edges,4);
p = 1;
for i = 2:levels-1
    for j = 1:num_nodes-1
            if nodes(j,3) == i
            edges(p,1)=sqrt((nodes(j+1,1)-nodes(j,1))^2+(nodes(j+1,2)-nodes(j,2))^2);       %Length of edge
            root_radius = root_radius*0.7;
            edges(p,2:4)=[  root_radius];
            p = p+1;
            end
    end
end
% Bruke levels på en eller annen måte?
edges
    
    
    