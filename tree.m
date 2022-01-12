clc;clear;
format short
axis(gca,'equal')
distance = 0.7;
levels = 4;
num_nodes = 2^(levels-1);
num_edges = num_nodes-1;
x1 = 0;
y1 = 0;
z1 = 0;
degrees = 90;
root_radius = 5;
x = [x1];
y = [y1];
z = [z1];
edges_start = [0];
radius_start = [0];
min_rate = 0.7;         % Hvor mye radiusen reduseres med for hvert nivå man går oppover

allnodes = newbranch(x,y,[levels],edges_start,radius_start,x1,y1,degrees,distance,levels,root_radius,min_rate)
%allnodes3D = newbranch3D(x_val,y_val,z_val,x1,y1,z1,degrees,distance,levels,root_radius);
%allnodes3D = allnodes3D(length(allnodes3D)+1-num_nodes:end,:)

% Matrise med oversikt over alle nodene. Nodene nummereres fra roten og
% oppover, levelvis, fra høyre til venstre. 
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
edges(1,:)=[sqrt((nodes(2,1)-nodes(1,1))^2+(nodes(2,2)-nodes(1,2))^2) 1 2 root_radius];
k=1;
p = levels;
for j = 1:levels
    for i = 2:num_nodes
        if allnodes(i,3) == p
            edges(k,1) = allnodes(i,4);
            edges(k,4) = allnodes(i,5);
            k = k+1;
        end
    end
    p = p-1;
end


    
    
    