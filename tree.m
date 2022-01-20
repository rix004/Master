function [nodes, edges]=tree(levels,root_radius,trunk_length,x1,y1,theta,r_rate,d_rate,delta_theta)
axis(gca,'equal')
num_nodes = 2^(levels-1);
num_edges = num_nodes-1;
x = [x1];
y = [y1];
edges_start = [0];
radius_start = [0];        % Hvor mye radiusen reduseres med for hvert nivå man går oppover

allnodes = newbranch(x,y,[levels],edges_start,radius_start,x1,y1,theta,trunk_length,levels,root_radius,r_rate,d_rate,delta_theta)
%allnodes3D = newbranch3D(x_val,y_val,z_val,x1,y1,z1,degrees,distance,levels,root_radius);

% Matrise med oversikt over alle nodene. Nodene nummereres fra roten og
% oppover, nivåvis, fra høyre til venstre.
% Kolonne 1 = nodenummer
% Kolonne 2 = x-koordinat
% Kolonne 3 = y-koordinat
% Kolonne 4 = nivå

nodes = zeros(length(allnodes),3);
p = levels;
k = 1;
for j = 1:levels
    for i = 1:num_nodes
        if allnodes(i,3) == p
            nodes(k,1) = k;
            nodes(k,2:3) = allnodes(i,1:2);
            nodes(k,4) = j;
            allnodes(i,6)=k;
            k = k+1;
        end
    end
    p = p-1;
end
% Matrise med oversikt over alle kantene. Kantene nummereres fra roten og
% oppver, nivåvis fra høyre til venstre.
% Kolonne 1 = lengde
% Kolonne 2 = node på nedside
% Kolonne 3 = node på overside
% Kolonne 4 = radius

edges = zeros(num_edges,4);
edges(1,:)=[sqrt((nodes(2,2)-nodes(1,2))^2+(nodes(2,3)-nodes(1,3))^2) 1 2 root_radius];
k=1;
p = levels;
for j = 1:levels
    for i = 2:num_nodes
        if allnodes(i,3) == p
            edges(k,1) = allnodes(i,4);
            edges(k,3) = k+1;
            edges(k,4) = allnodes(i,5);
            k = k+1;
        end
    end
    p = p-1;
end

g = 2;
f = 2;
for i = 1:num_edges/2
    edges(g,2) = f;
    edges(g+1,2) = f;
    g = g+2;
    f = f+1;
end
end


    
    
    