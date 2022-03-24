function [tree]=GetTree(v)
x = [v.StartPos(1)];
y = [v.StartPos(2)];
edges_start = [0];
radius_start = [0];

allnodes = MakeTree(x,y,[v.Levels],edges_start,radius_start,v.StartPos(1),v.StartPos(2),v.StartAngle,...
                    v.TrunkLength,v.Levels,v.TrunkRadius,v.RadiusRate,v.LengthRate,v.RotationAngle);


num_nodes = 2^(v.Levels-1);
num_edges = num_nodes-1;

% Matrise med oversikt over alle nodene. Nodene nummereres fra roten og
% oppover, nivåvis, fra høyre til venstre.
% Kolonne 1 = x-koordinat
% Kolonne 2 = y-koordinat
% Kolonne 3 = nivå
nodes = zeros(size(allnodes,1),3);
p = v.Levels;
k = 1;
for j = 1:v.Levels
    for i = 1:num_nodes
        if allnodes(i,3) == p
            nodes(k,1:2) = allnodes(i,1:2);
            nodes(k,3) = j;
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
edges(1,:)=[sqrt((nodes(2,2)-nodes(1,2))^2+(nodes(2,3)-nodes(1,3))^2) 1 2 v.TrunkRadius];
k=1;
p = v.Levels;
for j = 1:v.Levels
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
tree.nodes = nodes;
tree.edges = edges;
end


    
    
    