function [nodes,edges]= GetMatrices3D(levels,R_trunk,L_trunk,x0,y0,z0,theta,phi,R_rate,L_rate,delta_theta,delta_phi)
axis(gca,'equal');

num_nodes = 2^(levels-1);
num_edges = num_nodes-1;
x = [x0];
y = [y0];
z = [z0];
edges_start = [0];
radius_start = [0];        % Hvor mye radiusen reduseres med for hvert nivå man går oppover

allnodes = tree3D(levels,x,y,z,[levels],edges_start,radius_start,x0,y0,z0,theta,phi,L_trunk,R_trunk,R_rate,L_rate,delta_theta,delta_phi);


% Matrise med oversikt over alle nodene. Nodene nummereres fra roten og
% oppover, nivåvis, fra høyre til venstre.
% Kolonne 1 = nodenummer
% Kolonne 2 = x-koordinat
% Kolonne 3 = y-koordinat
% Kolonne 4 = z-koordinat
% Kolonne 5 = nivå

nodes = zeros(length(allnodes),4);
p = levels;
k = 1;
for j = 1:levels
    for i = 1:num_nodes
        if allnodes(i,4) == p
            nodes(k,1) = k;
            nodes(k,2:4) = allnodes(i,1:3);
            nodes(k,5) = j;
            allnodes(i,7)=k;
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
edges(1,:)=[sqrt((nodes(2,2)-nodes(1,2))^2+(nodes(2,3)-nodes(1,3))^2+(nodes(2,4)-nodes(1,4))^2) 1 2 R_trunk];
k=1;
p = levels;
for j = 1:levels
    for i = 2:num_nodes
        if allnodes(i,4) == p
            edges(k,1) = allnodes(i,5);
            edges(k,3) = k+1;
            edges(k,4) = allnodes(i,6);
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
    
    
    