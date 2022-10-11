function[NewNodes]=FindLevels(nodes,edges)
NewNodes = zeros(size(nodes,1),3);
NewNodes(:,1:2)=nodes(:,1:2);
G = graph(edges(:,2),edges(:,3));
D = distances(G);
NewNodes(:,3)=D(:,1);
end