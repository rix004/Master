function[nodes]=FindLevels(nodes,edges)
nodes(1,3)=1;
for i = 2:size(nodes,1)
    levels = 1;
    node = i;
    while node > 1
    e = find(node == edges(:,3));
    node = edges(e,2);
    levels = levels + 1;
    end
    nodes(i,3)=levels;
end
end