function [G, local2G, G2local] = makegraph(nodeconn, nsegment, nnode)
% Making a graph based on node connections. Returning the graph, and the
% transformations between a linear set of segments and the graph, e.g.
% additional values can be added to the graph by
% G.Edges.Feat = MyFeat(local2G);
%
if size(nodeconn,1) ~= nsegment
    warning('Multiple connections of at least one segment')
    nsegment = size(nodeconn,1);
end

v = cell(nnode,1);
for i = 1 : nnode
    v{i,1} = int2str(i);
end
NodeTable = table(v,'VariableNames',{'Name'});

index = (1:nsegment)';
edgetable = table(nodeconn, index, 'VariableNames',{'EndNodes', 'index'});
G = graph(edgetable, NodeTable);

% Transformations between graphs
local2G = G.Edges.index;    
[~,~, G2local] = intersect(index, G.Edges.index);