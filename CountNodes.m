function[NumNodes]=CountNodes(edges,MicroTermIndexes,MT)
Nodes2count = MicroTermIndexes;
for i = 1:length(MicroTermIndexes)
    Edge = find(MicroTermIndexes(i)==edges(:,3));
    NextNode = edges(Edge,2);
    if ismember(NextNode,Nodes2count)==0
         Nodes2count = [Nodes2count;NextNode];
    end
    while NextNode > MT
        Edge = find(NextNode==edges(:,3));
        NextNode = edges(Edge,2);
        if ismember(NextNode,Nodes2count)==0
            Nodes2count = [Nodes2count;NextNode];
        end
    end
    NumNodes = length(Nodes2count(1:end-1));
end