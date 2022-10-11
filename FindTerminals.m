function[TNinfo1,Ext_nodes] = FindTerminals(Tree)
nodes = Tree.nodes;
edges = Tree.edges;
Ne = size(edges,1);
Nn = size(nodes,1);
if Ne == 0
    Tn = Tree.nodes;
    Ext_nodes = 1;
else
    StartNodes = edges(:,2);
    EndNodes = edges(:,3);
    find_tn = ismember(EndNodes,StartNodes);
    Tn = [];
    Ext_nodes = [zeros(Nn,1)];
    for i = 1:length(EndNodes)
        if find_tn(i) == 0
            if sum(ismember(EndNodes,EndNodes(i)))<2
                Tn = [Tn;nodes(EndNodes(i),1:2),i];
                Ext_nodes(EndNodes(i))=1; 
            end
        end
    end
    TNinfo1 = Tn;
    count = 0;
    for i = 1:length(Ext_nodes)
        if Ext_nodes(i) == 1
            count = count + 1;
            TNinfo1(count,1:2)=nodes(i,1:2);
            edgenr = Tn(find(ismember(Tn(:,1:2),nodes(i,1:2),'rows')),3);
            TNinfo1(count,3)=edgenr;
        end
end

end