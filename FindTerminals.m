function[Tn,Ext_nodes] = FindTerminals(Tree)
nodes = Tree.nodes;
edges = Tree.edges;
Ne = length(edges);
Nn = length(nodes);
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
% if Tree.RootNodeIdx ~= 1
%     find_tn1 = ismember(StartNodes,EndNodes);
%     find_tn1(Tree.RootNodeIdx)=1;
%     for i = 1:length(StartNodes)
%         if find_tn1(i) == 0
%             if sum(ismember(StartNodes,StartNodes(i)))<2
%             Tn = [Tn;nodes(StartNodes(i),1:2),i];
%             Ext_nodes(StartNodes(i))=1;
%             end
%         end
%     end
% end
end