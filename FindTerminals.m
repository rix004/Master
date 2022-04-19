function[Tn,Te,Ext_nodes] = FindTerminals(nodes,edges)
Ne = length(edges);
Nn = length(nodes);
StartNodes = edges(:,2);
EndNodes = edges(:,3);
find_tn = ismember(EndNodes,StartNodes);
Tn = [];
Te = [];
Ext_nodes = [zeros(Nn,1)];
for i = 1:length(EndNodes)
    if find_tn(i) == 0
        Tn = [Tn;nodes(EndNodes(i),1:2),i];
        Te = [Te;i];
        Ext_nodes(EndNodes(i))=1; 
    end
end
end