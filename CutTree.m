function[newnodes,newedges]=CutTree(nodes,edges,pos)
CutNodes = find(nodes(:,3)>=pos);
newnodes = nodes;
newnodes(CutNodes,:)=[];
CutEdges = [];
for i = 1:length(CutNodes)
    CutEdges = [CutEdges;find(edges(:,3)==CutNodes(i))];
end
newedges = edges;
newedges(CutEdges,:)=[];
% Rearrange node numbers in edges
for i = 1:size(newedges,1)
        NewEdgeNr1=find(ismember(newnodes(:,1:2),nodes(newedges(i,2),1:2),"rows"));
        NewEdgeNr2=find(ismember(newnodes(:,1:2),nodes(newedges(i,3),1:2),"rows"));
        newedges(i,2) = NewEdgeNr1;
        newedges(i,3) = NewEdgeNr2;
end
end