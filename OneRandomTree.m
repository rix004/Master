function[tree]=OneRandomTree(tree,D,Np)

DomainSize = D(2)-D(1);
scale = DomainSize-D(2);
rng(131295);
coords = rand(Np*3,2)*DomainSize - scale;
%coords(end,:)=[mean([D(1) D(2)]) mean([D(3) D(4)])];
num_cells = 0;
i = 1;
while num_cells<Np
    NewTree = AddOnePoint(tree,coords(i,1),coords(i,2));
    tree.nodes = NewTree.nodes;
    tree.edges = NewTree.edges;
    [Tn,~] = FindTerminals(tree.nodes,tree.edges);
    num_cells = size(Tn,1);
    i = i+1;
end

end
