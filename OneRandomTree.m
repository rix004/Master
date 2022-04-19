function[tree]=OneRandomTree(tree,D,Np)

DomainSize = D(2)-D(1);
scale = DomainSize-D(2);
rng(141295);
coords = rand(Np,2)*DomainSize - scale;
coords(end,:)=[mean([D(1) D(2)]) mean([D(3) D(4)])];
for i = 1:Np
    NewTree = AddOnePoint(tree,coords(i,1),coords(i,2));
    tree.nodes = NewTree.nodes;
    tree.edges = NewTree.edges;
end

end
