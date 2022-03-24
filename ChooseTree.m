function[Tree]=ChooseTree(t,RandomTree,DT,Domain,Np)
if strcmp(t, 'Deterministic')
    Tree= GetTree(DT);
elseif strcmp(t, 'Random')
    RandomTree.nodes = [RootNode 1];
    RandomTree.edges = [];
    RandomTree = RRT_Tree(RandomTree,Domain,Np,0);
    Tree = RandomTree;
elseif strcmp(t, 'Combinated')
    DetTree= GetTree(DT);
    RandomTree.nodes = DetTree.nodes;
    RandomTree.edges = DetTree.edges;
    Tree = RRT_Tree(RandomTree,Domain,Np,0);
elseif strcmp(t, 'Half Deterministic')
    DetTree= GetTree(DT);
    pos_nodes = find(DetTree.nodes(:,1)>=0);
    neg_nodes = find(DetTree.nodes(:,1)<0);
    half_nodes = DetTree.nodes(pos_nodes,:);
    pos_edges = ismember(DetTree.edges(:,3),pos_nodes);
    half_edges = DetTree.edges(pos_edges==1,:);
    for i = 1:length(half_edges)
        ii = find(pos_nodes==half_edges(i,3));
        jj = find(pos_nodes==half_edges(i,2));
        half_edges(i,3)=ii;
        half_edges(i,2)=jj;
    end
    Tree.nodes = half_nodes;
    Tree.edges = half_edges;
end
end

