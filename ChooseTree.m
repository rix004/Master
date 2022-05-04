function[Tree]=ChooseTree(t,RandomTree,DT,Domain,Np)
if strcmp(t, 'Deterministic')
    Tree= GetTree(DT);
    DrawTree(Tree,10,[0.8500, 0.3250, 0.0980]);
elseif strcmp(t, 'Random')
    RandomTree.nodes = [RootNode 1];
    RandomTree.edges = [];
    RandomTree = RRT_Tree(RandomTree,Domain,Np,0);
    Tree = RandomTree;
    DrawTree(Tree,5,[0 0.4470 0.7410]);
elseif strcmp(t, 'Combinated')
    DetTree= GetTree(DT);

    % Find terminal nodes
    [Tn,~,~]=FindTerminals(DetTree.nodes,DetTree.edges);

    % Initiate random tree
    RandomTree.nodes = Tn;
    RandomTree.edges = [];
    RandomTree.TrunkRadius=DetTree.edges(end,4)*RandomTree.RadiusRate;
    RandomTree = RRT_Tree(RandomTree,Domain,Np);

    % Draw trees
    figure(1)
    DrawTree(DetTree,50,[0.8500, 0.3250, 0.0980]);
    DrawTree(RandomTree,1,[0 0.4470 0.7410]);

    % Define the combinated tree
    Tree.nodes=[DetTree.nodes;RandomTree.nodes];
    Tree.nodes = unique(Tree.nodes(:,1:2),"rows",'stable');
    Tree.edges(:,1)=[DetTree.edges(:,1);RandomTree.edges(:,1)];
    for i = 1:size(RandomTree.edges,1)
        NewEdgeNr1=find(ismember(Tree.nodes(:,1:2),RandomTree.nodes(RandomTree.edges(i,2),1:2),"rows"));
        NewEdgeNr2=find(ismember(Tree.nodes(:,1:2),RandomTree.nodes(RandomTree.edges(i,3),1:2),"rows"));
        RandomTree.edges(i,2) = NewEdgeNr1;
        RandomTree.edges(i,3) = NewEdgeNr2;
    end
    Tree.edges(:,2:4) = [DetTree.edges(:,2:4);RandomTree.edges(:,2:4)];
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
    DrawTree(Tree,10,[0.8500, 0.3250, 0.0980]);
end
end

