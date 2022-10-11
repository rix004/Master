function[Tree]=ChooseTree(t,RandomTree,DT,DLA,Domain)
if strcmp(t, 'Deterministic')
    Tree= GetTree(DT);
    Tree.nodes = FindLevels(Tree.nodes,Tree.edges);
    Tree.RootNodeIdx = 1;
    Tree.RadiusRate = DT.RadiusRate;
elseif strcmp(t, 'Random')
    RandomTree.nodes = [0.5 0.5];
    RandomTree.edges = [];
    RandomTree = RRT_Tree(RandomTree,Domain,Np);
    Tree = RandomTree;
elseif strcmp(t, 'Combinated')
    DetTree= GetTree(DT);

    % Find terminal nodes
    [Tn,~]=FindTerminals(DetTree);
    
    % Initiate random tree
    RandomTree.nodes = Tn(:,1:2);
    RandomTree.edges = [];
    Np = RandomTree.Ncells;
    RandomTree.TrunkRadius=DetTree.edges(end,4)*RandomTree.RadiusRate;
    RandomTree = RRT_Tree(RandomTree,Domain,Np);

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
    Tree.nodes = FindLevels(Tree.nodes,Tree.edges);
    Tree.RadiusRate = RandomTree.RadiusRate;
    [TNinfo,TNlogic]=FindTerminals(Tree);
    while size(TNinfo,1)>Np
        [Tree.nodes,Tree.edges]=CutTree(Tree.nodes,Tree.edges,max(Tree.nodes(:,3)));
        Tree.nodes = FindLevels(Tree.nodes,Tree.edges);
        [TNinfo,TNlogic]=FindTerminals(Tree);
        size(TNinfo,1);
    end
    while size(TNinfo,1)<Np
        [x_r,y_r] = RandomState(Domain);
        Tree.nodes = Tree.nodes(:,1:2);
        Tree = AddOnePoint(Tree,x_r,y_r);
        [TNinfo,TNlogic]=FindTerminals(Tree);
        Tree.nodes = FindLevels(Tree.nodes,Tree.edges);
        size(TNinfo,1);
    end
    Tree.RootNodeIdx = 1;
    
    % Fix edge radii

    % Define radii
    RootRadius = DT.TrunkRadius;
    TerminalRadius = RandomTree.TerminalRadius;
    
    % Make temporary tree
    TempTree.edges = Tree.edges(2:end,:);
    TempTree.nodes = Tree.nodes(2:end,:);
    TempTree.edges(:,2:3)=TempTree.edges(:,2:3)-1;
    TempTree.nodes(:,3)=TempTree.nodes(:,3)-1;
    %DrawTree(TempTree,100,'r',Domain);
    %text(TempTree.nodes(:,1),TempTree.nodes(:,2),num2str(TempTree.nodes(:,3)),'FontSize',15);

    [~,TNtemp]=FindTerminals(TempTree);
    MicroTermIndexes = find(TNtemp==1);
    [Graph, ~, ~]=makegraph(TempTree.edges(:,2:3),size(TempTree.edges,1),size(TempTree.nodes,1));
    NodeLevels = zeros(size(TempTree.nodes,1),1);
    for i = 1:length(MicroTermIndexes)
        SP = shortestpath(Graph,MicroTermIndexes(i),1);
        for j = 1:size(TempTree.nodes,1)
            if any(ismember(SP,j)) && NodeLevels(j) < find(ismember(SP,j))-1
                NodeLevels(j) = find(ismember(SP,j))-1;
            end
        end
    end

    % Set radii
    power = log(RootRadius/TerminalRadius)/(max(TempTree.nodes(:,3))-1);
    RadiusRate = exp(power);
    for i = 2:length(NodeLevels)
        edge = find(TempTree.edges(:,3)==i);
        edgerad = TerminalRadius*RadiusRate^(NodeLevels(i));
        TempTree.edges(edge,4)=edgerad;
    end

    TempTree.edges(:,2:3)=TempTree.edges(:,2:3)+1;
    Tree.edges = [Tree.edges(1,:);TempTree.edges];
    
elseif strcmp(t, 'DLA')
    setpath;
    DLAtree = npy2tree(DLA.Nparticles,DLA.version,DLA.RootRadius,DLA.TerminalRadius);
    Tree.edges = DLAtree.edges;
    Tree.nodes = DLAtree.nodes;
    Tree.RootNodeIdx = 1;

    % Increase level number before the new root node enters
    Tree.nodes(:,3)=Tree.nodes(:,3)+1;

    % Create new root node and add to nodes at the beginning
    NewRootNode = [Tree.nodes(1,1) Domain(3) 0];
    NewNodes = zeros(size(Tree.nodes,1)+1,3);
    NewNodes(2:end,:)=Tree.nodes;
    NewNodes(1,:)= NewRootNode;
    Tree.nodes = NewNodes;
    
    % Create new root edge and add to edges.
    Tree.edges(:,2)=Tree.edges(:,2)+1;
    Tree.edges(:,3)=Tree.edges(:,3)+1;

    NewEdges = zeros(size(Tree.edges,1)+1,4);
    NewEdgeLength = sqrt((Tree.nodes(1,1)-Tree.nodes(2,1))^2+(Tree.nodes(1,2)-Tree.nodes(2,2))^2);
    NewEdges(2:end,:)=Tree.edges;
    NewEdges(1,:)=[NewEdgeLength 1 2 max(Tree.edges(:,4))];
    Tree.edges = NewEdges;
end
end

