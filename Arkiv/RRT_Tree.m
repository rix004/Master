function[tree]=RRT_Tree(tree,D,Np,ForbTree)
% nodes: Matrix with all the nodes in the tree. 
% Column 1: x-coordinate
% Column 2: y-coordinate
% Column 3: level

% edges: Matrix with all the edges in the tree. 
% Column 1 = length
% Column 2 = from-node
% Column 3 = to-node
% Column 4 = radius

if ForbTree == 0
    for i = 1:Np
        [x_r,y_r] = RandomState(D);
        NewTree = AddOnePoint(tree,x_r,y_r);
        tree.nodes = NewTree.nodes;
        tree.edges = NewTree.edges;
    end
elseif ForbTree == 1
    StartNodes = tree.edges(:,2);
    EndNodes = tree.edges(:,3);
    find_tn = ismember(EndNodes,StartNodes);
    Tn = [];
    for i = 1:length(EndNodes)
        if find_tn(i) == 0
            Tn = [Tn;tree.nodes(EndNodes(i),:)];
        end
    end
    ForbNodesIndex = ismember(tree.nodes,Tn);
    ForbNodes = tree.nodes(ForbNodesIndex(:,1)==0,:);
    ForbiddenTree.nodes = ForbNodes;
    ForbiddenTree.edges = tree.edges;
    ForbNodes
    tree.edges
    for i = 1:Np
        [x_r,y_r] = RandomState(D);
        NewTree = AddOnePoint1(tree,ForbiddenTree,x_r,y_r);
        tree.nodes = NewTree.nodes;
        tree.edges = NewTree.edges;
    end
end

function[x_rand,y_rand]=RandomState(D)
    x_rand = D(1) + (D(2)-D(1)).*rand(1,1);
    y_rand = D(3) + (D(4)-D(3)).*rand(1,1);
end





end
