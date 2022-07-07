function[tree]=RRT_Tree(tree,D,Np)
% nodes: Matrix with all the nodes in the tree. 
% Column 1: x-coordinate
% Column 2: y-coordinate
% Column 3: level

% edges: Matrix with all the edges in the tree. 
% Column 1 = length
% Column 2 = from-node (parent)
% Column 3 = to-node   (child)
% Column 4 = radius

num_cells = 0;
while num_cells < Np
    [x_r,y_r] = RandomState(D);
    NewTree = AddOnePoint(tree,x_r,y_r);
    tree.nodes = NewTree.nodes;
    tree.edges = NewTree.edges;
    [Tn,~] = FindTerminals(tree.nodes,tree.edges);
    num_cells = size(Tn,1);
end
end
