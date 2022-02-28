function[tree]=RRT_Tree(tree,xMin,xMax,yMin,yMax,Np)
% nodes: Matrix with all the nodes in the tree. 
% Column 1: x-coordinate
% Column 2: y-coordinate
% Column 3: level

% edges: Matrix with all the edges in the tree. 
% Column 1 = length
% Column 2 = from-node
% Column 3 = to-node
% Column 4 = radius

    for i = 1:Np
        [x_r,y_r] = RandomState(xMax,xMin,yMax,yMin);
        NewTree = AddOnePoint(tree,x_r,y_r);
        tree.nodes = NewTree.nodes;
        tree.edges = NewTree.edges;
    end

function[x_rand,y_rand]=RandomState(xMax,xMin,yMax,yMin)
    x_rand = xMin + (xMax-xMin).*rand(1,1);
    y_rand = yMin + (yMax-yMin).*rand(1,1);
end





end
