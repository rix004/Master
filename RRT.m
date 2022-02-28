function[tree]=RRT_Tree(tree,xMin,xMax,yMin,yMax,Ni)
clc;clear;
% Number of iterations
%Nn = 500;

% Domain
%xMin = -2; xMax = 2; yMin = -2; yMax = 2;

% Start node
%x_init = 0; y_init = 0;

% nodes: Matrix with all the nodes in the tree. 
% Column 1: x-coordinate
% Column 2: y-coordinate
% Column 3: level

% edges: Matrix with all the edges in the tree. 
% Column 1 = length
% Column 2 = from-node
% Column 3 = to-node
% Column 4 = radius

RRT_nodes = [];
RRT_nodes(1,:) = [x_init y_init 1];
tree.nodes = RRT_nodes;
tree.edges = [];
tree.root_radius = 2;
tree.radius_ratio = 0.5;

    for i = 2:Nn
        [x_r,y_r] = RandomState(xMax,xMin,yMax,yMin);
        NewTree = AddOnePoint(tree,x_r,y_r);
        tree.nodes = NewTree.nodes;
        tree.edges = NewTree.edges;
    end
    
    DrawTree(tree)


function[x_rand,y_rand]=RandomState(xMax,xMin,yMax,yMin)
    x_rand = xMin + (xMax-xMin).*rand(1,1);
    y_rand = yMin + (yMax-yMin).*rand(1,1);
end





end
