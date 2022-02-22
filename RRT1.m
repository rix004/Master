clc;clear;
% Number of nodes
Nn = 1000;

% Domain
xMin = -2; xMax = 2; yMin = -2; yMax = 2;

% Start node
x_init = 0; y_init = 0;

%function[RRT_nodes,RRT_edges]=GenerateRRT(x_init,y_init,xMax,xMin,yMax,yMin,lMax, Nn)
RRT_nodes = [];
RRT_nodes(1,:) = [1 x_init y_init 1];
tree.nodes = RRT_nodes;
tree.edges = [];
tree.root_radius = 2;
tree.radius_ratio = 0.5;

% Matrise med oversikt over alle kantene. Kantene nummereres fra roten og
% oppover, nivåvis fra høyre til venstre.
% Kolonne 1 = lengde
% Kolonne 2 = fra-node
% Kolonne 3 = til-node
% Kolonne 4 = radius
    for i = 2:Nn
        [x_r,y_r] = RandomState(xMax,xMin,yMax,yMin);
        NewTree = AddOnePoint(tree,x_r,y_r);
        tree.nodes = NewTree.nodes;
        tree.edges = NewTree.edges;
        tree.nodes;
        tree.edges;
    end
for i = 1:size(tree.edges,1)
    x1 = tree.nodes(tree.edges(i,2),2);
    x2 = tree.nodes(tree.edges(i,3),2);
    y1 = tree.nodes(tree.edges(i,2),3);
    y2 = tree.nodes(tree.edges(i,3),3);
    line([x1 x2],[y1,y2],'LineWidth',tree.edges(i,4),'Color',[0.8500, 0.3250, 0.0980])
    hold on
    axis([xMin xMax yMin yMax])
end

function[x_rand,y_rand]=RandomState(xMax,xMin,yMax,yMin)
    x_rand = xMin + (xMax-xMin).*rand(1,1);
    y_rand = yMin + (yMax-yMin).*rand(1,1);
end






