% % An example of rapidly-exploring random trees and path planning in 2-D
% % Ref: "Rapidly-Exploring Random Trees: A New Tool for Path Planning",
% % Steven M. LaValle, 1998
%~~~~
% Code can also be converted to function with input format
% [tree, path] = RRT(K, xMin, xMax, yMin, yMax, xInit, yInit, xGoal, yGoal, thresh)
% K is the number of iterations desired.
% xMin and xMax are the minimum and maximum values of x
% yMin and yMax are the minimum and maximum values of y
% xInit and yInit is the starting point of the algorithm
% xGoal and yGoal are the desired endpoints
% thresh is the allowed threshold distance between a random point the the
% goal point.
% Output is the tree structure containing X and Y vertices and the path
% found obtained from Init to Goal
%~~~~ 
% Written by: Omkar Halbe, Technical University of Munich, 31.10.2015
%~~~~ 
clear all; close all;
K=2000;
xMin=0; xMax=100;
yMin=0; yMax=100;
xInit=0; yInit=0; %initial point for planner
xGoal=100; yGoal=100; %goal for planner
thresh=5;           %acceptable threshold for solution
tree.vertex(1).x = xInit;
tree.vertex(1).y = yInit;
tree.vertex(1).xPrev = xInit;
tree.vertex(1).yPrev = yInit;
tree.vertex(1).dist=0;
tree.vertex(1).ind = 1; tree.vertex(1).indPrev = 0;
xArray=xInit; yArray = yInit;
figure(1); hold on; grid on;
plot(xInit, yInit, 'ko', 'MarkerSize',10, 'MarkerFaceColor','k');
plot(xGoal, yGoal, 'go', 'MarkerSize',10, 'MarkerFaceColor','g');
for iter = 2:K
    xRand = (xMax-xMin)*rand;
    yRand = (yMax-yMin)*rand;
    dist = Inf*ones(1,length(tree.vertex));
    for j = 1:length(tree.vertex)
        dist(j) = sqrt( (xRand-tree.vertex(j).x)^2 + (yRand-tree.vertex(j).y)^2 );
    end
    [val, ind] = min(dist);
       
    tree.vertex(iter).x = xRand; tree.vertex(iter).y = yRand;
    tree.vertex(iter).dist = val;
    tree.vertex(iter).xPrev = tree.vertex(ind).x;
    tree.vertex(iter).yPrev = tree.vertex(ind).y;
    tree.vertex(iter).ind = iter; tree.vertex(iter).indPrev = ind;
    
    if sqrt( (xRand-xGoal)^2 + (yRand-yGoal)^2 ) <= thresh
        plot([tree.vertex(iter).x; tree.vertex(ind).x],[tree.vertex(iter).y; tree.vertex(ind).y], 'r');
        break
    end
    
    plot([tree.vertex(iter).x; tree.vertex(ind).x],[tree.vertex(iter).y; tree.vertex(ind).y], 'r');
    pause(0);
end
if iter < K
    path.pos(1).x = xGoal; path.pos(1).y = yGoal;
    path.pos(2).x = tree.vertex(end).x; path.pos(2).y = tree.vertex(end).y;
    pathIndex = tree.vertex(end).indPrev;
    j=0;
    while 1
        path.pos(j+3).x = tree.vertex(pathIndex).x;
        path.pos(j+3).y = tree.vertex(pathIndex).y;
        pathIndex = tree.vertex(pathIndex).indPrev;
        if pathIndex == 1
            break
        end
        j=j+1;
    end
    path.pos(end+1).x = xInit; path.pos(end).y = yInit;
    for j = 2:length(path.pos)
        plot([path.pos(j).x; path.pos(j-1).x;], [path.pos(j).y; path.pos(j-1).y], 'b', 'Linewidth', 3);
    %     plot([tree.vertex(i).x; tree.vertex(ind).x],[tree.vertex(i).y; tree.vertex(ind).y], 'r');
    %     pause(0);
    end
else
    disp('No path found. Increase number of iterations and retry.');
end