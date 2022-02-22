clc;clear;
% Number of nodes
Nn = 1000;

% Domain
xMin = -2; xMax = 2; yMin = -2; yMax = 2;

% Maximal length of edge in tree
lMax = sqrt((xMax-xMin)^2+(yMax-yMin)^2)/30;

% Start node
x_init = 0; y_init = 0;

%function[RRT_nodes,RRT_edges]=GenerateRRT(x_init,y_init,xMax,xMin,yMax,yMin,lMax, Nn)
RRT_nodes = zeros(Nn,2);
RRT_nodes(1,1:2) = [x_init y_init];

% Matrise med oversikt over alle kantene. Kantene nummereres fra roten og
% oppver, nivåvis fra høyre til venstre.
% Kolonne 1 = lengde
% Kolonne 2 = fra-node
% Kolonne 3 = til-node
% Kolonne 4 = radius
RRT_edges= zeros(Nn-1,4);
    for i = 2:Nn
        [x_r,y_r] = RandomState(xMax,xMin,yMax,yMin)
        [x_near,y_near,k] = NearestNeighbour(x_r,y_r,RRT_nodes);
        [x_new,y_new] = NewConf(x_near,y_near,x_r,y_r,lMax,RRT_edges,RRT_nodes);
        RRT_nodes(i,1:2) = [x_new,y_new];
        RRT_edges(i-1,2:3)=[k,i];
        line([x_near x_new],[y_near,y_new],'LineWidth',1,'Color',[0.8500, 0.3250, 0.0980]);
        hold on
        plot(x_new,y_new,'.b')
        hold on
        axis([xMin xMax yMin yMax])
        pause(1);
    end

function[x_rand,y_rand]=RandomState(xMax,xMin,yMax,yMin)
    x_rand = xMin + (xMax-xMin).*rand(1,1);
    y_rand = yMin + (yMax-yMin).*rand(1,1);
end

function[x_near,y_near,index]=NearestNeighbour(x,y,nodes)
dist = sqrt((x-nodes(1,1))^2+(y-nodes(1,2))^2);
x_near = nodes(1,1);
y_near = nodes(1,2);
index = 1;
    for i = 2:length(nodes)
        if sqrt((x-nodes(i,1))^2+(y-nodes(i,2))^2) < dist
            dist = sqrt((x-nodes(1,1))^2+(y-nodes(1,2))^2);
            x_near = nodes(i,1);
            y_near = nodes(i,2);
            index = i;
        end
    end
end

function[x_new,y_new]=NewConf(x_near,y_near,x_rand,y_rand,lMax,edges,nodes)
dist = sqrt((x_rand-x_near)^2+(y_rand-y_near)^2);
theta = acosd(abs(x_rand-x_near)/dist);
    if dist > lMax
        x_new = x_near + sign(x_rand-x_near)*lMax*cosd(theta);
        y_new = y_near + sign(y_rand-y_near)*lMax*sind(theta);
    elseif dist <= lMax
        x_new = x_rand;
        y_new = y_rand;
    end

% Check wether new edge intersect with any other edges

    for i = 1:length(edges)
        if edges(i,2) ~= 0
        x1 = nodes(edges(i,2),1);
        x2 = nodes(edges(i,3),1);
        y1 = nodes(edges(i,2),2);
        y2 = nodes(edges(i,3),2);
        newEdge = [x_new x_near;y_new y_near];
        oldEdge = [x1 x2;y1 y2];
        intrsct = LineIntersection(newEdge,oldEdge);
        p = 0;
            while intrsct == 1
                p = p+1
                x_new = x_new*0.9 %+ sign(x_rand-x_near)*lMax*0.5^p*cosd(theta)
                y_new = y_new*0.9; %+ sign(y_rand-y_near)*lMax*0.5^p*sind(theta)
                newEdge = [x_new x_near;y_new y_near]
                intrsct = LineIntersection(newEdge,oldEdge);
            end
        end
    end
end



