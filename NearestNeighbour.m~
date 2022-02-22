% NearestNeighbour(x,y,nodes) finner den noden fra matrisen nodes som ligger nærmest punktet (x,y).
% x_near,y_near er koordinatene til den næremste noden fra nodes.
% index angir hvilken plassering den nærmeste noden har i nodes.
% distance angir avstand fra punktet (x,y) til punktet (x_near, y_near).


function[x_near,y_near,index,distance]=NearestNeighbour(x,y,nodes)
distance = sqrt((x-nodes(1,2))^2+(y-nodes(1,3))^2);
x_near = nodes(1,2);
y_near = nodes(1,3);
index = 1;
if length(nodes(:,1)) == 1
    distance = sqrt((x-nodes(1,2))^2+(y-nodes(1,3))^2);
    x_near = nodes(1,2);
    y_near = nodes(1,3);
    index = 1;
else
    for i = 2:length(nodes(:,1))
        if sqrt((x-nodes(i,2))^2+(y-nodes(i,3))^2) < distance
            distance = sqrt((x-nodes(i,2))^2+(y-nodes(i,3))^2);
            x_near = nodes(i,2);
            y_near = nodes(i,3);
            index = nodes(i,1);
        end
    end
end
end