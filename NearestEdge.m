function[edge_nr,d] = NearestEdge(x,y,nodes,edges)
p = [x y];

% Check first edge
x1 = [nodes(edges(1,2),1) nodes(edges(1,2),2)];
x2 = [nodes(edges(1,3),1) nodes(edges(1,3),2)];
d1 = norm(p-x1);
d2 = norm(x2-x1);
d3 = norm(p-x2);
if d1^2 < d3^2+d2^2 && d3^2 < d1^2+d2^2
    theta = acosd((d1^2+d2^2-d3^2)/(2*d1*d2));
    d = abs(sind(theta)*d1);
    edge_nr = 1;
else
    d = 10E10;
    edge_nr = 0;
end


% Check if any other edge is closer
for i = 2:size(edges,1)
    X1 = [nodes(edges(i,2),1) nodes(edges(i,2),2)];
    X2 = [nodes(edges(i,3),1) nodes(edges(i,3),2)];
    D1 = norm(p-X1);
    D2 = norm(X2-X1);
    D3 = norm(p-X2);
    if abs(sind(acosd((D1^2+D2^2-D3^2)/(2*D1*D2)))*D1) < d && D1^2 < D3^2+D2^2 && D3^2 < D1^2+D2^2
        theta = acosd((D1^2+D2^2-D3^2)/(2*D1*D2));
        d = abs(sind(theta)*D1);
        edge_nr = i;
        x1 = X1;
        x2 = X2;
        d1 = D1;
        d2 = D2;
        d3 = D3;
    end
end
end