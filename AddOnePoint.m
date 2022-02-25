function[NewTree]=AddOnePoint(Tree,x,y)
% Find the node in Tree closest to (x,y). This point has node number k, and the distance to it from (x,y) is d1.
[x_near,y_near,k,d1] = NearestNode(x,y,Tree.nodes);

% Check if there are any edges. If it is not, a new edge is made from (x,y)
% to the closest (and only) point.
% If there are any edges, check if any edge is  closer to (x,y) than the 
% closest node. If so, the new edge is made as a normal on this edge. 
% if no edge is closer to (x,y) than the nearest node, 
% Kolonne 1: nodenummeret til den andre noden som danner kanten.
% Kolonne 2: Lengde p� kant
% Kolonne 3: Kantnummer
% Dersom kant-matrisen er tom, dvs. vi har kun en startnode, blir den nye
% kanten tegnet fra startnoden til punktet (x,y).
if size(Tree.edges,1) > 0
    [edge_nr,DistToEdge] = NearestEdge(x,y,Tree);
    if DistToEdge < d1
% Find point
            X1 = [Tree.nodes(Tree.edges(edge_nr,2),2) Tree.nodes(Tree.edges(edge_nr,2),3)];
            X2 = [Tree.nodes(Tree.edges(edge_nr,3),2) Tree.nodes(Tree.edges(edge_nr,3),3)];
            p = [x y];
            D1 = norm(p-X1);
            D2 = norm(X2-X1);
            D3 = norm(p-X2);
            theta = acosd((D1^2+D2^2-D3^2)/(2*D1*D2));
            ratio = D1*cosd(theta)/D2;
            x_new = X1(1) + ratio*(X2(1)-X1(1));
            y_new = X1(2) + ratio*(X2(2)-X1(2));
            l1 = sqrt((x_new-Tree.nodes(Tree.edges(edge_nr,2),2))^2+(y_new-Tree.nodes(Tree.edges(edge_nr,2),3))^2);
            l2 = Tree.edges(edge_nr,1)-l1;
            NewTree.edges = zeros(size(Tree.edges,1)+2,size(Tree.edges,2));
            NewTree.edges(1:edge_nr,:) = Tree.edges(1:edge_nr,:);
            NewTree.edges(edge_nr,:) = [l1 Tree.edges(edge_nr,2) Tree.nodes(end,1)+1 Tree.edges(edge_nr,4)];
            NewTree.edges(edge_nr+1,:) = [l2 Tree.nodes(end,1)+1 Tree.edges(edge_nr,3) Tree.edges(edge_nr,4)];
            NewTree.edges(edge_nr+2:end-1,:) = Tree.edges(edge_nr+1:end,:);
            NewTree.edges(end,:) = [DistToEdge Tree.nodes(end,1)+1 Tree.nodes(end,1)+2 Tree.edges(edge_nr,4)*Tree.radius_ratio];
            NewTree.nodes = [Tree.nodes;Tree.nodes(end,1)+1 x_new y_new Tree.nodes(Tree.edges(edge_nr,2),4); Tree.nodes(end,1)+2 x y Tree.nodes(Tree.edges(edge_nr,2),4)+1];
    elseif DistToEdge >= d1
        check_edges = [];
        for i = 1:size(Tree.edges,1)
            if Tree.edges(i,2) == k
                n2 = Tree.edges(i,3);
                d2 = Tree.edges(i,1);
                d3 = sqrt((x-Tree.nodes(n2,2))^2+(y-Tree.nodes(n2,3))^2);
                theta = acosd((d1^2+d2^2-d3^2)/(2*d1*d2));
                check_edges = [check_edges;n2 d2 i theta];
            elseif Tree.edges(i,3) == k
                n2 = Tree.edges(i,2);
                d2 = Tree.edges(i,1);
                d3 = sqrt((x-Tree.nodes(n2,2))^2+(y-Tree.nodes(n2,3))^2);
                theta = acosd((d1^2+d2^2-d3^2)/(2*d1*d2));
                check_edges = [check_edges;n2 d2 i theta];
            end
        end

        % Sjekker hvilken kant som er n�rmest ved � sjekke hvilken kant som danner
        % minste vinkel med vektoren fra (x,y) til (x_near,y_near).
        angles = check_edges(:,4);
        theta = min(angles);
        index = find(angles==theta);
        edge_nr = check_edges(index,3);
        d2 = check_edges(index,2);
        n2_koord = check_edges(index,1);
        ratio = d1*cosd(theta)/d2;

        % Dersom den minste vinkelen er mindre enn 90 grader blir den nye
        % kanten en normal p� n�rmeste kant. Dersom vinkelen er st�rre enn 90
        % grader betyr det at n�rmeste node er en terminalnode, og den nye
        % kanten blir en ny "terminalkant". 
        if theta < 90
            x_new = x_near + ratio*(Tree.nodes(n2_koord,2)-x_near);
            y_new = y_near + ratio*(Tree.nodes(n2_koord,3)-y_near);
            l1 = sqrt((x_new-Tree.nodes(Tree.edges(edge_nr,2),2))^2+(y_new-Tree.nodes(Tree.edges(edge_nr,2),3))^2);
            l2 = Tree.edges(edge_nr,1)-l1;
            NewTree.edges = zeros(size(Tree.edges,1)+2,size(Tree.edges,2));
            NewTree.edges(1:edge_nr,:) = Tree.edges(1:edge_nr,:);
            NewTree.edges(edge_nr,:) = [l1 Tree.edges(edge_nr,2) Tree.nodes(end,1)+1 Tree.edges(edge_nr,4)];
            NewTree.edges(edge_nr+1,:) = [l2 Tree.nodes(end,1)+1 Tree.edges(edge_nr,3) Tree.edges(edge_nr,4)];
            NewTree.edges(edge_nr+2:end-1,:) = Tree.edges(edge_nr+1:end,:);
            NewTree.edges(end,:) = [d1*sind(theta) Tree.nodes(end,1)+1 Tree.nodes(end,1)+2 Tree.edges(edge_nr,4)*Tree.radius_ratio];
            NewTree.nodes = [Tree.nodes;Tree.nodes(end,1)+1 x_new y_new Tree.nodes(Tree.edges(edge_nr,2),4); Tree.nodes(end,1)+2 x y Tree.nodes(Tree.edges(edge_nr,2),4)+1];
        else
            NewTree.nodes = [Tree.nodes;Tree.nodes(end,1)+1 x y Tree.nodes(Tree.edges(edge_nr,3),4)+1];
            NewTree.edges = [Tree.edges;d1 k Tree.nodes(end,1)+1 Tree.edges(edge_nr,4)];
        end
    else
        disp('shady ting skjer 1')
    end
elseif size(Tree.edges,1) == 0
        NewTree.nodes = [Tree.nodes;Tree.nodes(end,1)+1 x y Tree.nodes(end,4)+1];
        NewTree.edges = [Tree.edges;d1 k k+1 Tree.root_radius];
else
    disp('shady ting skjer 2')
end

end