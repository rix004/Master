function[NewTree]=AddOnePoint(Tree,x,y)
% Find the node in Tree closest to (x,y). This point has node number k, and the distance to it from (x,y) is d1.
[x_near,y_near,k,d1] = NearestNode(x,y,Tree.nodes);

% Check if there are any edges. If not, a new edge is made from (x,y)
% to the closest node.

% If there are any edges, check if any edge is  closer to (x,y) than the 
% closest node. If so, the new edge is made as a normal to this edge. 
if size(Tree.edges,1) > 0
    [edge_nr,DistToEdge] = NearestEdge(x,y,Tree.nodes,Tree.edges);
    if DistToEdge < d1
            % Find point on edge
            X1 = [Tree.nodes(Tree.edges(edge_nr,2),1) Tree.nodes(Tree.edges(edge_nr,2),2)];
            X2 = [Tree.nodes(Tree.edges(edge_nr,3),1) Tree.nodes(Tree.edges(edge_nr,3),2)];
            p = [x y];
            D1 = norm(p-X1);
            D2 = norm(X2-X1);
            D3 = norm(p-X2);
            theta = acosd((D1^2+D2^2-D3^2)/(2*D1*D2));
            ratio = D1*cosd(theta)/D2;
            x_new = X1(1) + ratio*(X2(1)-X1(1));
            y_new = X1(2) + ratio*(X2(2)-X1(2));
            
            % Find lenght of the new edges (a new node on an already existing edge will split up this edge in two)
            l1 = sqrt((x_new-Tree.nodes(Tree.edges(edge_nr,2),1))^2+(y_new-Tree.nodes(Tree.edges(edge_nr,2),2))^2);
            l2 = Tree.edges(edge_nr,1)-l1;
            
            % Update tree
            NewTree.edges = zeros(size(Tree.edges,1)+2,size(Tree.edges,2));
            NewTree.edges(1:edge_nr,:) = Tree.edges(1:edge_nr,:);
            NewTree.edges(edge_nr,:) = [l1 Tree.edges(edge_nr,2) size(Tree.nodes,1)+1 Tree.edges(edge_nr,4)];
            NewTree.edges(edge_nr+1,:) = [l2 size(Tree.nodes,1)+1 Tree.edges(edge_nr,3) Tree.edges(edge_nr,4)];
            NewTree.edges(edge_nr+2:end-1,:) = Tree.edges(edge_nr+1:end,:);
            NewTree.edges(end,:) = [DistToEdge size(Tree.nodes,1)+1 size(Tree.nodes,1)+2 Tree.edges(edge_nr,4)*Tree.RadiusRate];
            NewTree.nodes = [Tree.nodes;x_new y_new;x y];
    elseif DistToEdge >= d1
    % If no edge is closer to (x,y) than the nearest node, we check all the
    % edges belonging to the nearest node. The data of these edges are stored
    % in the matrix check_edges, structured in the following way: 
    % Column 1: The node number of the other node making up the edge.
    % Column 2: Length of edge
    % Column 3: Edge number
        check_edges = [];
        for i = 1:size(Tree.edges,1)
            if Tree.edges(i,2) == k
                n2 = Tree.edges(i,3);
                d2 = Tree.edges(i,1);
                d3 = sqrt((x-Tree.nodes(n2,1))^2+(y-Tree.nodes(n2,2))^2);
                theta = acosd((d1^2+d2^2-d3^2)/(2*d1*d2));
                check_edges = [check_edges;n2 d2 i theta];
            elseif Tree.edges(i,3) == k
                n2 = Tree.edges(i,2);
                d2 = Tree.edges(i,1);
                d3 = sqrt((x-Tree.nodes(n2,1))^2+(y-Tree.nodes(n2,2))^2);
                theta = acosd((d1^2+d2^2-d3^2)/(2*d1*d2));
                check_edges = [check_edges;n2 d2 i theta];
            end
        end
        if size(check_edges,1)>0
            % Check which edge is closer by checking which one makes the smallest angle with the vector from (x,y) to (x_near,y_near).
            angles = check_edges(:,4);
            theta = min(angles);
            index = find(angles==theta);
            edge_nr = check_edges(index,3);

            % If the smallest angle is smaller than 90 degrees, the new edge
            % becomes a normal on the closest edge. If the angle is bigger than
            % 90 degrees, then the clostest node is a terminal node, and the
            % new edge becomes a "terminal edge", and (x,y) becomes the new
            % terminal node. 
            if theta < 90
                % Find point on edge
                d2 = check_edges(index,2);
                n2_koord = check_edges(index,1);
                ratio = d1*cosd(theta)/d2;
                x_new = x_near + ratio*(Tree.nodes(n2_koord,1)-x_near);
                y_new = y_near + ratio*(Tree.nodes(n2_koord,2)-y_near);

                % Find lenght of the new edges (a new node on an already existing edge will split up this edge in two)
                l1 = sqrt((x_new-Tree.nodes(Tree.edges(edge_nr,2),1))^2+(y_new-Tree.nodes(Tree.edges(edge_nr,2),2))^2);
                l2 = Tree.edges(edge_nr,1)-l1;

                % Update tree
                NewTree.edges = zeros(size(Tree.edges,1)+2,size(Tree.edges,2));
                NewTree.edges(1:edge_nr,:) = Tree.edges(1:edge_nr,:);
                NewTree.edges(edge_nr,:) = [l1 Tree.edges(edge_nr,2) size(Tree.nodes,1)+1 Tree.edges(edge_nr,4)];
                NewTree.edges(edge_nr+1,:) = [l2 size(Tree.nodes,1)+1 Tree.edges(edge_nr,3) Tree.edges(edge_nr,4)];
                NewTree.edges(edge_nr+2:end-1,:) = Tree.edges(edge_nr+1:end,:);
                NewTree.edges(end,:) = [d1*sind(theta) size(Tree.nodes,1)+1 size(Tree.nodes,1)+2 Tree.edges(edge_nr,4)*Tree.RadiusRate];
                NewTree.nodes = [Tree.nodes;x_new y_new;x y];
            else
                % Update tree
                NewTree.nodes = [Tree.nodes;x y];
                NewTree.edges = [Tree.edges;d1 k size(Tree.nodes,1)+1 Tree.edges(edge_nr,4)*Tree.RadiusRate];
            end
        elseif size(check_edges,1) == 0
            % The nearest node is a terminal node (there are edges in the
            % Tree, but not on this node, this will happen when building
            % the combinated tree).
            % Update tree
            NewTree.nodes = [Tree.nodes;x y];
            NewTree.edges = [Tree.edges;d1 k size(Tree.nodes,1)+1 Tree.TrunkRadius];
        else
            disp('Check edges feil')
        end
    else
        disp('Shady ting skjer 1')
    end
elseif size(Tree.edges,1) == 0
        % There are no edges at all in the Tree. 
        % Update tree
        NewTree.nodes = [Tree.nodes;x y];
        NewTree.edges = [Tree.edges;d1 k size(Tree.nodes,1)+1 Tree.TrunkRadius];
else
    disp('shady ting skjer 2')
end
NewTree.RadiusRate = Tree.RadiusRate;
end