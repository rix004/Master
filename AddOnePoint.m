function[NewTree]=AddOnePoint(Tree,x,y)
% Finner det punktet i Tree som ligger nærmest punktet (x,y). Dette punkter
% har nodenummer k.
[x_near,y_near,k,d1] = NearestNeighbour(x,y,Tree.nodes);

% Sjekker hvilke kanter (dersom kantmatrisen ikke er tom) som hører til den nærmeste noden, lagrer de i
% check_edges. Denne matrisen består av: 
% Kolonne 1: nodenummeret til den andre noden som danner kanten.
% Kolonne 2: Lengde på kant
% Kolonne 3: Kantnummer

if size(Tree.edges,1) > 0
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
    % Sjekker hvilken kant som er nærmest ved å sjekke hvilken kant som danner
    % minste vinkel med vektoren fra (x,y) til (x_near,y_near).
    angles = check_edges(:,4);
    theta = min(angles);
    index = find(angles==theta);
    edge_nr = check_edges(index,3);
    d2 = check_edges(index,2);
    n2_koord = check_edges(index,1);
    ratio = d1*cosd(theta)/d2;
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
        %line([x_new x],[y_new,y],'LineWidth',Tree.edges(edge_nr,4)*Tree.radius_ratio,'Color',[0.8500, 0.3250, 0.0980])
        %hold on
    else
        NewTree.nodes = [Tree.nodes;Tree.nodes(end,1)+1 x y Tree.nodes(Tree.edges(edge_nr,3),4)+1];
        NewTree.edges = [Tree.edges;d1 k Tree.nodes(end,1)+1 Tree.edges(edge_nr,4)];
        %line([x_near x],[y_near,y],'LineWidth',Tree.edges(edge_nr,4),'Color',[0.8500, 0.3250, 0.0980])
        %hold on
    end
else
        NewTree.nodes = [Tree.nodes;Tree.nodes(end,1)+1 x y Tree.nodes(end,4)+1];
        NewTree.edges = [Tree.edges;d1 k k+1 Tree.root_radius];
        %line([x_near x],[y_near,y],'LineWidth',Tree.root_radius,'Color',[0.8500, 0.3250, 0.0980])
        %hold on
end

end