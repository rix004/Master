function[v]=FindNearestNodes(Tree)
% For alle noder, finn alle noder som er nærmere node i enn rotnoden, men
% som er lenger unna rotnoden enn node i er. 
% Output: Vektor v med like mange kolonner som noder, hver rad gir et tall på
% hvor mange andre noder node i "forsyner".
D = distances(graph(Tree.edges(:,2),Tree.edges(:,3)));
RootNodeDistances = D(Tree.RootNodeIdx,:);
v = zeros(size(Tree.nodes,1),1);
for i = 1:size(D,1)
    for j = 1:size(D,2)
        if D(i,j)<RootNodeDistances(j)  && RootNodeDistances(j) > RootNodeDistances(i) && 
            v(i)= v(i)+1;
        elseif i == Tree.RootNodeIdx
            v(i)=size(Tree.nodes,1)-1;
        end
    end
end