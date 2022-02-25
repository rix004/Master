function[NewTree]=BinaryTree(Tree)
if size(Tree.edges,1) > 0
    StartNodes = Tree.edges(:,2);
    EndNodes = Tree.edges(:,3);
    Tn = []; Te = [];
    for i = 1:length(EndNodes)
        v = ismember(EndNodes,StartNodes);
        if v(i) == 0
            Tn = [Tn;EndNodes(i)];
            Te = [Te;i];
        end
    end
    for i = 1:length(Tn)
        x0 = Tree.nodes(Tn(i),2);
        y0 = Tree.nodes(Tn(i),3);
        x1 = x0+cosd(Tree.RotAngle)*Tree.edges(i,1)*Tree.LengthRate;
        y1 = y0+sind(Tree.RotAngle)*Tree.edges(i,1)*Tree.LengthRate;
        x2 = x0+cosd(Tree.RotAngle)*Tree.edges(i,1)*Tree.LengthRate;
        y2 = y0+sind(Tree.RotAngle)*Tree.edges(i,1)*Tree.LengthRate;
        NewTree.nodes=[Tree.nodes;
                        Tree.nodes(end,1)+1,x1,y1,Tree.nodes(Tn(i),4)+1;
                        Tree.nodes(end,1)+2,x2,y2,Tree.nodes(Tn(i),4)+1]; % Plusser på nodenummer og nivå
        NewTree.edges=[Tree.edges;
                        sqrt((x1-x0)^2+(y1-y0)^2) Tn(i) Tree.nodes(end,1)+1 Tree.edges(Te(i),4)*Tree.RadiRate;
                        sqrt((x2-x0)^2+(y2-y0)^2) Tn(i) Tree.nodes(end,1)+2 Tree.edges(Te(i),4)*Tree.RadiRate;];
    Tree.nodes = NewTree.nodes;
    Tree.edges = NewTree.edges; 
    end
    
else
    x0 = Tree.nodes(1,2);
    y0 = Tree.nodes(1,3);
    x1 = x0+cosd(Tree.RotAngle)*Tree.TrunkLength;
    y1 = y0+sind(Tree.RotAngle)*Tree.TrunkLength;
    x2 = x0-cosd(Tree.RotAngle)*Tree.TrunkLength;
    y2 = y0+sind(Tree.RotAngle)*Tree.TrunkLength;
    NewTree.nodes=[Tree.nodes;
                    Tree.nodes(1,1)+1,x1,y1,Tree.nodes(1,4)+1;
                    Tree.nodes(1,1)+2,x2,y2,Tree.nodes(1,4)+1];
    NewTree.edges=[Tree.edges;
                    sqrt((x1-x0)^2+(y1-y0)^2) Tree.nodes(1,1) Tree.nodes(1,1)+1 Tree.TrunkRadi;
                    sqrt((x2-x0)^2+(y2-y0)^2) Tree.nodes(1,1) Tree.nodes(1,1)+2 Tree.TrunkRadi;];
                
    Tree.nodes = NewTree.nodes;
    Tree.edges = NewTree.edges;
end
NewTree.nodes = Tree.nodes;
NewTree.edges = Tree.edges;
end