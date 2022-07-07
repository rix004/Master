function DrawTree(Tree,linewidth,color,D)
for i = 1:size(Tree.edges,1)
    x1 = Tree.nodes(Tree.edges(i,2),1);
    x2 = Tree.nodes(Tree.edges(i,3),1);
    y1 = Tree.nodes(Tree.edges(i,2),2);
    y2 = Tree.nodes(Tree.edges(i,3),2);
    line([x1 x2],[y1,y2],'LineWidth',Tree.edges(i,4)*linewidth,'Color',color)
    hold on
    %plot(Tree.nodes(:,1),Tree.nodes(:,2),'.','MarkerSize',10,'Color',color)
    axis(D)
end
axis equal
end