function DrawTree(Tree)
figure(2)
for i = 1:size(Tree.edges,1)
    x1 = Tree.nodes(Tree.edges(i,2),2);
    x2 = Tree.nodes(Tree.edges(i,3),2);
    y1 = Tree.nodes(Tree.edges(i,2),3);
    y2 = Tree.nodes(Tree.edges(i,3),3);
    line([x1 x2],[y1,y2],'LineWidth',Tree.edges(i,4),'Color',[0.8500, 0.3250, 0.0980])
    hold on
    %pause(0.1)
    axis([min(Tree.nodes(:,2)) max(Tree.nodes(:,2)) min(Tree.nodes(:,3)) max(Tree.nodes(:,3))])
end
axis equal
end