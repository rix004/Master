function DrawTree(Tree)
figure(2)
for i = 1:size(Tree.edges,1)
    x1 = Tree.nodes(Tree.edges(i,2),1);
    x2 = Tree.nodes(Tree.edges(i,3),1);
    y1 = Tree.nodes(Tree.edges(i,2),2);
    y2 = Tree.nodes(Tree.edges(i,3),2);
    line([x1 x2],[y1,y2],'LineWidth',Tree.edges(i,4)*10,'Color',[0.8500, 0.3250, 0.0980])
    hold on
    %pause(0.1)
    axis([min(Tree.nodes(:,1)) max(Tree.nodes(:,1)) min(Tree.nodes(:,2)) max(Tree.nodes(:,2))])
end
axis equal
end