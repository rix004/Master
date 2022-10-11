function [cells, vertices,c] = VoronoiDiagram(Tn,bs_ext)
[v,c,xy]=VoronoiLimit(Tn(:,1),Tn(:,2),'bs_ext',bs_ext,'figure','off');
for i = 1:size(c,1)
    cell_index = find(ismember(xy,Tn(i,1:2),'rows'));
    if any(cell_index)
        cells(i)=c(cell_index(1));
    else
        disp('Terminalnode utenfor ytre grense');
    end
end
cells = cells';
vertices = v;