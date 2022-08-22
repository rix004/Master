function [valids] = genvalidsegments(adjindsegment, n, validn)

% Valid segments   
valids = true(n,1);
adjindsegment = unique(cell2mat(adjindsegment(~validn)));
valids(adjindsegment) = false;
% ind = find(~validn);
% for j = 1 : numel(ind)
%     notvalid = adjindsegment{ind(j)};
%     valids(notvalid) = false;
% end
