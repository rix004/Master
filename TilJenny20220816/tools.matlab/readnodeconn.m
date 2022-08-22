% Fast version of updating connectivity
function [adjind, adjindsegment] = readnodeconn(nodeconn, nnode)
nsegment = size(nodeconn,1);
adjind = cell(nnode,1);
adjindsegment = cell(nnode,1);

% Generate one long vector of each connection
segnum = repmat((1:nsegment)', 2, 1);
% nc = nodeconn(:);
nc = [nodeconn;fliplr(nodeconn)];

% Sort according to increasing nodes
[sortval, sortind] = sort(nc(:,1));
segnum = segnum(sortind);
nc = nc(sortind,:);

% Find the unique values in nodes, to split the array at right places
[inters, ia, ~] = unique(sortval);
ia = [ia;2*nsegment+1];
split = ia(2:end) - ia(1:end-1);

% Assign the values to the cell arrays    
c = mat2cell(nc(:,2), split);

% Assign to cell array
adjind(inters) = c;

% Assign the values to the cell arrays    
c = mat2cell(segnum, split);
adjindsegment(inters) = c;

% return
% 
% % Old code
% adjind = cell(nnode,1);
% adjindsegment = cell(nnode,1);
% for i = 1 : 2
%     
%     ic = 3 - i;
%     
%     % Segment numbering
%     segnum = (1:nsegment)';
% 
%     % Sort according to increasing nodes
%     [sortval, sortind] = sort(nodeconn(:,i));
%     
%     % Also sort the segment number and connections
%     segnum = segnum(sortind);
%     conn = nodeconn(sortind, ic);
%     
%     % Find the unique values in nodes, to split the array at right places
%     [inters, ia, ~] = unique(sortval);
%     ia = [ia;nsegment+1];
%     split = ia(2:end) - ia(1:end-1);
%     
%     % Assign the values to the cell arrays    
%     c = mat2cell(conn, split);
%     adjindh = cell(nnode,1);
%     adjindh(inters, i) = c;
%     
%     c = mat2cell(segnum, split);
%     adjindsegmenth = cell(nnode,1);
%     adjindsegment(inters, i) = c;
%     
%     % Sum the cell arrays
%     adjind = cellfun(@(x1,x2) ([x1;x2]), adjind, adjindh, 'UniformOutput', false);
%     adjindsegment = cellfun(@(x1,x2) ([x1;x2]), adjindsegment, adjindsegmenth, 'UniformOutput', false);
% end
% 
% 
