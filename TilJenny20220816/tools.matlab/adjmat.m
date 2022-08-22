% Make adjacency matrix
function [adj] = adjmat(nodeconn,nnode,fval)

if isempty(nodeconn)
    adj = zeros(nnode,nnode);
    return;
end
fval = fval(:);
% nodeconn = cell2mat(nodeconn);
nodeconn = [nodeconn; fliplr(nodeconn)];
adj = sparse(nodeconn(:,1), nodeconn(:,2), repmat(fval, 2, 1), nnode,nnode);
