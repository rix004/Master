function [node] = relabelnodes(node)

node.isinterior = full(sum(node.adjmat, 2) > 1);

node.indinterior = find(node.isinterior);
node.indterminal = find(node.isterminal);
node.indroot = find(node.isroot);

node.ninterior = numel(node.indinterior);
node.nterminal = numel(node.indterminal);
node.nroot = numel(node.indroot);
