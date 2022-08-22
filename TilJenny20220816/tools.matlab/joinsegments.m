function [tree, modifieds] = joinsegments(tree, indremove, modifieds, prm)

n = numel(indremove);

msg = ['Removing ' int2str(n) ' nodes in ' mfilename ' due to mid of a vessel'];
disp(msg);

validn = true(tree.node.n,1);
valids = true(tree.segment.n,1);
for i = 1 : n
    
    % Current node
    curnode = indremove(i);
    
    % Dont remove a root node
    if tree.node.isroot(curnode)
        continue;
    end
    
    % Neighbor nodes
    neighnode = tree.node.adjind{curnode};    
    conns = tree.node.adjindsegment{curnode};
    
    targetconn = conns(1);
    removeconn = conns(2);
    
    % Copy the segment and the node indices to another segment
    tree.segment.ind{targetconn} = [tree.segment.ind{targetconn}; ...
                                    tree.node.ind{curnode}; ...
                                    tree.segment.ind{removeconn}];
                                    
    % Copy the segment and node coordinates to another segment
    tree.segment.c{targetconn} = [tree.segment.c{targetconn}; ...
                                    tree.node.c{curnode}; ...
                                    tree.segment.c{removeconn}];
    % Update number of coordinates
    tree.segment.ncoord(targetconn) = numel(tree.segment.ind{targetconn});

    % Also copy distrind
%     if prm.tree.assignvoxels
%         tree.segment.distrind{targetconn} = [tree.segment.distrind{targetconn}; ...
%                                     tree.node.ind{curnode}; ...
%                                     tree.segment.distrind{removeconn}];
%     end
    
    % Update segment connections    
    tree.segment.nodeconn(targetconn,:) = neighnode;            
    
    % Update node connections neighbor 1
    bw = tree.node.adjind{neighnode(1)} == curnode;
    tree.node.adjind{neighnode(1)}(bw) = neighnode(2);
    tree.node.adjindsegment{neighnode(1)}(bw) = targetconn;
    
    % Update node connections neighbor 2
    bw = tree.node.adjind{neighnode(2)} == curnode;
    tree.node.adjind{neighnode(2)}(bw) = neighnode(1);
    tree.node.adjindsegment{neighnode(2)}(bw) = targetconn;

    % Central node
    tree.node.adjind{curnode}(1) = [];
    tree.node.adjindsegment{curnode}(1) = [];
    
    % Mark the nodes and segments as invalid
    validn(curnode) = false;
    valids(removeconn) = false;
    
    % Label the connection as updated
    modifieds(targetconn) = true;
    
end

% Crop the tree
tree = croptree(tree, prm, validn, valids, '', false);
modifieds = modifieds(valids);
% flag = verifytree(tree)
