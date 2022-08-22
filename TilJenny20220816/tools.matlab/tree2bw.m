function [bw] = tree2bw(tree, prm)

bw = false(prm.dim);
bw(cell2mat(tree.segment.ind)) = true;
try
    bw(cell2mat(tree.segment.distrind)) = true;
catch
    msg = ['Could not assign distrind in ' mfilename];
    disp(msg)   
end
bw(cell2mat(tree.node.ind)) = true;
