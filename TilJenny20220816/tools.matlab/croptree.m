function [tree] = croptree(tree, prm, validn, valids, labelterminal, genbw)
disp('Cropping tree');
node = tree.node;
segment = tree.segment;
dim = prm.dim;

if sum(~validn) == 0
    if sum(~valids) == 0
        disp('No cropping, terminating');
        return;
    end
end


if ~prm.tree.savememory
    root.bw = false(dim);
    terminal.bw = false(dim);
    segment.faser = zeros(prm.dim); 
    node.faser = zeros(dim);
end

if ~isfield(prm, 'verbose')
    prm.verbose = false;
end

% Number of nodes
node.n = sum(validn);
node.num = (1:node.n)';

% Number of segments
segment.n = sum(valids);
segment.num = segment.numfnc((1:segment.n)');

% Transformations for nodes
full2redn = NaN(tree.node.n,1);
full2redn(validn) = node.num;

%----------------------------
% Update segments in the tree
%----------------------------

% Update node connection
nodeconn = segment.nodeconn(valids,:);
nodeconn(:,1) = full2redn(nodeconn(:,1));
nodeconn(:,2) = full2redn(nodeconn(:,2));
segment.nodeconn = nodeconn;
clear nodeconn;


% Segment fields with no ordering
fn = {'ind', 'c', 'dist', 'L', 'avrad', 'ncoord', 'geodesdist', ...
    'distrind', 'phi', 'distrvol', 'distrn', 'distrind', 'kN','qN'};
for i = 1 : numel(fn)
    fnh = fn{i};
    if prm.verbose
        msg = ['Restricting field ' fnh];
        disp(msg);
    end
    try
        segment.(fnh) = segment.(fnh)(valids);
    catch
        if prm.verbose
            msg = ['Could not restrict field ' fnh];
            disp(msg);
        end
    end
end

% % Generate images
% if ~prm.tree.savememory
%     segment.bw = false(dim);
%     for i = 1 : segment.n
%         % Skeleton and faser 
%         ind = segment.ind{i};   
%         segment.bw(ind) = true;
%         segment.faser(ind) = i;        
%     end    
% end

% try
%     for i = 1 : segment.n
%         % Skeleton and faser 
%         ind = cell2mat(segment.distrind{i});
%         bw(ind) = true;
%     end    
% catch
%     disp('Could not restrict distrind');
% end

% Update summarizing fields in struct
segment.nvar = nnz(cell2mat(segment.ind));

%----------------------------
% Update nodes in the tree
%----------------------------

% Make adjacency matrices of connectivity
fval = ones(segment.n,1);
node.adjmat = adjmat(segment.nodeconn,node.n,fval);

% Read node connections
[node.adjind, node.adjindsegment] = readnodeconn(segment.nodeconn, node.n);

% Update fields with no ordering
fn = {'ind', 'c', 'isterminal', ...
    'isroot', 'avc', 'avind', ...
    'ncoord', 'avx', 'nodelevel', ...
    'isinterior', 'distrind', 'phi', ...
    'distrvol', 'distrn', 'avctrue', ...
    'belongroot', 'shp', 'kT', ...
    'cumgraphvar', 'graphvar', ...
    'diam', 'graphnum', 'pN', ...
    'group', 'distrind','shpnum'};
for i = 1 : numel(fn)
    fnh = fn{i};
    try        
        node.(fnh) = node.(fnh)(validn,:);
    catch
        if prm.verbose
            msg = ['Could not restrict field ' fnh];
            disp(msg);
        end
    end
end
% Special treatment for splitnum to remove spaces
try
    [~,I] = sort(node.graphnum);
    [~,J] = sort(I);
    s = (1:node.n)';
    node.graphnum = s(J);
    clear I J s
catch
    msg = ['Could not update order of graphnum'];
    disp(msg);
end

if ~isempty(labelterminal)
    % Update terminal and interior labeling based on new connectivity
    % Special case of one node only, then its both a root and a terminal
    if node.n == 1
        node.isterminal = true;
        node.isroot = true;
        node.isinterior = false;
    else
        node.(labelterminal) = full(sum(node.adjmat,2) == 1 | sum(node.adjmat,2) == 0);
        node.isinterior = full(sum(node.adjmat,2) > 1);           
    end
end
node.isinterior(node.isroot) = false;

% Relabel nodes
node = relabelnodes(node);
  
% % Generate images
% if ~prm.tree.savememory
%     node.bw = false(dim);
%     for i = 1 : node.n    
%         ind = node.ind{i};    
%         node.bw(ind) = true;    
%         node.faser(ind) = i; 
%     end    
% end

% if prm.tree.assignvoxels
%     try
%         tree.bw = false(prm.dim);
%         ind = cell2mat(segment.distrind); 
%         tree.bw(ind) = true;
%         ind = cell2mat(node.ind); 
%         tree.bw(ind) = true;
%         clear ind;
%     catch
%        disp('Could not assign voxels to bw'); 
%     end
% end


% if ~prm.tree.savememory
%     for i = 1 : node.nterminal
%         ind = node.ind{node.indterminal(i)};
%         terminal.bw(ind) = true; 
%     end
%     tree.terminal = terminal;clear terminal
%     for i = 1 : node.nroot
%         ind = node.ind{node.indroot(i)};
%         root.bw(ind) = true; 
%     end
%     tree.root = root;clear root;
%     tree.skel = segment.faser > 0 | node.faser > 0;
% end

% Assign
tree.node = node;clear node
tree.segment = segment;clear segment

% Regenerate tree
% try
%     
% catch
%    disp('Could not delete bw in tree'); 
% end
% tree = rmfield(tree, 'bw');
if genbw    
    tree.bw = tree2bw(tree, prm); 
end


