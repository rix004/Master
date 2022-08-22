function [tree] = simfullbraintreefun(varargin)

% Default parameters
prm.aiftype = '';
prm.showaif = false;
prm.setaif = false;

bwtree = varargin{1};
phi = varargin{2};
bwbrain = varargin{3};
treepar = varargin{4};
prmin = varargin{5};
if ~isempty(prm)
    prm = mergeinputpar(prm, prmin); 
end
clear varargin;

dim = size(bwtree);
if numel(dim) == 2
    dim = [dim,1];
end
% % Assign tree
% tree.bw = bw;

if dim(3) == 1
    prm.conn = 8;
else
   prm.conn =  26;
end

% Datatyper
prm.indfnc = getintfnc(prod(prm.dim));
prm.coordfnc = getintfnc(max(prm.dim));

%--------------------------------------------------------
% Generate skeleton
%--------------------------------------------------------
[bwskeleton, bwnode, bwterminal, bwtree] = simfullbrainskel3DMatlab(bwtree, prm);

% Ensure nothing goes outside the mask
bwtree(bwskeleton) = true;

% Save data
tree.bw = bwtree; clear bwtree;
tree.skel.bw = bwskeleton; clear bwskeleton;
tree.node.bw = bwnode; clear bwnode
tree.terminal.bw = bwterminal; clear bwterminal;

%--------------------------------------------------------
% Find root points
%--------------------------------------------------------
disp('Find AIF points')
tic
% Convert to matrix coordinates
aifpos = real2coord(treepar.aifpos, prm.h, prm.dim, true);

% Get AIF positions and generate root matrix
[tree.root.bw, aifpos] = ...
            getaifpos(bwbrain, tree.terminal.bw | tree.node.bw, aifpos, treepar.aifoption); 
clear bwbrain;
disp(['Number of root points: ' int2str(size(aifpos, 1))])
clear faser;
toc

if prm.showaif
    H = figure;plotisosurface(H,imdilate(root.bw,ones(2,2,2)),[1,1,1],'blue',0.9,0); grid on
    hold on;plotisosurface(H,bwtree,[1,1,1],'red',0.9,0); grid on
    pause
end
    
if prm.setaif
    H = figure;plotisosurface(H,imdilate(terminal.bw,ones(2,2,2)),[1,1,1],'blue',0.9,0); grid on
    hold on;plotisosurface(H,bwtree,[1,1,1],'red',0.9,0); grid on
    pause
end

% ----------------------------------------------------------------
% ----- Only keep objects connected spatially to root points -----
% ----------------------------------------------------------------
msg = ['Only keep objects that are spatially connected to any root'];
disp(msg);
faser = bwconncomp(tree.bw, prm.conn);
faser.PixelIdxList = (cellfun(prm.indfnc,faser.PixelIdxList,'UniformOutput', false))';       
indroot = find(tree.root.bw);%prm.indfnc(cell2mat(tree.node.indroot));
isvalid = false(faser.NumObjects,1);
for i = 1 : faser.NumObjects
    isvalid(i) = any(ismember(indroot, faser.PixelIdxList{i}));        
end
msg = ['Removing ' int2str(faser.NumObjects - sum(isvalid)) ' objects'];
disp(msg);

if sum(isvalid) < faser.NumObjects    
    disp('Regenerate tree skeleton and update binary tree');
    tic
    faser.PixelIdxList = faser.PixelIdxList(isvalid);
    faser.NumObjects = sum(isvalid);
    tree.bw = false(dim);
    tree.bw(cell2mat(faser.PixelIdxList)) = true;
    % Update binary masks    
    tree.skel.bw = tree.skel.bw & tree.bw;
    tree.node.bw = tree.node.bw & tree.bw;
    tree.terminal.bw = tree.terminal.bw & tree.bw;
    toc
end
clear faser isvalid;

% This can happen if the tree is only one point
if nnz(tree.node.bw) == 0 && nnz(tree.terminal.bw) == 0   
   tree.terminal.bw = tree.skel.bw;
end 

% Nodes are branching points, tips and roots
tree.node.bw = tree.node.bw | tree.terminal.bw | tree.root.bw;

% Generate node info
tree.node = nodefun(tree.node, tree.root, tree.terminal, prm);

% Generate segment info
tree.segment.bw = tree.skel.bw;
tree.segment.bw(tree.node.bw) = false;
tree.node = rmfield(tree.node, 'bw');
if prm.tree.savememory
    tree = rmfield(tree, 'root');
    tree = rmfield(tree, 'terminal');
    tree = rmfield(tree, 'skel');
end
tree.segment = segmentfun(tree.segment, tree.node, tree.bw, prm);

% Make adjacency matrices of connectivity
fval = ones(tree.segment.n,1);
tree.node.adjmat = adjmat(tree.segment.nodeconn,tree.node.n,fval);clear fval;

% Read node connections
[tree.node.adjind, tree.node.adjindsegment] = readnodeconn(tree.segment.nodeconn, tree.node.n);
    
% % Make sure that segments are sorted according to a user-defined input, if
% % given
% nodeconn = [];
% if ~isempty(nodeconn)        
%     if tree.segment.n ~= size(nodeconn,1)
%         error('Wrong number of segments');
%     end    
%     % Find the match of segment-endpoints to identify the segments
%     match = false(segment.n, segment.n);
%     for i = 1 : segment.n
%         for j = 1 : size(nodeconn,1)
%             if isequal(segment.nodeconn(i,:), nodeconn(j,:)) || ...
%                     isequal(segment.nodeconn(i,:), flip(nodeconn(j,:)))
%                 match(i,j) = true;
%             end
%         end        
%     end
%     % Find the new order of segments    
%     order = zeros(segment.n,1);
%     for i = 1 : segment.n
%         order(i) = find(match(:,i) > 0);
%     end
%     % Reassign faser of segments
%     faser = zeros(dim);
%     for i = 1 : segment.n
%        reghere = segment.faser == order(i);
%        faser(reghere) = i;              
%     end
%     segment.faser = faser;
%     % Reassign other fields
%     segment.nodeconn = segment.nodeconn(order,:);
%     segment.ind = segment.ind(order);
%     segment.c = segment.c(order);
% end                                                       
                                                                 
% --------------------------------------------------------------------
% Assign correct root belonging. Must find matching coordinates of the
% root points since they might not be exactly equal
% --------------------------------------------------------------------
c0 = aifpos;clear aifpos;
c1 = double(tree.node.avc(tree.node.indroot, :));
d = pdist2(c1, c0);
[assign,~] = munkres(d);
tree.node.belongroot = repmat({'none'}, tree.node.n, 1);
tree.node.belongroot(tree.node.indroot) = prm.paths.belongroot(assign)';

%---------------------------------------------------------------
% Assign voxels to segments 
%---------------------------------------------------------------
% (NB must be done before we remove any edges, or we will get voxels not 
% belonging to the right edges)
if prm.tree.assignvoxels
    msg = ['Assigning voxels to skeleton'];
    disp(msg);
    tic
    tree.segment = assignvoxels(tree.bw, tree.node, tree.segment, prm);
    tree.bw = tree2bw(tree, prm);
    toc;
end        
clear faser indstart indstop c;

%---------------------------------------------------------------
% Remove nodes that have only 2 connections - these are bumps in the road
%---------------------------------------------------------------

% Label the modified segments - to be used later
modifieds = false(tree.segment.n, 1);

% Number of connections
nconn = full(sum(tree.node.adjmat, 2));

% These are the nodes that are midlines nodes and that we want to remove
indremove = find(nconn == 2 & ~tree.node.isroot);

% Remove the given nodes
[tree, modifieds] = joinsegments(tree, indremove, modifieds, prm);

% Reorder coordinates
tree.segment = reordercoord(tree.segment, modifieds);

% Update vessel length and geodesic distance
tree.segment = compvessellength(tree.segment, tree.node, prm, modifieds);

% ----------------------------------------------------                            
% -------------- Remove self-loops -------------------
% ----------------------------------------------------
% NB This must not be done after assigning voxels, then we will end up with
% removing parts of the bniary tree
[G, ~, ~] = makegraph(tree.segment.nodeconn, tree.segment.n, tree.node.n);
G = rmedge(G, 1:numnodes(G), 1:numnodes(G));
validn = true(tree.node.n,1);
valids = false(tree.segment.n,1);
valids(G.Edges.index) = true;
tree = croptree(tree, prm, validn, valids, '', false);    

%----------------------------------------------------------------------
% Remove edges that are shorter than wide, and which have a large radius
%----------------------------------------------------------------------
% Doesnt this code produce disconnected parts of the tree?
% valids = tree.segment.L./tree.segment.avrad > 0.1 & tree.segment.avrad < 12*mean(prm.h);
% % Dont remove any edges connected to a root
% valids(cell2mat(tree.node.adjindsegment(tree.node.indroot))) = true;
% validn = true(tree.node.n,1);
% msg = ['Removing ' int2str(nnz(valids==0)) ' segments due to suspected shortcuts'];
% disp(msg);
% tree = croptree(tree, prm, validn, valids, '', false);

% -----------------------------------------------------
% Check that all nodes have a connection to a root node
% -----------------------------------------------------

% Find shortest path to root
msg = 'Find shortest path to any root';
disp(msg);
[G, ~, ~] = makegraph(tree.segment.nodeconn, tree.segment.n, tree.node.n);
tic
nodelevel = distances(G, tree.node.indroot)';
toc

% Valid nodes
validn = sum(isinf(nodelevel), 2) < tree.node.nroot;
% Crop away invalid nodes
if any(validn == 0)
    % Valid segments
    valids = genvalidsegments(tree.node.adjindsegment, tree.segment.n, validn);

    % Crop tree
    tree = croptree(tree, prm, validn, valids, [], true);
end

%---------------------------------------------------------------
% Extract porosity of individual nodes and segments?
%---------------------------------------------------------------
if ~isempty(phi)
    [tree.node,tree.segment] = assignphi(tree.node, tree.segment, phi);
    [tree.node,tree.segment] = assigndistrvol(tree.node,tree.segment,prm.h);    
end

% % Check for node position
% if ~isempty(prm.tree.nodepos)
%     tree = reassignposnode(tree, prm.tree);
% end

end

%------------------------------------------------

function [node] = nodefun(node, root, terminal, prm)        

%--------------------------------------------
% Connected components of nodes
%--------------------------------------------
disp('Generate connected components of nodes');
faser = bwconncomp(node.bw, prm.conn);
faser.PixelIdxList = (cellfun(@prm.indfnc,faser.PixelIdxList,...
    'UniformOutput', false))';
node.ind = faser.PixelIdxList;
node.n = numel(faser.PixelIdxList);
clear faser;
node.numfnc = getintfnc(node.n);

% Node numbering
node.num = node.numfnc((1:node.n)');

% If only one root, its also a terminal
if node.n == 1    
    node.isterminal = node.isroot;    
end

%--------------------------------------------
% Generate node coordinates and indices
%--------------------------------------------
disp('Assign nodes coordinates');
tic
c = cell(node.n,1);
avc = prm.coordfnc(zeros(node.n,3));
avind = prm.indfnc(zeros(node.n,1));
isterminal = false(node.n,1);
isroot = false(node.n,1);
ncoord = NaN(node.n,1);
for i = 1 : node.n
    
    % This node   
    indh = node.ind{i};
         
    % Coordinate
    n = numel(indh);
    ch = prm.coordfnc(zeros(n,3));   
    [ch(:,1), ch(:,2), ch(:,3)] = ind2sub(prm.dim, indh);
    c{i} = ch;
   
    % Number of coordinates in node
    node.ncoord(i) = numel(indh);
    
    % Average coordinate   
    mid = round(numel(indh)/2);
    avc(i,:) = c{i}(mid,:);
    avind(i) = sub2ind(prm.dim,avc(i,1),avc(i,2),avc(i,3));         
      
    
    if sum(root.bw(indh)) > 0
        isroot(i) = true;
    end
    if sum(terminal.bw(indh)) > 0
        isterminal(i) = true;
    end   
end
node.ncoord = ncoord;clear ncoord;
node.isroot = isroot;clear isroot;
node.isterminal = isterminal;clear isterminal;
node.c = c;clear c;
node.avc = avc;clear avc;
node.avind = avind;clear avind;
clear c;
toc

% Indexes    
if ~isfield(node, 'avx')
    % Compute average matrix coordinates of nodes
    node.avx = coord2real(node.avc, prm.h);
end

% Interior nodes 
node.isinterior = imcomplement(node.isterminal | node.isroot);    
   
% Linear, local indexes
node.indterminal = find(node.isterminal);
node.indroot = find(node.isroot);
node.indinterior = find(node.isinterior);

% Number of terminals and nodes
node.nterminal = nnz(node.isterminal);
node.nroot = nnz(node.isroot);
node.ninterior = nnz(node.isinterior);

% Node numbering
node.num = (1:node.n)';

% Print properties
msg = ['There are ' int2str(node.n) ' nodes'];
disp(msg);

msg = ['There are ' int2str(node.nroot) ' root nodes'];
disp(msg);

msg = ['There are ' int2str(node.nterminal) ' terminal nodes'];
disp(msg);

msg = ['There are ' int2str(node.ninterior) ' interior nodes'];
disp(msg);
    
end

%------------------------------------------------

function [node, segment] = assignphi(node, segment, im)

eps = 0.05;
% Assign to nodes
node.phi = eps*ones(node.n,1);
for i = 1 : node.n
   ind = node.ind{i};
   val = mean(im(ind));
%    val = 1;
   node.phi(i) = max(val,eps);     
end
% node.phi(1:2) = 1;
% node.phi(4) = 1;
% node.phi(:) = 1;

% node.phi(node.phi == 0) = eps;

% Assign to segments
segment.phi = cell(segment.n,1);
for i = 1 : segment.n
    segment.phi{i} = NaN(segment.ncoord(i),1);
    for j = 1 : segment.ncoord(i) 
        ind = segment.distrind{i}{j};
        val = mean(im(ind));  
%         % Take away
%         val = ones(size(val));
        segment.phi{i}(j) = max(val,eps);         
    end
    % In case the segment has zero voxels...
    segment.phi{i}(isnan(segment.phi{i})) = eps;
end
end

%------------------------------------------------

function [segment] = assignvoxels5(bw, node, segment, prm)
fnc = getintfnc(segment.n);
labelim = fnc(false(prm.dim));
for i = 1 : segment.n
   labelim(segment.ind{i}) = i;
end

a = 2;
end
%------------------------------------------------


%------------------------------------------------

function [segment] = assignvoxels3(bw, segment, prm)

% Generate image for segmentation
fnc = getintfnc(segment.n);

% Generate markers
skel = false(prm.dim);
skel(cell2mat(segment.ind)) = true;

% Impose markers
wat = fnc(imcomplement(bw));
wat = imimposemin(wat, skel);clear skel;

% Do segmentation
wat = vincent_soille_watershed(wat, 8);
wat(bw == 0) = 0;clear bw;

% Generate mapping
mapping = zeros(segment.n,1);
for i = 1 : segment.n    
    mapping(i) = mode(wat(segment.ind{i}));
end

% Extract values to a cell array and apply mapping to get the right order

segment.distrind = wat(mapping)';

end

%------------------------------------------------

function [node,segment] = assignvoxels0(bw, node,segment,prm)

nbw = nnz(bw);
bwc = imcomplement(bw);clear bw;

% Generate label image of the skeleton
fnc = getintfnc(segment.n);
faser = fnc(false(prm.dim));
for i = 1 : segment.n
   ind = segment.ind{i};
   faser(ind) = i; 
end

% valmax = segment.n + 1;
ind0 = find(faser > 0);
val0 = faser(ind0);
c = 0;
while 1
    c = c + 1;
    imtemp = fnc(zeros(prm.dim));
    for i = -1:1
       for j = -1:1
           for k = -1:1
               % Move image
               imtemp = max(imtemp, transim(faser,i,j,k));
           end
       end
    end
    % Remove background
    imtemp(bwc) = 0;

    % Remore previous values of faser
    imtemp(faser > 0) = 0;
        
    % Find the assigned values
    bwh = imtemp > 0;
    
    % Assign them to faser
    faser(bwh) = imtemp(bwh);clear bwh;
    
    n = nnz(faser > 0);
    msg = ['Dilation ' int2str(c) ', number of pixels assigned: ' int2str(n) ' (out of ' int2str(nbw) ')'];
    disp(msg);        

    % We have covered all values in bw
    if n == nbw
        break;
    end    

end
faser(ind0) = val0;

% clear faser assingedbw indbw;
indglobal = prm.indfnc(find(faser>0));
indlabel = faser(indglobal); 

% Sort the array to count number of occurrences
[indlabel, indsort] = sort(indlabel);
indglobal = indglobal(indsort);

% Count number of occurrence of each value
uval = unique(indlabel);
n = numel(uval);
N = histcounts(indlabel, n);
segment.distrind = mat2cell(indglobal, N);


end
%------------------------------------------------
% Old code
function [node,segment] = assignvoxels2(bw,node,segment,prm)
        
    indbw = find(bw);
    nbw = numel(indbw);    
    cbw = single(zeros(nbw,3));
    [cbw(:,1),cbw(:,2),cbw(:,3)] = ind2sub(prm.dim, indbw);

    
    % Coordinates of nodes    
    labn = cell(node.n,1);
    ncoord = node.ncoord;
    num = node.num;
    parfor i = 1 : node.n
        nc = ncoord(i);
        labn{i,1} = prm.indfnc(num(i)*ones(nc,1));
    end
    cn = single(cell2mat(node.c));
    labn = cell2mat(labn);
    
    % Coordinates of segments    
    labsegment = cell(segment.n,1);
    labnumber = cell(segment.n,1);
    
    ncoord = segment.ncoord;
    parfor i = 1 : segment.n 
        nc = ncoord(i);       
        labsegment{i,1} = prm.indfnc(i*ones(nc,1));
        labnumber{i,1} = prm.indfnc((1:nc)');
    end
    cs = single(cell2mat(segment.c));
    labs = cell2mat(labsegment);
    numbers = cell2mat(labnumber);   
    
    % Initialize
    node.distrind = cell(node.n,1);
    segment.distrind = cell(segment.n,1); 
    ncoord = segment.ncoord;
    for i = 1 : segment.n   
       nc = ncoord(i);
       segment.distrind{i} = cell(nc,1);
    end

    % ----------------- Check each voxel for belonging--------------
    
    for i = 1 : nbw

        ch = single(cbw(i,:));
        ind = indbw(i);

        % Distance to nodes        
        dist = ch - cn;
        dist = sqrt(sum(dist.^2,2));
        [distn,idxn] = min(dist);        

        % Distance to segment
        if ~isempty(cs)
            dist = bsxfun(@minus,ch,cs);    
            dist = sqrt(sum(dist.^2,2));
            [dists,idxs] = min(dist);
        else
            dists = Inf;
        end
        
        % Are we closest to node or segment?
        if distn <= dists    
            % Close to node        
            lab = labn(idxn);
            node.distrind{lab} = [node.distrind{lab};ind];
        elseif distn > dists
            % Close to segment            
            lab = labs(idxs);
            number = numbers(idxs);
            segment.distrind{lab}{number} = [segment.distrind{lab}{number};ind];
        end   
    end
end

%-----------------------------------------

function [node,segment] = assigndistrvol(node,segment,h)

    % Voxel volume
    voxelvol = prod(h);
    
    % Measure distribution volume
    node.distrvol = zeros(node.n,1);
    node.distrn = zeros(node.n,1);
    for i = 1 : node.n
        node.distrn(i) = max(numel(node.distrind{i}),1);
        node.distrvol(i) = node.distrn(i)*voxelvol;                 
    end
    segment.distrvol = cell(segment.n,1);    
    segment.distrn = cell(segment.n,1);
    for i = 1 : segment.n
        ncoord = segment.ncoord(i);
        segment.distrvol{i} = zeros(ncoord,1);
        segment.distrn{i} = zeros(ncoord,1);
        for j = 1 : ncoord            
            segment.distrn{i}(j) = max(numel(segment.distrind{i}{j}),1);
            segment.distrvol{i}(j) = segment.distrn{i}(j)*voxelvol; 
        end
    end
        
end

%------------------------------------------

function [segment] = segmentfun(segment, node, bwtree, prm)
    
%-----------------------------------------
% Assign segments
%-----------------------------------------
disp('Generating connected components of vessels');
tic
faser = bwconncomp(segment.bw, prm.conn);
segment = rmfield(segment, 'bw');
toc

%-----------------------------------------
% Find the indices of segments
%-----------------------------------------
disp('Making segment indices');
tic
segment.ind = faser.PixelIdxList';
clear faser;
segment.n = numel(segment.ind);
segment.numfnc = getintfnc(segment.n);
segment.c = cell(segment.n,1);
segment.ncoord = zeros(segment.n,1);
coordfnc = prm.coordfnc;
indfnc = prm.indfnc;
for i = 1 : segment.n
    segment.ind{i} = indfnc(segment.ind{i});
    n = numel(segment.ind{i});
    ch = coordfnc(zeros(n,3));
    [ch(:,1), ch(:,2), ch(:,3)] = ind2sub(prm.dim, segment.ind{i});
    segment.c{i} = ch;    
    segment.ncoord(i) = numel(segment.ind{i});
end
toc

%-----------------------------------------
% Order coordinates in consecutive order
%-----------------------------------------
modifieds = true(segment.n,1);
segment = reordercoord(segment, modifieds);clear modifieds;

% -------------------------------------------------
% --- Find node connections using sliding window --
% -------------------------------------------------
disp('Find node connections using sliding window');
tic
% Generate node label image
node.faser = zeros(prm.dim, char(node.numfnc));
ind = node.ind;
for i = 1 : node.n
    node.faser(ind{i}) = i;
end
clear ind;

% Extract start and end point of the vessels
start = prm.indfnc(zeros(segment.n, 3));
stop = prm.indfnc(zeros(segment.n, 3));
for i = 1 : segment.n
    start(i,:) = single(segment.c{i}(1,:));
    stop(i,:) = single(segment.c{i}(end,:));
end

% Make a moving window
indstart = sub2ind(prm.dim, start(:,1), start(:,2), start(:,3));
indstop = sub2ind(prm.dim, stop(:,1), stop(:,2), stop(:,3));
% clear start stop;
vals = [-1,0,1];
nc = zeros(segment.n,2);
idx{1} = 1:prm.dim(1);
idx{2} = 1:prm.dim(2);
idx{3} = 1:prm.dim(3);
for i = vals    
    for j = vals        
        for k = vals            
            % Move the label image
            idxh{1} = idx{1} + i;
            idxh{2} = idx{2} + j;
            idxh{3} = idx{3} + k;            
            idxh = limitidx(idxh, prm.dim);
            % This requires less RAM than transim
            % faser = transim(node.faser, i, j, k);            
            faser = node.faser(idxh{1}, idxh{2}, idxh{3});            
            
            % Read out node-intersection with the segment end
            val1 = double(faser(indstart));        
            val2 = double(faser(indstop));
        
            % Ensure that we dont assign the same node to both the starting and
            % the end point! This means that we cut the graph and get multiple
            % connections to the same node. It happended with the old code for
            % segments with length 1 pixel.
            
            % Valid nodes for assignment must
            % 1. Be different from the other node 
            % 2. Have a zero value (unassigned)
            valid = (val1 ~= nc(:,2)) & (nc(:,1) == 0);
            nc(valid,1) = val1(valid);

            valid = (val2 ~= nc(:,1)) & (nc(:,2) == 0);
            nc(valid,2) = val2(valid);
        end
    end
end
node = rmfield(node, 'faser');

%-----------------------------------------------------------------------
% ---- Remove segments with less than two identified connections -------
%-----------------------------------------------------------------------
notvalid = any(nc==0, 2);
msg = ['Number of zero valued connections: ' int2str(sum(notvalid))];    
disp(msg);
toc
segment.nodeconn = nc;clear nc;

if sum(notvalid) > 0
    segment.ind(notvalid) = [];
    segment.c(notvalid) = [];
    segment.ncoord(notvalid) = [];
    segment.n = numel(segment.ind);
    segment.nodeconn(notvalid,:) = [];
end

%     %-----------------------------------------
%     % Read out long-distance node connections
%     %-----------------------------------------
% 
%     tic
%     disp('Read out node connections');
%     tic;    
%     avc = single(node.avc);
%     for i = 1 : segment.n
%         if ~notvalid(i)
%             continue;
%         end
%         % Find node connections
%         dist = sqrt(sum((start(i,:) - avc).^2,2)); 
%         [~, indmin1] = min(dist);     
%         dist = sqrt(sum((stop(i,:) - avc).^2,2)); 
% 
%         % Ensure that we dont connect the segment to the same node twice    
%         dist(indmin1) = Inf;
%         [~, indmin2] = min(dist);
% 
%         % Assign
%         nc(i,:) = [indmin1, indmin2];         
% 
%     end
%     toc

%-----------------------------------------------------
% Computing vessel length and geodesic distane
%-----------------------------------------------------
msg = 'Computing vessel length and geodesic distance';
disp(msg);
update = true(segment.n,1);
segment = compvessellength(segment, node, prm, update);

% ---------------------------------------------------
% ------- Computing vessel radius--------------------
% ---------------------------------------------------

msg = ['Computing vessel radius by distance map'];
disp(msg);

% Make a cell image to save RAM for large datasets

% Crop size
d0 = 400;

% Mean pixel size as an approximation
meanh = mean(prm.h);

if prm.dim(3) == 1
    radim = prm.coordfnc(round(bwdist(imcomplement(bwtree))));
else
    % Find a slicing for patches
    v = cell(prm.ndim,1);
    nv = zeros(prm.ndim,1);
    for i = 1 : prm.ndim
       v{i} = 1:d0:prm.dim(i); 
       if ~any(v{i} == prm.dim(i))
           v{i} = [v{i}, prm.dim(i)];
       end
       nv(i) = numel(v{i}) - 1;
       v{i}(1) = 0;
    end
    nv(nv == 0) = 1;
    % Compute bwdist patchwise
    radim = cell(nv(1), nv(2), nv(3));
    msg = ['Computing bwdist for ' int2str(prod(nv)) ' patches'];
    disp(msg);
    tic
    for i = 1 : nv(1)
        for j = 1 : nv(2)
            for k = 1 : nv(3)
                start = [v{1}(i)+1,v{2}(j)+1,v{3}(k)+1];
                stop = [v{1}(i+1),v{2}(j+1),v{3}(k+1)];
                maskh = imcomplement(bwtree(start(1):stop(1), start(2):stop(2),start(3):stop(3)));
                radim{i,j,k} = prm.coordfnc(round(bwdist(maskh)));
            end
        end
    end
    radim = cell2mat(radim);
    clear mask;
end
clear bwtree maskh;
toc         

% Extract the radius
segment.avrad = NaN(segment.n,1);
for i = 1 : segment.n
    ind = segment.ind{i};
    val = mean(radim(ind));
    if isempty(val)
        val = meanh;
    end
    segment.avrad(i) = val*meanh;
end
clear radim v nv maskh start stop;    

% % Override the radius by user-provided values?
% if ~isempty(avrad)
%     segment.avrad = avrad;
% end
% if ~isempty(vessellength)
%     segment.L = vessellength;
% end
    
end

%-------------------------------------------
