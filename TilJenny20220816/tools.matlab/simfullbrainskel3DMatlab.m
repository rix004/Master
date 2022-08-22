function [bwskeleton, bwnode, bwterminal, bwtree] = simfullbrainskel3DMatlab(bwtree, prm)

if prm.tree.loadskeleton
    msg = ['Loading skeleton as ' prm.tree.pathskeleton];
    disp(msg);
    D = load(prm.tree.pathskeleton);
    bwskeleton = false(prm.dim);
    bwskeleton(D.ind) = true;
    clear D;
else
    % Thin image
%     bwtree = bwthin(bwtree, prm.tree.thiniter);
    
    disp('Generate tree skeleton');                
    msg = ['Skeletonize with minimum branch length: ' int2str(prm.tree.minbranchlength)];
    disp(msg)
    tic
    bwskeleton = bwskel(bwtree,'MinBranchLength',prm.tree.minbranchlength);
    toc
end

% Ensure skeleton is in tree - there is a bug in the scipy algorithm
% bwtree(bwskeleton) = true;

% Remove single voxel skeleton points
disp('Removing single voxel objects in skeleton');
bwskeleton = bwareaopen(bwskeleton, 2, prm.conn);

if prm.tree.saveskeleton
    msg = ['Saving skeleton as ' prm.tree.pathskeleton];
    disp(msg);
    ind = find(bwskeleton);    
    dim = prm.dim;
    save(prm.tree.pathskeleton, 'ind', 'dim', '-v7.3');
    clear ind;    
end

% Find interior node points
disp('Find branching points')
tic
bwnode = bwmorph3(bwskeleton, 'branchpoints');
toc

% Find leaves
disp('Find leaves');
tic
bwterminal = bwmorph3(bwskeleton, 'endpoints');
toc
