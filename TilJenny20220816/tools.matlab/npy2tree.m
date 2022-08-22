% Sett stien, kjøre filen
% setpath
clear;
% Deretter kjør filen
% npy2tree;

% Les inn data
bw = readNPY('DLAtre.npy') > 0;
% bw = bw(230:270, 230:270);
% bw = imfill(bw, 'holes');

% ------------- Parametre ---------------

% Pixel størrelse
prm.h = [1,1,1]*2;

% Posisjon for roten/røttene i fysiske størrelser relatert til pixel
% størrelsen
% prm.tree.aifpos = [26, 20, 1];
prm.tree.aifpos = [262, 251, 1]*2;

% Ikke still på disse
prm.dim = [size(bw), 1];
prm.tree.minbranchlength = 1;
prm.tree.loadskeleton = false;
prm.tree.saveskeleton = false;
prm.tree.aifoption = 'userdefined';
prm.tree.savememory = false;
prm.paths.belongroot = {'arterial'};
prm.tree.assignvoxels = true;

% Lag tree
tree = simfullbraintreefun(bw, [], bw > -Inf, prm.tree, prm);

% Fjerne doble kanter. Må først sortere slik at unique finner de doble
% kantene
tree.segment.nodeconn = sort(tree.segment.nodeconn, 2);
[~, idx] = unique(tree.segment.nodeconn, 'rows');
fn = fieldnames(tree.segment);
for i = 1 : numel(fn)
   try
       tree.segment.(fn{i}) = tree.segment.(fn{i})(idx,:);
   catch
      msg = ['Coult not restrict ' fn{i}];
      disp(msg);
   end
end

% % Få over på eget format
DLAedges(:,1) = tree.segment.L;
DLAedges(:,2:3)=tree.segment.nodeconn;
G = graph(DLAedges(:,2),DLAedges(:,3));
D = distances(G);
levels = D(:,tree.node.isroot)+1;
DLAnodes(:,1:3)=[tree.node.avx(:,1:2) levels];


% Sett radius
startradius = 1;
radiusrate = 0.95;
for i = 1:size(DLAedges,1)
    nodelevel = DLAnodes(DLAedges(i,3),3);
    radius = startradius*radiusrate^(nodelevel-1);
    DLAedges(i,4)=radius;
end
tree.segment.avrad = DLAedges(:,4);
    
% Plot tre
p = plotGraph(tree,'2D', 'nodelabel', false, 'edgelabel', false,'linewidth',DLAedges(:,4)*5);
highlight(p, tree.node.isroot, ...
    'Nodecolor', 'r', ...
    'MarkerSize', 8)

% Lagre treet
[~, ~, ~] = mkdir('data');
pathsave = fullfile('data', 'DLA-tree.mat');
save(pathsave, 'tree', 'prm');