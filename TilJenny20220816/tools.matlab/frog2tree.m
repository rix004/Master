% Sett stien, kjøre filen
% setpath

% Deretter kjør filen
% npy2tree;

% Les inn data
D = load('frog.mat');
bw = false(D.dim);
bw(D.ind) = true;

% ------------- Parametre ---------------

% Pixel størrelse
prm.h = D.h;

clear D
% Posisjon for roten/røttene i fysiske størrelser relatert til pixel
% størrelsen
% prm.tree.aifpos = [26, 20, 1];
prm.tree.aifpos = [632, 209,1;
                    632, 308,1];
prm.tree.aifpos = prm.tree.aifpos.*prm.h;

% Ikke still på disse
prm.dim = [size(bw), 1];
prm.tree.minbranchlength = 1;
prm.tree.loadskeleton = false;
prm.tree.saveskeleton = false;
prm.tree.aifoption = 'userdefined';
prm.tree.savememory = false;
prm.paths.belongroot = {'arterial', 'arterial'};
prm.tree.assignvoxels = false;

% Lag tree
tree = simfullbraintreefun(bw, [], bw > -Inf, prm);

% Plot tre
p = plotGraph(tree,'2D', 'nodelabel', false, 'edgelabel', false);
highlight(p, tree.node.isroot, ...
    'Nodecolor', 'r', ...
    'MarkerSize', 8)

% Lagre treet
[~, ~, ~] = mkdir('data');
pathsave = fullfile('data', 'DLA-tree.mat');
save(pathsave, 'tree', 'prm');