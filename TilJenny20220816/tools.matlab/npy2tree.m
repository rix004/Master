function[DLAtree]=npy2tree(n_particles,version,RootRadius,TerminalRadius)
% Les inn data
pathload = ['Tree' int2str(n_particles) 'particles_v' int2str(version) '.npy'];
bw = readNPY(pathload)>0;
% bw = bw(230:270, 230:270);
% bw = imfill(bw, 'holes');

% ------------- Parametre ---------------

% Pixel størrelse
prm.h = [1,1,1]/500;

% Posisjon for roten/røttene i fysiske størrelser relatert til pixel
% størrelsen
% prm.tree.aifpos = [26, 20, 1];
prm.tree.aifpos = [262, 251, 1]/500;

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

% Translate to my format
DLAedges(:,1) = tree.segment.L;
DLAedges(:,2:3)=tree.segment.nodeconn;
G = graph(DLAedges(:,2),DLAedges(:,3));
D = distances(G);
levels = D(:,tree.node.isroot);
DLAnodes(:,1:3)=[tree.node.avx(:,1:2) levels];

newDLAnodes = sortrows(DLAnodes,3);

% Correction of directions of edges
e = DLAedges(:,2:3);
for i = 1:size(e,1)
    if DLAnodes(e(i,1),3) > DLAnodes(e(i,2),3)
        DLAedges(i,2)=e(i,2);
        DLAedges(i,3)=e(i,1);
    end
end

% Correction of order of edges
newDLAedges = zeros(size(DLAedges,1),3);
for i = 1:size(DLAedges,1)
    NewEdgeNr1=find(ismember(newDLAnodes(:,1:2),DLAnodes(DLAedges(i,2),1:2),"rows"));
    NewEdgeNr2=find(ismember(newDLAnodes(:,1:2),DLAnodes(DLAedges(i,3),1:2),"rows"));
    newDLAedges(i,:) = [DLAedges(i,1) NewEdgeNr1 NewEdgeNr2];
end

% Sett radius
newDLAedges(:,4)=0;
DLAtree.nodes = newDLAnodes;
DLAtree.edges = newDLAedges;
DLAtree.edges(:,4)=1;

DLAtree.RootNodeIdx = find(newDLAnodes(:,3)==0);
[~,TermLog]=FindTerminals(DLAtree);
TermIndexes = find(TermLog==1);

[Graph, ~, ~]=makegraph(DLAtree.edges(:,2:3),size(DLAtree.edges,1),size(DLAtree.nodes,1));
NodeLevels = zeros(size(DLAtree.nodes,1),1);
for i = 1:length(TermIndexes)
    SP = shortestpath(Graph,TermIndexes(i),1);
    for j = 1:size(DLAtree.nodes,1)
        if any(ismember(SP,j)) && NodeLevels(j) < find(ismember(SP,j))-1
            NodeLevels(j) = find(ismember(SP,j))-1;
        end
    end
end


power = log(RootRadius/TerminalRadius)/(max(DLAtree.nodes(:,3))-1);
RadiusRate = exp(power);
for i = 2:length(NodeLevels)
    edge = find(DLAtree.edges(:,3)==i);
    edgerad = TerminalRadius*RadiusRate^(NodeLevels(i));
    DLAtree.edges(edge,4)=edgerad;
end



end