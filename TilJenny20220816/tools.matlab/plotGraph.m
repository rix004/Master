function [H, G, local2G, G2local] = plotGraph(varargin)

tree = varargin{1};
type = varargin{2};
fn = 'avx';
nodelabel = true;
edgelabel = true;
for i = 3 : 2 : nargin
   switch lower(varargin{i})
       case 'h'
           H = varargin{i+1};
       case 'type'
           type = varargin{i+1};
       case 'lims'
           lims = varargin{i+1};
       case 'cols'
           cols = varargin{i+1};   
       case 'linewidth'
          lw = varargin{i+1};
       case 'coordvar'
           fn = varargin{i+1};
       case 'nodelabel'
           nodelabel = varargin{i+1};
       case 'edgelabel'
           edgelabel = varargin{i+1};           
       otherwise
      error('Wrong option');     
   end
end
if ~exist('lw', 'var')
    lw = ones(tree.segment.n,1);
end
% lims = [min(data.master.tree.node.avx, [], 1);max(data.master.tree.node.avx, [], 1)]

[G, local2G, G2local] = makegraph(tree.segment.nodeconn, tree.segment.n, tree.node.n);

% Reorder correctly
lw = lw(local2G);

if ~exist('H', 'var')
    H = figure;
else    
    hold on;
end

if ~exist('lims', 'var')
    lims = [min(tree.node.(fn), [], 1);max(tree.node.(fn), [], 1)];
end
flaglims = true;
if any(lims(2,:) - lims(1,:) == 0)
    flaglims = false;
end

if ~exist('cols', 'var')
    cols = {'r','k'};
end

if isequal(type, 'force')  
%     H = plot(G, 'Linewidth', lw);
    H = plot(G,'layout', 'force', 'Linewidth', lw);
elseif isequal(type, '2D')
    H = plot(G, 'Linewidth', lw, ...
               'XData',tree.node.(fn)(:,1), ...
               'YData',tree.node.(fn)(:,2));
        if flaglims
            xlim(lims(:,1));
            ylim(lims(:,2));       
        end
elseif isequal(type, '3D')
    H = plot(G, 'Linewidth', lw, ...
               'XData',tree.node.(fn)(:,1), ...
               'YData',tree.node.(fn)(:,2), ...
               'ZData',tree.node.(fn)(:,3), ...
               'EdgeAlpha',0.3);
       xlim(lims(:,1));
       ylim(lims(:,2));
       zlim(lims(:,3));
       view(-60, 55)
end
if islogical(nodelabel)
    if nodelabel
        set(H, 'NodeLabel', G.Nodes.Name);
    else
        set(H, 'NodeLabel', []);
    end
else
    set(H, 'NodeLabel', nodelabel);
end

if islogical(edgelabel)
    if edgelabel
        labeledge(H,G.Edges.EndNodes(:,1), G.Edges.EndNodes(:,2), G.Edges.index) 
   end
else
   labeledge(H,G.Edges.EndNodes(:,1), G.Edges.EndNodes(:,2), edgelabel(G.Edges.index))  
end
grid on
hold on;
% axis equal;
axis image
if all(tree.node.avc(:,3) == 1)
    view([90 90])
end
axis tight;

% if hl
% %     highlight(H,tree.node.indterminal, 'MarkerSize', 4);
% %     hold on;    
% %     highlight(H,tree.node.indroot, 'MarkerSize', 6);
%     highlight(H,tree.node.indroot, 'NodeColor', cols{1}, 'MarkerSize', 4);
%     hold on;
%     highlight(H,tree.node.indterminal, 'NodeColor', cols{2}, 'MarkerSize', 4);
%     hold on;    
%     
% end

