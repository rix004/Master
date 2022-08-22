function [segment] = compvessellength(segment, node, prm, update)

% segment.dist = cell(segment.n,1);
updateL = false;
if ~isfield(segment, 'L')
    updateL = true;
    segment.L = NaN(segment.n,1);     
end
updategeodes = false;
if ~isfield(segment, 'geodesdist')
    updategeodes = true;
    segment.geodesdist = cell(segment.n,1);        
end

tic
for i = 1 : segment.n
    
    if update(i) == 0
        continue;
    end
    % Node connection
    nodeconn = segment.nodeconn(i,:);

    % Coordinate to start from
    cnode1 = node.avc(nodeconn(1,1),:);
    cnode2 = node.avc(nodeconn(1,2),:);

    % The coordinates of this segment
    c = segment.c{i};

    % Combined with nodes
    ccomb = [cnode1;c;cnode2];

    % Real coordinates of cell centers
    x = coord2real(double(ccomb), prm.h);

    % Find cumulative distance                
    d = x(2:end,:,:) - x(1:end-1,:,:);
    d = sqrt(sum(d.^2,2));                              
    dist = [0;d];        
    cumdist = cumsum(dist);                   

    % Length of segment
    L = cumdist(end);

    % Geodesic distance
    if updategeodes
        segment.geodesdist{i} = single(cumdist);
    end

    % Lenght of segment
    if updateL
        if segment.ncoord(i) == 0
            % Length is length between node centers
            dist = node.avx(nodeconn(1),:) - node.avx(nodeconn(2),:);
            dist = sqrt(sum(dist.^2));
            segment.L(i) = dist;           
        else
            segment.L(i) = L;            
        end
    end

    
end

num = 1;
segment.num = cell(segment.n,1);
for i = 1 : segment.n
    % Numbering of segments
    numnew = num + segment.ncoord(i);
    segment.num{i} = (num:(numnew-1))';
    num = numnew;
end    
toc