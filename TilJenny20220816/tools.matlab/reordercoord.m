function [segment] = reordercoord(segment, modifieds)

    disp('Reorder segment coordinates');
    tic
    c = segment.c;
    ind = segment.ind;    
    for i = 1 : segment.n
           
        ch = c{i};
        indh = ind{i};
        ncoord = size(ch,1);

        if ncoord == 1 || ncoord == 2
            continue;
        end

        % Only compute for the segments which need reordering
        if modifieds(i)

            % Distance
            D = pdist(single(ch));

            % Not neighbours
            D(D > sqrt(3)) = Inf;
            D = squareform(D);
            adjmat = D < Inf;    
            endpts = find(sum(adjmat, 2) == 2);
            nendpts = numel(endpts);

            % Make graph
            G = graph(D);
            
            if nendpts < 2
                endpts(1) = 1;
                endpts(2) = ncoord;
            elseif nendpts > 2
                % Find the two points of maximal distance
                D2 = distances(G);
                D2 = D2(endpts, endpts);
                [nx, ny] = size(D2);
                [ix, iy] = ndgrid(1:nx, 1:ny);
                D2 = D2(:);
                ix = ix(:);
                iy = iy(:);
                [~, maxind] = max(D2);
                idxmax = [ix(maxind), iy(maxind)];
                endpts = endpts(idxmax);
            end

            % Find the shortest path
            sp = shortestpath(G, endpts(1), endpts(2));
            segment.ncoord(i) = numel(sp);

            % Assign
            c{i} = ch(sp,:);
            ind{i} = indh(sp);    

        end
        
%         % Start and stop coordinates
%         start(i,:) = single(c{i}(1,:));
%         stop(i,:) = single(c{i}(end,:));
        
    end        
    segment.c = c;clear c;
    segment.ind = ind;clear ind;
    toc        
        
end
