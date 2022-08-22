function [bwroot, aifpos] = getaifpos(bwbrain, bwterminal, aifpos, aifoption)

dim = size(bwterminal);
ndim = numel(dim);
if ndim == 2
    conn = 8;
elseif ndim == 3
    conn = 26;
end
% Default
bwroot = false(dim);
indroot = [];

if isequal(aifoption, 'mask')
    msg = ['Checking AIF candidates within maskfile ' aifpath];
    disp(msg);
    D = load(aifpath);
    bwaif = D.bw;
    clear D;
    bwroot = bwterminal & bwaif;    
    
elseif isequal(aifoption,'userdefined')    

    % Find the closest candidates
    faser = bwconncomp(bwterminal, conn);
    format = getintfnc(prod(dim));
    faser.PixelIdxList = (cellfun(format,faser.PixelIdxList,'UniformOutput', false))';       
    ind = faser.PixelIdxList;
    % Take the first indice, this is an ok approximation
    ind = cellfun(@(x)x(1),ind);    
    L = faser.NumObjects;
    meanc = zeros(L,3);
    [meanc(:,1), meanc(:,2), meanc(:,3)] = ind2sub(dim, ind); 
    clear faser;
    naif = size(aifpos,1);    
    bwroot = false(dim);    
    for i = 1 : naif
       dist = bsxfun(@minus, meanc, aifpos(i,:));
       dist = sqrt(sum(dist.^2,2));
    
       % Ensure not to pick the same terminal more times
       [~, minind] = min(dist);        
       indh = ind(minind);
       bwroot(indh) = true;    
       aifpos(i,:) = meanc(minind,:);
       
       % Remove this point
       ind(minind) = [];
       meanc(minind,:) = [];        
    end
    
elseif isequal(aifoption, 'outsidebrain')
    
    % Outside brain
    bwroot = bwterminal & ~bwbrain;

elseif isequal(aifoption, 'minx')
    
    % The terminal with the smallest x coordinate
    [faser,L] = bwlabeln(bwterminal, conn);
    cm = zeros(L,3);
    for i = 1 : L
        ind = find(faser == i);    
        clear c;    
        [c(:,1), c(:,2), c(:,3)] = ind2sub(dim, ind);
        cm(i,:) = mean(c,1);        
    end
    [~, ind] = min(cm(:,1));
    bwroot = faser == ind;
   
elseif isequal(aifoption, 'maxx')
    
    % The terminal with the largest x coordinate
    [faser,L] = bwlabeln(bwterminal, conn);
    cm = zeros(L,3);
    for i = 1 : L
        ind = find(faser == i);    
        clear c;    
        [c(:,1), c(:,2), c(:,3)] = ind2sub(dim, ind);
        cm(i,:) = mean(c,1);        
    end
    [~, ind] = max(cm(:,1));
    bwroot = faser == ind;   
    
% elseif isequal(option, 'userdefined')
%     bwroot0 = bwroot;
%     bwroot = false(dim);
%     [faser,L] = bwlabeln(bwroot0);
%     for i = 1 : L
%         reghere = faser == i;
%         
%     end
%     bwroot = bwterminal;
    
end