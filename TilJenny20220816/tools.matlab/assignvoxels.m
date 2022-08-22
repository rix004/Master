function [segment] = assignvoxels(bw, node, segment, prm)

fnc = getintfnc(segment.n);
faser = fnc(zeros(prm.dim, 'logical'));
for i = 1 : segment.n
   faser(segment.ind{i}) = i; 
end
se = ones(3,3,3);
vals = -1:1;
cmain = 0;
while 1
    cmain = cmain + 1;
    msg = ['Dilation number ' int2str(cmain)];
    disp(msg);
    bwf = faser > 0;
    diff = imdilate(bwf, se);
    diff(bwf) = false;clear bwf;    
    diff = diff & bw;
    ind0 = find(diff);clear diff;
    n = numel(ind0);
    msg = ['Number of voxels to assign: ' int2str(n)];
    disp(msg);
    if n == 0
        break;
    end
    clear c;
%     if cmain > 100
%         a = 2;
%     end
    [c(:,1), c(:,2), c(:,3)] = ind2sub(prm.dim, ind0);
    nind = numel(ind0);
    neigh = fnc(zeros(nind, 1));
    count = 0;
    for i = vals
        for j = vals
            for k = vals
                count = count + 1;
                c1 = c(:,1) + i;
                c2 = c(:,2) + j;
                c3 = c(:,3) + k;
                c1(c1 < 1) = 1;c1(c1 > prm.dim(1)) = prm.dim(1);
                c2(c2 < 1) = 1;c2(c2 > prm.dim(2)) = prm.dim(2);
                c3(c3 < 1) = 1;c3(c3 > prm.dim(3)) = prm.dim(3);
                ind = sub2ind(prm.dim, c1, c2, c3);
                % Only assign the voxels that have not already been
                % assigned
                bwh = neigh == 0;
%                 neigh = max(neigh(bw), faser(ind(bw)));
                neigh(bwh) = faser(ind(bwh));
            end
        end        
    end
    faser(ind0) = neigh;        
end
clear bw;
faser(cell2mat(node.ind)) = 0;
segment.distrind = label2idx(faser)';
fnc = getintfnc(prod(prm.dim));
for i = 1 : segment.n
   segment.distrind{i} = fnc(segment.distrind{i}); 
end

end