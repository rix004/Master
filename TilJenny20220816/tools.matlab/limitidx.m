function [idx] = limitidx(idx, dim)

for i = 1 : numel(idx)
   idx{i}(idx{i} < 1) = 1;
   idx{i}(idx{i} > dim(i)) = dim(i);
end