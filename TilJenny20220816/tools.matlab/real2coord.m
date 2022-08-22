function [c] = real2coord(x,h,dim,flag)

if isempty(x)
    c = [];
    return;
end
n = size(x,1);
h = repmat(h,n,1);
c = x./h+0.5;
if flag
    c = round(c);
end
c(c<1) = 1;
c = min(c, dim);
