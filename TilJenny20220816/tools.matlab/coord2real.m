function [creal] = coord2real(c,h)

c = double(c);
n = size(c,1);
h = repmat(h,n,1);
creal = c.*h - h/2;