function [fnc] = getintfnc(maxval)

nval = zeros(4,1);
fnc = cell(4,1);

nval(1) = intmax('uint8');
fnc{1} = @uint8;

nval(2) = intmax('uint16');
fnc{2} = @uint16;

nval(3) = intmax('uint32');
fnc{3} = @uint32;

nval(4) = intmax('uint64');
fnc{4} = @uint64;

ind = find(nval > maxval, 1);
fnc = fnc{ind};
