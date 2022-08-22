% SHOWALL Show one plane awt a time for stack.
%
%   SHOWALL(IM) Shows all planes in a 3D stack succesively when hitting ENTER. 
%
%   SHOWALL(IM1,IM2,...) There can be several input images, then they are 
%   plottet in different figures.  They have to have the same number of 
%   planes. If the last argument is 'colorbar', then the images are plotted 
%   with the colorbar option. 
%

function [] = showall(varargin)

numim = 0;
col = 0;
c = 0;
clim = [];
while c < nargin
    c = c + 1;
    
    if isequal(varargin{c},'colorbar')
        col = 1;
        continue;
    end
    
    if isequal(varargin{c},'clim')
        clim = varargin{c+1};
        c = c + 1;
        continue;
    end
    
    numim = numim + 1;
    f{numim} = (varargin{c});
    
end
clear varargin;
        



F = f{1};
[M N O P] = size(F);

if P > 1
    t = 1;
    while 1        
        for j = 1 : numim
            F = full(f{j});
            if isempty(clim)
                show(F(:,:,:,t),j);            
            else
                show(F(:,:,:,t),j, clim);
            end
        end
        msg = ['Time point ' int2str(t)];
        disp(msg);
        
        t = input('Continue: 0, Time: timepoint, Next timepoint: Enter ');

        if isequal(t,0)
            return
        elseif isempty(t)
            t = t + 1;
        else
            t = t;
        end                        

        if t > O
            break;
        end

    end
else


    niter = 1;
    while 1

        for j = 1 : numim
            F = full(f{j});

            if isempty(clim)
                show(F(:,:,niter),j)        
            else
                show(F(:,:,niter), j, clim);
            end
            if col == 1
                colorbar;
            end
        end
        disp(sprintf('Plane %i',niter))

        c = input('Continue: 0, Plane: planenumber, Next plane: Enter ');

        if isequal(c,0)
            return
        elseif isempty(c)
            niter = niter + 1;
        else
            niter = c;
        end                        

        if niter > O
            break;
        end
    end
end
    