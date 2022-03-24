% TEST TPFA
close all
%clc;
clear;

k = @(x,y) 1;
p_exact = @(x,y) x.*(x-1).*y.*(y-1);
nx = 64;
ny = nx;
dx=1/(nx-1);
dy=1/(ny-1);
h = sqrt(dx^2+dy^2);

[vertices, cells]=GridGeneration(nx,ny);
figure
for i = 1:size(cells,1)
        coords=vertices(cells{i},:);
        pgon = polyshape(coords(:,1),coords(:,2));
        pg = plot(pgon);
        hold on
end

f = @(x,y) 2.*y -2.*y^2 +2.*x -2.*x.^2; 
[LHS,RHS,cell_center,poly_edges] = TPFA(cells,vertices,f,k,1,0);
p = LHS\RHS;

X1 = reshape(cell_center(:,1),nx-1,ny-1);
X2 = reshape(cell_center(:,2),nx-1,ny-1);
P = reshape(p,nx-1,ny-1);
P_exact = reshape(p_exact(cell_center(:,1),cell_center(:,2)),nx-1,ny-1);

figure
subplot(1,2,1)
surf(X1,X2,P_exact);
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
title('Analytic solution','FontSize',18)
axis([0 1 0 1]);
subplot(1,2,2)
surf(X1,X2,P-P_exact);
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
title('Error','FontSize',18)
set(gcf,'Position',[100 100 1200 500]);
axis([0 1 0 1]);

error = p-p_exact(cell_center(:,1),cell_center(:,2));
l2_error = sqrt(sum(error.^2))*sqrt(dx*dy)

