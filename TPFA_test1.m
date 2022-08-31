% TEST TPFA
close all
%clc;
clear;

K_D = 2;
p_exact = @(x,y) x.*(x-1).*y.*(y-1);
nx = 12;
ny = nx;
dx=1/(nx-1);
dy=1/(ny-1);
h = sqrt(dx^2+dy^2)

[vertices, cells]=GridGeneration(nx,ny,[0 1 0 1]);
figure
for i = 1:size(cells,1)
        coords=vertices(cells{i},:);
        pgon = polyshape(coords(:,1),coords(:,2));
        pg = plot(pgon);
        hold on
end

f = @(x,y) (2.*y -2.*y^2)*K_D +(2.*x -2.*x.^2)*K_D; 
[Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells,bv_out] = TPFA(cells,vertices,f,K_D,1,0);
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
zlabel('u','FontSize',15)
title('Analytic solution','FontSize',18)
axis([0 1 0 1]);
subplot(1,2,2)
surf(X1,X2,P-P_exact);
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
zlabel('error','FontSize',15)
title('u_{analytic} - u_h','FontSize',18)
set(gcf,'Position',[100 100 1200 500]);
axis([0 1 0 1]);
error = p-p_exact(cell_center(:,1),cell_center(:,2));

l2_error = sqrt(sum(error.^2))*sqrt(dx*dy)
PEX = p_exact(cell_center(:,1),cell_center(:,2));
max(PEX)
% h_vect = [0.4714 0.2020 0.0943 0.0456 0.0224]';
% l2_vect = [0.0052 9.8575e-04 2.1690e-04 5.0907e-05 1.2333e-05]';

% figure
% loglog(h_vect,l2_vect,'*-','LineWidth',2)
% hold on
% loglog([0.25/2 0.25/4 0.25/8 0.25/16 0.25/32],[0.05 0.025/2 0.025/8 0.025/32 0.025/(32*4)],'-','LineWidth',2)
% hold on
% xlabel('h','FontSize',15)
% ylabel('Error','FontSize',15)
% lgd = legend('Convergence on square grid','2nd order reference');
% lgd.FontSize=(14);
% legend('Location','northwest')