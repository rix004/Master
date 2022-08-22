close all
clear;
Ncells = 100; 
iterations = 5;
iterations1 = 1;
L2_error = zeros(iterations,iterations1);
h = zeros(iterations,iterations1);

% Deterministic tree data
RootNode = [0.5 0];
DT.Levels = 4;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.05;    % mm
DT.RadiusRate = 0.7;
DT.TrunkLength = 1/sqrt(2)/2; %mm
DT.LengthRate = 1/sqrt(2);

% Random tree data
RandomTree.TrunkRadius = DT.TrunkRadius*DT.RadiusRate^DT.Levels;
RandomTree.RadiusRate = 0.7;

% Domain
D = [0 1 0 1];
D_area = (D(2)-D(1))*(D(4)-D(3));

for iter1 = 1:iterations1
    for iter2 = 1:iterations
        trees = {'Deterministic','Random','Combinated','Half Deterministic'};
        flag.case = trees{3};
        Tree = ChooseTree(flag.case,RandomTree,DT,D,Ncells);
        DrawTree(Tree,150,'b',D);
        nodes = Tree.nodes; edges = Tree.edges;
        
        % Find terminal nodes
        [TNinfo,TNlogic]=FindTerminals(nodes,edges);
        
        %%% Voronoi diagram %%%
        [cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');
        
        f = @(x,y) 2.*y -2.*y^2 +2.*x -2.*x.^2;
        K_D = 1;
        p_exact = @(x,y) x.*(x-1).*y.*(y-1);
        
        [Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells] = TPFA(cells,vertices,f,K_D,1,0,0);
        p_num = LHS\RHS;
        
        error = p_num-p_exact(cell_center(:,1),cell_center(:,2));
        l2_error = 0;
        for i =1:length(error)
            l2_error = l2_error + error(i)^2*cell_area(i)/D_area;
        end

        L2_error(iter1,iter2) = sqrt(sum(l2_error));
        h(iter1,iter2) = sqrt(1/size(cells,1));
        Ncells=Ncells*2;
    end
end

figure
loglog([0.5 0.25 0.25/2 0.25/4 0.25/8],[0.02 0.01 0.005 0.0025 0.0025/2],'-','LineWidth',3,'Color','r')
hold on

for iter3 = 1:size(L2_error,1)
    loglog(h(iter3,:),L2_error(iter3,:),'-.','LineWidth',1.5,'Color',[0    0.5686    0.7157])
    hold on
end

mean_vect = zeros(1,size(L2_error,2));
mean_h = zeros(1,size(L2_error,2));
for iter4 = 1:size(L2_error,2)
    mean_vect(iter4)=mean(L2_error(:,iter4));
    mean_h(iter4)=mean(h(:,iter4));
end

std_dev = std(L2_error);
errorbar(mean_h,mean_vect,std_dev,'x-','LineWidth',2.5,'Color',[0    0.1882    0.9059])
hold on

xlabel('h','FontSize',16)
ylabel('Error','FontSize',16)
lgd = legend('1st order reference','','','','','','','','','','Error','Average error');
lgd.FontSize=(14);
legend('Location','northwest')