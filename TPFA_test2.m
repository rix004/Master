close all
clc; clear;
RootNode = [0.5 0];

% Deterministic tree data
DT.Levels = 4;
DT.StartPos = RootNode;
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.01;    
DT.RadiusRate = 0.5;
DT.TrunkLength = 0.4;
DT.LengthRate = 0.65;

% Random tree data
RandomTree.TrunkRadius = 0.0001;
RandomTree.RadiusRate = 0.5;
D = [0 1 0 1];

trees = {'Deterministic','Random','Combinated','Half Deterministic'};
flag.case = trees{3};

numb_tests = 10;
numb_doubles = 6;
L2_error = zeros(numb_tests,numb_doubles);
for iter1 = 1:numb_tests
    for iter2 = 1:numb_doubles
        Np = 10*2^(iter2);
        Tree = ChooseTree(flag.case,RandomTree,DT,D,Np);
        nodes = Tree.nodes; edges = Tree.edges;
        
        % Find terminal nodes
        Ne = length(edges);
        Nn = length(nodes);
        StartNodes = edges(:,2);
        EndNodes = edges(:,3);
        find_tn = ismember(EndNodes,StartNodes);
        Tn = [];
        Te = [];
        Bc_nodes = [1;zeros(Nn-1,1)];
        for i = 1:length(EndNodes)
            if find_tn(i) == 0
                Tn = [Tn;nodes(EndNodes(i),1:2),i];
                Te = [Te;i];
                Bc_nodes(EndNodes(i))=1;
            end
        end
        
        %%% Set boundary values %%%
        Bc_vals = [1;zeros(Nn-1,1)];
        for i = 2:length(Bc_vals)
            Bc_vals(i)=Bc_nodes(i)*0;
        end
        
        %%% VORONOI DIAGRAM %%%
        [cells, vertices] = VoronoiDiagram(Tn,[0 0 1 1 0;0 1 1 0 0]');
        
        
        f = @(x,y) 2.*y -2.*y^2 +2.*x -2.*x.^2;
        k = @(x,y) 1;
        p_exact = @(x,y) x.*(x-1).*y.*(y-1);
        
        [Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells] = TPFA(cells,vertices,f,k,1,0,0);
        p_num = LHS\RHS;
        
        error = p_num-p_exact(cell_center(:,1),cell_center(:,2));
        l2_error = 0;
        for i =1:length(error)
            l2_error = l2_error + error(i)^2*cell_area(i);
        end

        L2_error(iter1,iter2) = sqrt(sum(l2_error));
        h(iter1,iter2) = sqrt(1/size(cells,1));
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