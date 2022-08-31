close all
clear;
iterations = 6;
iterations1 = 10;
L2_error = zeros(iterations1,iterations);
h = zeros(iterations1,iterations);

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
RandomTree.RadiusRate = 0.8;

% Domain
D = [0 1 0 1];
D_area = (D(2)-D(1))*(D(4)-D(3));

% Test problem
K_D = 2;
f = @(x,y) (2.*y -2.*y^2)*K_D +(2.*x -2.*x.^2)*K_D; 
p_exact = @(x,y) x.*(x-1).*y.*(y-1);

for iter1 = 1:iterations1
    Ncells = 50; 
    for iter2 = 1:iterations
        flag.case = 'Combinated';
        Tree = ChooseTree(flag.case,RandomTree,DT,D,Ncells);
        %DrawTree(Tree,150,'b',D);
        nodes = Tree.nodes; edges = Tree.edges;

        % Fix edge radii
        rel = edges(:,1)./edges(:,4);
        fix = find(rel<20);
        edges(fix,4)=edges(fix,1)/100;
        fix = find(edges(:,4)<1E-6);
        edges(fix,4)=1E-6;

        % Find terminal nodes
        [TNinfo,TNlogic]=FindTerminals(nodes,edges);
        
        %%% Voronoi diagram %%%
        [cells, vertices] = VoronoiDiagram(TNinfo,[D(1) D(1) D(2) D(2) D(1);D(3) D(4) D(4) D(3) D(3)]');
       
        
        [Grad_D,LHS,D_bvs,RHS,cell_center,cell_edges,cell_area,boundary_cells] = TPFA(cells,vertices,f,K_D,1,0,0);
        p_num = LHS\RHS;
        
        error = p_num-p_exact(cell_center(:,1),cell_center(:,2));
        l2_error = 0;
        for i =1:length(error)
            l2_error = l2_error + error(i)^2*cell_area(i);
        end

        L2_error(iter1,iter2) = sqrt(sum(l2_error));
        h(iter1,iter2) = sqrt(1/size(cells,1));
        Ncells=Ncells*2;
    end
    disp(iter1)
end

save('L2errorData','L2_error','h');
Plots('TPFAConvergenceVoronoi')