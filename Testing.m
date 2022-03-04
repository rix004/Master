% Testing
clc; clear;     
DT.Levels = 8;
DT.StartPos = [0 0];
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.1;
DT.RadiusRate = 0.5;
DT.TrunkLength = 1;      
DT.LengthRate = 0.7;

DetTree= GetTree(DT);
DrawTree(DetTree)
nodes = DetTree.nodes;
edges = DetTree.edges;

RandomTree.nodes = nodes;
RandomTree.edges = edges;
RandomTree.TrunkRadius = 2;
RandomTree.RadiusRate = 0.5;
rrt_tree = RRT_Tree(RandomTree,-2,2,-2,2,1000);
DrawTree(rrt_tree)


%%% VORONOIDIAGRAM %%%
%bs_int=[.2 .8 .8 .2;.6 .6 .2 .2]';
% bs_ext=[-.8 1.80 .5 -.8;-.05 -.05 1.7 -.05]';
% [X,Y] = meshgrid(-2:.5:2, -2:.5:2);
% X=X(:);Y=Y(:);
% [V,C,XY]=VoronoiLimit(X,Y,'figure','on');
% 
% figure()
% for i = 1:length(C)
%     coords = V(C{i},:);
%     pgon = polyshape(coords(:,1),coords(:,2));
%     pg = plot(pgon);
%     hold on
% end





% tree.nodes = nodes;
% tree.edges = edges;
% x = -0.6;
% y = 1.4;
% [edgenr,distance,onEdge,Point] = NearestEdge(x,y,tree)
% plot(x,y,'*')
% hold on
% plot(Point(1),Point(2),'.')


%%% Eksperimentering med ny kode for  binært tre

% values.Levels = 9;
% values.Angle = 90;
% values.RotationAngle = 90;
% values.Radius = 0.1;    % Radius of trunk (mm)
% values.RadiusRate = 0.5;
% values.Length = 1;      % Length of trunk (mm)
% values.LengthRate = 0.7;
% tree.xval = 0; 
% tree.yval = 0;
% tree.lev = [values.Levels];
% tree.edges = [0];
% tree.radius = [0];

%allnodes = MakeTree1(tree,0,0,values)

% TestTree.LengthRate = L_rate;
% TestTree.RadiRate = R_rate;
% TestTree.RotAngle = rotation_angle;
% TestTree.TrunkRadi = R_trunk;
% TestTree.TrunkLength = L_trunk;
% TestTree.nodes = [1 0 0 1]
% TestTree.edges=[];
% for i = 1:levels
%     NewTree = BinaryTree(TestTree);
%     TestTree.nodes = NewTree.nodes;
%     TestTree.edges = NewTree.edges;
% end
