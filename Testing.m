% Testing
clc; clear;     
DT.Levels = 4;
DT.StartPos = [0 0];
DT.StartAngle = 90;
DT.RotationAngle = 90;
DT.TrunkRadius = 0.1;    % Radius of trunk (mm)
DT.RadiusRate = 0.5;
DT.TrunkLength = 1;      % Length of trunk (mm)
DT.LengthRate = 0.7;

DetTree= GetTree(DT);

nodes = DetTree.nodes;
edges = DetTree.edges;
%[x,y,i,d]=NearestNode(-1.7615,-0.0700,nodes)
%[edge_nr,dist]=NearestEdge(-1.7615,-0.0700,DetTree)
testtree.nodes = nodes;
testtree.edges = edges;
testtree.TrunkRadius = 2;
testtree.RadiusRate = 0.5;
rrt_tree = RRT_Tree(testtree,-2,2,-2,2,100);
DrawTree(rrt_tree)

% levels = 5;
% R_trunk = 0.1;   
% L_trunk = 1;      
% R_rate = 0.5;
% L_rate = 0.5;
% rotation_angle = 90;
% [nodes, edges]= GetTree(levels,R_trunk,L_trunk,0,0,90,R_rate,L_rate,rotation_angle);
% 
% tree.nodes = nodes;
% tree.edges = edges;
% x = -0.6;
% y = 1.4;
% [edgenr,distance,onEdge,Point] = NearestEdge(x,y,tree)
% plot(x,y,'*')
% hold on
% plot(Point(1),Point(2),'.')


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


%NewTree = AddOnePoint(TestTree,-0.6,0.7);
%NewTree.nodes
%NewTree.edges
%line([NewTree.nodes(end-1,2) NewTree.nodes(end,2)],[NewTree.nodes(end-1,3),NewTree.nodes(end,3)],'LineWidth',NewTree.edges(end,4),'Color',[0.8500, 0.3250, 0.0980]);

% DrawTree(TestTree);
