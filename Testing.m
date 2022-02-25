% Testing
clc; clear;        
levels = 5;
R_trunk = 0.1;   
L_trunk = 1;      
R_rate = 0.5;
L_rate = 0.5;
rotation_angle = 90;
[nodes, edges]= GetTree(levels,R_trunk,L_trunk,0,0,90,R_rate,L_rate,rotation_angle);

tree.nodes = nodes;
tree.edges = edges;
x = -0.6;
y = 1.4;
[edgenr,distance,onEdge,Point] = NearestEdge(x,y,tree)
plot(x,y,'*')
hold on
plot(Point(1),Point(2),'.')

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
