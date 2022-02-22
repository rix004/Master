% Testing
clc; clear;        
levels = 8;
R_trunk = 0.1;   
L_trunk = 1;      
x0 = 0;
y0 = 0;
theta = 90;
R_rate = 0.5;
L_rate = 0.7;
rotation_angle = 90;
[nodes, edges]= GetTree(levels,R_trunk,L_trunk,x0,y0,theta,R_rate,L_rate,rotation_angle);

TestTree.nodes = nodes;
TestTree.edges = edges;
TestTree.radius_ratio = R_rate;

%TestTree.nodes = [1 0 0 1];
%TestTree.edges = [1 1 2 1;1 1 3 1];
[x,y,i,d]= NearestNeighbour(-0.6,0.6,nodes)
NewTree = AddOnePoint(TestTree,-0.6,0.7);
%NewTree.nodes
%NewTree.edges
line([NewTree.nodes(end-1,2) NewTree.nodes(end,2)],[NewTree.nodes(end-1,3),NewTree.nodes(end,3)],'LineWidth',NewTree.edges(end,4),'Color',[0.8500, 0.3250, 0.0980]);
