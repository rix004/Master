clc;clear;
axis(gca,'equal')
distance = 0.7;
levels = 7;
num_nodes = 2^(levels-1);
x1 = 0;
y1 = 0;
z1 = 0;
degrees = 90;
root_radius = 5;
x_val = [x1];
y_val = [y1];
z_val = [z1];


allnodes = newbranch1(x_val,y_val,[levels],x1,y1,degrees,distance,levels,root_radius)
%allnodes3D = newbranch3D(x_val,y_val,z_val,x1,y1,z1,degrees,distance,levels,root_radius);
%allnodes3D = allnodes3D(length(allnodes3D)+1-num_nodes:end,:)

