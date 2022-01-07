function [nodes] = newbranch3D(x_val,y_val,z_val,x1,y1,z1,theta, distance,levels,root_radius)
    if levels ~= 1
        x2 = x1+cosd(theta)*distance;
        y2 = y1+sind(theta)*distance;
        z2 = z1+distance;
        x_val = [x_val x2];
        y_val = [y_val y2];
        z_val = [z_val z2];
        nodes = [x_val' y_val' z_val'];
        line([x1 x2],[y1,y2],[z1,z2],'LineWidth',root_radius);
        leftside = newbranch3D(x_val,y_val,z_val,x2,y2,z2,theta+20,distance*0.5,levels-1,root_radius*0.7);
        x = leftside(:,1)';
        y = leftside(:,2)';
        z = leftside(:,3)';
        rightside = newbranch3D(x,y,z,x2,y2,z2,theta-20,distance*0.5,levels-1,root_radius*0.7);
        x_val = [x_val rightside(:,1)'];
        y_val = [y_val rightside(:,2)'];
        z_val = [z_val rightside(:,3)'];
        nodes = [x_val' y_val' z_val'];
        %pause(0.01);
        view(3);
        xlabel('x');
        ylabel('y');
        zlabel('z');
    end
    nodes = [x_val' y_val' z_val'];
     
end