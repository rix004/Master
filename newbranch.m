function [nodes] = newbranch(x_val,y_val,lev,x1,y1,theta, distance,levels,root_radius)
    if levels ~= 1
        x2 = x1+cosd(theta)*distance;
        y2 = y1+sind(theta)*distance;
        l2 = levels-1;
        x_val = [x_val x2];
        y_val = [y_val y2];
        lev = [lev l2];
        nodes = [x_val' y_val' lev'];       % lev er en kolonne som holder styr på hvilket nivå noden er på.
        line([x1 x2],[y1,y2],'LineWidth',root_radius);
        leftside = newbranch(x_val,y_val,lev,x2,y2,theta+20,distance*0.5,levels-1,root_radius*0.7);
        x = leftside(:,1)';
        y = leftside(:,2)';
        z = leftside(:,3)';
        rightside = newbranch(x,y,z,x2,y2,theta-20,distance*0.5,levels-1,root_radius*0.7);
        x_val = [rightside(:,1)'];
        y_val = [rightside(:,2)'];
        lev = [rightside(:,3)'];
        nodes = [x_val' y_val' lev'];
        pause(0.01);
        xlabel('x');
        ylabel('y');
    end
    nodes = [x_val' y_val' lev'];
     
end