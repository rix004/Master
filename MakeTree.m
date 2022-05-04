function [nodes] = MakeTree(x_val,y_val,lev,edges,r_vect,x1,y1,theta,distance,levels,root_radius,r_rate,d_rate,delta_theta)
    if levels ~= 1
        x2 = x1+cosd(theta)*distance;
        y2 = y1+sind(theta)*distance;
        dist = sqrt((x2-x1)^2+(y2-y1)^2);
        x_val = [x_val x2];
        y_val = [y_val y2];
        lev = [lev levels-1];
        edges = [edges dist];
        r_vect = [r_vect root_radius];
        nodes = [x_val' y_val' lev' edges' r_vect'];
        
        leftside = MakeTree(x_val,y_val,lev,edges,r_vect,x2,y2,theta-(delta_theta+1*(-1)^(round(rand()))*3.75*rand()),distance*d_rate,levels-1,root_radius*r_rate+0.00*rand(),r_rate,d_rate,delta_theta);
        
        x = leftside(:,1)';
        y = leftside(:,2)';
        z = leftside(:,3)';
        e = leftside(:,4)';
        r = leftside(:,5)';
        
        rightside = MakeTree(x,y,z,e,r,x2,y2,theta+(delta_theta+1*(-1)^(round(rand()))*3.75*rand()),distance*d_rate,levels-1,root_radius*r_rate+0.00*rand(),r_rate,d_rate,delta_theta);
        
        x_val = [rightside(:,1)'];
        y_val = [rightside(:,2)'];
        lev = [rightside(:,3)'];
        edges = [rightside(:,4)'];
        r_vect = [rightside(:,5)'];
        nodes = [x_val' y_val' lev' edges' r_vect'];
    end
    nodes = [x_val' y_val' lev' edges' r_vect'];
end