function [nodes] = newbranch3D(x_val,y_val,z_val,lev,edges,r_vect,x1,y1,z1,theta, distance,levels,root_radius,rate)
    if levels ~= 1
        x2 = x1+cosd(theta)*distance;
        y2 = y1+sind(theta)*distance;
        z2 = z1+distance;
        l2 = levels-1;
        dist = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
        x_val = [x_val x2];
        y_val = [y_val y2];
        lev = [lev l2];
        edges = [edges dist];
        r_vect = [r_vect root_radius];
        nodes = [x_val' y_val' lev' edges' r_vect'];       % lev er en kolonne som holder styr på hvilket nivå noden er på.
        line([x1 x2],[y1,y2],'LineWidth',root_radius);
        
        leftside = newbranch3D(x_val,y_val,z_val,lev,edges,r_vect,x2,y2,z2,theta+20,distance*0.5,levels-1,root_radius*rate,rate);
        
        x = leftside(:,1)';
        y = leftside(:,2)';
        z = leftside(:,3)';
        e = leftside(:,4)';
        r = leftside(:,5)';
        
        rightside = newbranch3D(x,y,l,z,e,r,x2,y2,z2,theta-20,distance*0.5,levels-1,root_radius*rate,rate);
        
        x_val = [rightside(:,1)'];
        y_val = [rightside(:,2)'];
        lev = [rightside(:,3)'];
        edges = [rightside(:,4)'];
        r_vect = [rightside(:,5)'];
        nodes = [x_val' y_val' lev' edges' r_vect'];
        pause(0.01);
        xlabel('x');
        ylabel('y');
    end
    nodes = [x_val' y_val' lev' edges' r_vect'];
     
end