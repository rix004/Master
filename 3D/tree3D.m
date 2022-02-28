function [nodes] = tree3D(levels,x_val,y_val,z_val,lev,edges,r_vect,x1,y1,z1,theta,phi,distance,root_radius,r_rate,d_rate,delta_theta,delta_phi)
    if levels ~= 1
        x2 = x1+sind(phi)*cosd(theta)*distance;
        y2 = y1+sind(phi)*sind(theta)*distance;
        z2 = z1+cosd(phi)*distance;
        dist = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
        x_val = [x_val x2];
        y_val = [y_val y2];
        z_val = [z_val z2];
        lev = [lev levels-1];
        edges = [edges dist];
        r_vect = [r_vect root_radius];
        %nodes = [x_val' y_val' z_val' lev' edges' r_vect'];       % lev er en kolonne som holder styr på hvilket nivå noden er på.
        line([x1 x2],[y1,y2],[z1,z2],'LineWidth',root_radius*20);
        leftside = tree3D(levels-1,x_val,y_val,z_val,lev,edges,r_vect,x2,y2,z2,theta-delta_theta*(0.5 + 1*rand()),phi+delta_phi,distance*d_rate,root_radius*r_rate,r_rate,d_rate,delta_theta,90*rand());
        
        x = leftside(:,1)';
        y = leftside(:,2)';
        z = leftside(:,3)';
        l = leftside(:,4)';
        e = leftside(:,5)';
        r = leftside(:,6)';
        
        
        rightside = tree3D(levels-1,x,y,z,l,e,r,x2,y2,z2,theta+delta_theta*(0.5 + 1*rand()),phi+delta_phi,distance*d_rate,root_radius*r_rate,r_rate,d_rate,delta_theta,90*rand());
        
        x_val = [rightside(:,1)'];
        y_val = [rightside(:,2)'];
        z_val = [rightside(:,3)'];
        lev = [rightside(:,4)'];
        edges = [rightside(:,5)'];
        r_vect = [rightside(:,6)'];
        %nodes = [x_val' y_val' z_val' lev' edges' r_vect'];
        pause(0.001);
        xlabel('x');
        ylabel('y');
        zlabel('z');
    end
    nodes = [x_val' y_val' z_val' lev' edges' r_vect'];
    %*(0.5 + 1*rand())
end