function [nodes] = tree(x_val,y_val,lev,edges,r_vect,x1,y1,theta,distance,levels,root_radius,r_rate,d_rate,delta_theta)
    if levels ~= 1
        x2 = x1+cosd(theta)*distance;
        y2 = y1+sind(theta)*distance;
        dist = sqrt((x2-x1)^2+(y2-y1)^2);
        x_val = [x_val x2];
        y_val = [y_val y2];
        lev = [lev levels-1];
        edges = [edges dist];
        r_vect = [r_vect root_radius];
        nodes = [x_val' y_val' lev' edges' r_vect'];       % lev er en kolonne som holder styr på hvilket nivå noden er på.
        line([x1 x2],[y1,y2],'LineWidth',root_radius*20,'Color',[0.8500, 0.3250, 0.0980]);
        hold on
        
        leftside = tree(x_val,y_val,lev,edges,r_vect,x2,y2,theta-delta_theta,distance*d_rate,levels-1,root_radius*r_rate,r_rate,d_rate,delta_theta);
        
        x = leftside(:,1)';
        y = leftside(:,2)';
        z = leftside(:,3)';
        e = leftside(:,4)';
        r = leftside(:,5)';
        
        rightside = tree(x,y,z,e,r,x2,y2,theta+delta_theta,distance*d_rate,levels-1,root_radius*r_rate,r_rate,d_rate,delta_theta);
        
        x_val = [rightside(:,1)'];
        y_val = [rightside(:,2)'];
        lev = [rightside(:,3)'];
        edges = [rightside(:,4)'];
        r_vect = [rightside(:,5)'];
        nodes = [x_val' y_val' lev' edges' r_vect'];
        pause(0.001);
        xlabel('x');
        ylabel('y');
    end
    nodes = [x_val' y_val' lev' edges' r_vect'];
    %*(0.5 + 1*rand())
end