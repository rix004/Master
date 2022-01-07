function [nodes] = newtree(x1,y1,theta, distance,levels)
    x2 = x1+cosd(theta)*distance
    y2 = y1+sind(theta)*distance
    %nodes_init_x=[nodes(1,:) x2];
    %nodes_init_y=[nodes(2,:) y2];
    %nodes = [nodes_init_x;nodes_init_y]
    if levels ~= 1
        line([x1 x2],[y1,y2],'LineWidth',2);
        newtree(x2,y2,theta+20,distance*0.5,levels-1);
        newtree(x2,y2,theta-20,distance*0.5,levels-1);
    end
end