function [nodes] = MakeTree1(ls,x1,y1,val)
axis(gca,'equal')
    if val.Levels ~= 1
        x2 = x1+cosd(val.Angle)*val.Length;
        y2 = y1+sind(val.Angle)*val.Length;
        dist = sqrt((x2-x1)^2+(y2-y1)^2);
        ls.xval = [ls.xval x2]
        ls.yval = [ls.yval y2]
        ls.lev = [ls.lev val.Levels-1]
        ls.edges = [ls.edges dist]
        ls.radius = [ls.radius val.Radius]
        nodes = [ls.xval' ls.yval' ls.lev' ls.edges' ls.radius'];       % lev er en kolonne som holder styr på hvilket nivå noden er på.
        line([x1 x2],[y1,y2],'LineWidth',val.Radius*20,'Color',[0.8500, 0.3250, 0.0980]);
        hold on
        
        val.Angle=val.Angle - val.RotationAngle;
        val.Trunk = val.Length*val.LengthRate;
        val.levels = val.Levels-1;
        val.Radius = val.Radius*val.RadiusRate;
        
        leftside = MakeTree1(ls,x2,y2,val);
        
        rs.xval = leftside(:,1)';
        rs.yval = leftside(:,2)';
        rs.lev = leftside(:,3)';
        rs.edges = leftside(:,4)';
        rs.radius = leftside(:,5)';
        
        val.Angle=val.Angle + val.RotationAngle;
        val.Length = val.Length*val.LengthRate;
        val.levels = val.Levels-1;
        val.Radius = val.Radius*val.RadiusRate;
        
        rightside = MakeTree1(rs,x2,y2,val);
        
        ls.xval = [rightside(:,1)'];
        ls.yval = [rightside(:,2)'];
        ls.lev = [rightside(:,3)'];
        ls.edges = [rightside(:,4)'];
        ls.radius = [rightside(:,5)'];
        nodes = [ls.xval' ls.yval' ls.lev' ls.edges' ls.radius'];
        %pause(0.001);
        xlabel('x');
        ylabel('y');
    end
    nodes = [ls.xval' ls.yval' ls.lev' ls.edges' ls.radius'];
    %*(0.5 + 1*rand())
end