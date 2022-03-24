function[distance]=DistanceToEdge(p,X1,X2)
            D1 = norm(p-X1);
            D2 = norm(X2-X1);
            D3 = norm(p-X2);
            theta = acosd((D1^2+D2^2-D3^2)/(2*D1*D2));
            ratio = D1*cosd(theta)/D2;
            distance = D1*sind(theta);
end