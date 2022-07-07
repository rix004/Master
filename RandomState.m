function[x_rand,y_rand]=RandomState(D)
    x_rand = D(1) + (D(2)-D(1)).*rand(1,1);
    y_rand = D(3) + (D(4)-D(3)).*rand(1,1);
end