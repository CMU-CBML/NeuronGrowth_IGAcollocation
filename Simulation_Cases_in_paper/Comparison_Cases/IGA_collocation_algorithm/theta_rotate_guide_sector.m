function [theta] = theta_rotate_guide_sector(lenu,lenv,max_x,max_y,rotate)
theta = zeros(lenu,lenv);

for i = 1:lenu
    for j = 1:lenv
        x=i-max_y;
        y=j-max_x;
        r = sqrt(x^2+y^2);
        if (r<10) % size of guiding sector
            theta(i,j) = (atan2(y,x));
            if isnan(theta(i,j))
                theta(i,j) = 1;
            end
            theta(i,j) = theta(i,j)+rotate;
            if(theta(i,j)>pi)
                theta(i,j) = theta(i,j) - 2*pi;
            elseif (theta(i,j)<-pi)
                theta(i,j) = theta(i,j) + 2*pi;
            end
        end
    end
end

for i = max_y-2:max_y+2
    for j = max_x-2:max_x+2
            theta(i,j) = 3;
    end
end

theta(abs(theta)<=2.8) = 0;
theta(abs(theta)>2.8) = 1;
