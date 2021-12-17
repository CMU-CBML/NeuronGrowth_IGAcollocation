function [theta] = highlightZone(lenu,lenv,Max_x,Max_y,size_Max)
theta = zeros(lenu,lenv);
ttt = zeros(lenu,lenv);
for l=1:size_Max
    max_x = Max_x(l);
    max_y = Max_y(l);
    for i = max_y-2:max_y+2
        for j = max_x-2:max_x+2
                ttt(i,j) = 1;
        end
    end
    theta = theta+ttt;    
end
theta(abs(theta)>0) = 1;