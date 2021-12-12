function [X,Y] = ashleeExpDist(old_max_x,old_max_y,max_x,max_y,r)

X = [];
Y = [];
for i = 1:min(length(old_max_x),length(max_x))
    x_diff = max_x(i)-old_max_x(i);
    y_diff = max_y(i)-old_max_y(i);
    growth_angle = atan2d(y_diff,x_diff);
%     g_diff = sqrt(y_diff^2+x_diff^2);

% %     if g_diff>5
% %         g_diff = 2;
% %     end
    g_diff = 10;
    aa = growth_angle + r;

    y_increment = sind(aa)*g_diff;
    x_increment = cosd(aa)*g_diff;

    X = [X,round(max_x(i)+x_increment)];
    Y = [Y,round(max_y(i)+y_increment)];
end
