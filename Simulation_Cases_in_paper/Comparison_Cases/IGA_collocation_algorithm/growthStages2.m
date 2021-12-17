function [theta_ori] = growthStages2(phi_plot)

    [lenu,lenv] = size(phi_plot);
    tip = sum_filter(phi_plot,0);
    regionalMaxima = imregionalmax(full(tip));
    [Max_y,Max_x] = find(regionalMaxima);
    size_Max = length(Max_x);
    [theta_ori] = theta_rotate(lenu,lenv,Max_x,Max_y,size_Max);
end