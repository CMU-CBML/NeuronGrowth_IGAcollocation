function [theta_ori] = cuedGrowthStages3(phi_plot,CX,CY)
%#codegen
    [lenu,lenv] = size(phi_plot);
    phi_id = round(phi_plot);
    tip = sum_filter(phi_id,1);
    regionalMaxima = imregionalmax(full(tip));
    L = bwconncomp(regionalMaxima,4);
    S = regionprops(L,'Centroid');

    centroids = floor(cat(1,S.Centroid));
%     % Matlab MEX cannot directly access S.centroid, need to loop
%     numRegions = numel(S);
%     centroids = zeros(numRegions, numel(S(1).Centroid)); % pre-initialize to set size
%     for region = 1:numRegions
%         centroids(region, :) = floor(S(region).Centroid); % write one row at a time
%     end
    
    max_x = [];
    max_y = [];
    numCue = length(CX);
    dist= zeros(lenu,lenv);
    for k = 1:numCue
        for i = 1:lenu
            for j = 1:lenv
                dist(i,j) = sqrt(bsxfun(@times,(i-CX(k)),(i-CX(k)))+bsxfun(@times,(j-CY(k)),(j-CY(k))));
            end
        end

        tips = zeros(1,L.NumObjects);
        for i = 1:L.NumObjects
            tips(i) = dist(centroids(i,1),centroids(i,2));
        end
        [~,tip_want_ind] = max(tips);
%         tip_want_ind = find(tips>=0.95*max(tips));
        max_x = [max_x,centroids(tip_want_ind,1).'];
        max_y = [max_y,centroids(tip_want_ind,2).'];
    end

    size_Max = length(max_x);
    [theta_ori] = theta_rotate(lenu,lenv,max_x,max_y,size_Max);

end