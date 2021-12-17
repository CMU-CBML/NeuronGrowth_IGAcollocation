function [theta_ori] = growthStages3(phi_id)
%#codegen
    [lenu,lenv] = size(phi_id);
    phi_id = round(phi_id);

    L = bwconncomp(phi_id,4);
    S = regionprops(L,'Centroid');

    centroids = floor(cat(1,S.Centroid));
%     % Matlab MEX cannot directly access S.centroid, need to loop
%     numRegions = numel(S);
%     centroids = zeros(numRegions, numel(S(1).Centroid)); % pre-initialize to set size
%     for region = 1:numRegions
%         centroids(region, :) = floor(S(region).Centroid); % write one row at a time
%     end
    
    ID = zeros(size(phi_id));
    dist= zeros(lenu,lenv,L.NumObjects);

    max_x = [];
    max_y = [];
    for k = 1:L.NumObjects
        ID(L.PixelIdxList{k}) = k;
%         ID(L.RegionIndices(k))= k;
        for i = 1:lenu
            for j = 1:lenv
                dist(i,j,k) = (ID(i,j) == k)*sqrt((i-centroids(k,2))^2+(j-centroids(k,1))^2);
            end
        end

        dist_k = reshape(dist(:,:,k),lenu*lenv,1);
        [~,max_index] = max(dist_k);
        max_x = [max_x,ceil(max_index/lenu)];
        max_y = [max_y,rem(max_index,lenu)];

    end
    size_Max = length(max_x);
    [theta_ori] = theta_rotate(lenu,lenv,max_x,max_y,size_Max);