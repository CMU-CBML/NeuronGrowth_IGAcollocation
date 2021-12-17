function [phi_sum] = sum_filter(phi,method)
%#codegen
% THis function takes phi as input and output a phi_sum variable that has
% maximum value at tips (0~1)

phi = full(phi);
% get size of input
[Nx,Ny] = size(phi);
% round phi -> discrete
phi = round(phi);
% initialize
phi_sum = zeros(Nx,Ny);

L = bwconncomp(phi,4);
ID = zeros(size(phi));
for i = 1:L.NumObjects
    ID(L.PixelIdxList{i}) = i;  %for regular code
%     ID(L.RegionIndices(i)) = i; %for mex
end

% loop through and calculate sum of phi values around i,j
for i = 5:Nx-4
    for k = 5:Ny-4
        baseState = ID(i,k);
        for j = k-4:k+4
            s_4 = isequal(baseState,ID(i-4,j));
            s_3 = isequal(baseState,ID(i-3,j));
            s_2 = isequal(baseState,ID(i-2,j));
            s_1 = isequal(baseState,ID(i-1,j));
            s_0 = isequal(baseState,ID(i-0,j));
            s1_ = isequal(baseState,ID(i+1,j));
            s2_ = isequal(baseState,ID(i+2,j));
            s3_ = isequal(baseState,ID(i+3,j));
            s4_ = isequal(baseState,ID(i+4,j));

            phi_sum(i,k) = phi_sum(i,k) + sum(s_4*phi(i-4,j) + s_3*phi(i-3,j) + ...
                s_2*phi(i-2,j) + s_1*phi(i-1,j) + s_0*phi(i,j) + s1_*phi(i+1,j) + ...
                s2_*phi(i+2,j) + s3_*phi(i+3,j) + s4_*phi(i+4,j));
        end
    end
end

%% remember to simplify (using imlocalmin, then this section is not needed)
% scaling for better identification
phi_sum = phi./phi_sum;
phi_sum_max = max(max(phi_sum));
phi_sum = phi_sum./phi_sum_max;

if method == 1
    phi_sum_temp = phi_sum;
    TF = isoutlier(phi_sum_temp);
    phi_sum_temp(TF) = 0;
    phi_sum_temp(isnan(phi_sum_temp))=0;
    cutoff = prctile(reshape(phi_sum_temp,(Nx)^2,1),99.5);
    
%     phi_sum(isnan(phi_sum))=0;
%     cutoff = 0.7;
elseif method == 0
    phi_sum_temp = phi_sum;
    TF = isoutlier(phi_sum_temp);
    phi_sum_temp(TF) = 0;
    phi_sum_temp(isnan(phi_sum_temp))=0;
    cutoff = prctile(reshape(phi_sum_temp,(Nx)^2,1),99.97);
else
    cutoff = 1; % dummy value for mex code (mex gives error when variable is not fully defined on all path
end
phi_sum(isnan(phi_sum))=0;
phi_sum(phi_sum<cutoff)=0;