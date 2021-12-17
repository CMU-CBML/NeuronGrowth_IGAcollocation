function [initial_LAPpk_out] = ...
    kqExpandDomainForInitialLAPpk(sz,initial_LAPpk)

len = length(initial_LAPpk);
M = sqrt(len);

initial_LAPpk = reshape(initial_LAPpk,M,M);

[M, ~] = size(initial_LAPpk);
outSz = sqrt(sz);

initial_LAPpk_out = zeros(outSz);

i_off = floor(outSz/2-M/2);
j_off = floor(outSz/2-M/2);

initial_LAPpk_out(3+i_off:M-2+i_off,3+j_off:M-2+j_off) = initial_LAPpk(3:M-2,3:M-2);
initial_LAPpk_out =sparse(reshape(initial_LAPpk_out,outSz^2,1));
