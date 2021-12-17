function [M_flip, M_id] = extract_diags(M)
% Extract the diagonals from the banded matrix M and store in a format 
% to be used by banded_dot_star.m
[M_diags_trans, M_id] = spdiags(M');
M_flip = fliplr(M_diags_trans);

end

