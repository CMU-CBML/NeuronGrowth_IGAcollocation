function [A] = banded_dot_star(t, M_flip, M_id)
%#codegen
%computes the dot star product A = t.*M of a vector t and a banded sparse
% matrix M by extracting the diagonals of M
n = size(M_flip,1);
M_t_diags = bsxfun(@times,t,M_flip);
A = spdiags(fliplr(M_t_diags), M_id, n, n)';
end

