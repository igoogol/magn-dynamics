function [H1] = generate_block_diag(A,M)
% #========================================================================
% # A; M integers 
% # generate block diagonal matrix with diagonal element matrix A. In total
% # M units of A-blocks.
% #     H0 = {A  0  0  0 ...
% #           0  A  0  0 ...
% #           0  0  A  0 ...
% #           0  0  0  A ...
% #           ...  ...   ...}_{M*M} 
% #========================================================================
for k = 1:M
  AA{k} = A;
end
H1 = blkdiag(AA{:});