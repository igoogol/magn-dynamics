function [H0] = generate_block_tridiag(A,B,M)
% #========================================================================
% # A,B square matrices; M integers
% # generate block tri-diagonal matrix with diagonal element matrix A, and 
% # off-diagonal element matrix B. In total M units of A-blocks.
% #     H0 = {A  B  0  0 ...
% #           B  A  B  0 ...
% #           0  B  A  B ...
% #           0  0  B  A ...
% #           ...  ...   ...}_{M*M} 
% #========================================================================
dim = length(A);
Tmp1 = generate_block_diag(B,M-1);
Tmp2 = generate_block_diag(B',M-1);
zeros_1 = zeros(dim*(M-1),dim);
zeros_2 = zeros(dim,dim*(M-1));
zeros_3 = zeros(dim,dim);    
H1 = [[zeros_1,Tmp1]; [zeros_3,zeros_2]];   
H2 = [[zeros_3,zeros_2]; [Tmp2,zeros_1]];
H0 = generate_block_diag(A,M) + H1 + H2;