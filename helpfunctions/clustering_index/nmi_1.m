function MIhat = nmi_1(A, B)
%NMI Normalized mutual information
% A, B: 1*N;
if length(A) ~= length(B)
    error('length( A ) must == length( B)');
end
N = length(A);
A_id = unique(A);
K_A = length(A_id);
B_id = unique(B);
K_B = length(B_id);
% Mutual information
A_occur = double (repmat( A, K_A, 1) == repmat( A_id', 1, N ));
B_occur = double (repmat( B, K_B, 1) == repmat( B_id', 1, N ));
AB_occur = A_occur * B_occur';
P_A= sum(A_occur') / N;
P_B = sum(B_occur') / N;
P_AB = AB_occur / N;
MImatrix = P_AB .* log(P_AB ./(P_A' * P_B)+eps);
MI = sum(MImatrix(:));
% Entropies
H_A = -sum(P_A .* log(P_A + eps),2);
H_B= -sum(P_B .* log(P_B + eps),2);
%Normalized Mutual information
MIhat = MI / sqrt(H_A*H_B);
end