function [eigA, eigJ] = eigSCC(A, J, S, C)
%  [eigA, eigJ] = eigSCC(A, J, S, C)
% This function returns the eigenvalues of the adiacency and jacobian
% matrices for each strongly connected component.
%
% A     IN: Adjacency (N,N)
% J     IN: Jacobian (N,N)
% S     IN: number of strongly connected components
% C     IN: id of the strongly connected components for each state (1,N)
% eigA  OUT: eigenvalues of A for each strongly connected component, cell(1,S)
% eigJ  OUT: eigenvalues of J for each strongly connected component, cell(1,S)

% eigenvalues
eigA = cell(1, S);
eigJ = cell(1, S);
for isi = 1:S
    idinode = find(C == isi);
    eigA{isi} = eig(A(idinode, idinode));
    eigJ{isi} = eig(J(idinode, idinode));
end

end