function [A, J] = MO2AJ(M, O)
% This function buld the A and J matrices starting from M and O.
% Assumptions:
% - cell division of the type ii --> jj+jj are not possible
% - in O only state transitions
% - in M just one option (transition or division) for each state
% - no multiple transitions in O and between M and O
% 
% M     IN: Deterministic matrix (N,N)
% O     IN: Stochastic matrix (N,N)
% A     OUT: Adjacency (N,N)
% J     OUT: Jacobian matrix (N,N)

% number of states
N = size(M,1); 

% Adjacency matrix
A = M + O;

% Jacobian
sumM = sum(M > 0); sumM(sumM == 0) = 1;
if any(sumM > 2)
    error('simCellState:MO2TL:WrongInputs', 'More than one option (transition or division) in input matrix M')
end
D = diag(sum(O) + sum(M./repmat(sumM, N, 1)));
J = A-D;

end