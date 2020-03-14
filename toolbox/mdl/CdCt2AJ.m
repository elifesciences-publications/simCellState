function [A, J] = CdCt2AJ(Cd, Ct)
%  [A, J] = CdCt2AJ(Cd, Ct)
% This function buld the adiacency and jacoian matrices starting from the
% division and transition matrices.
% Assumptions:
% - in Ct all the state transitions
% - in Cd only cell division, ii --> jj+jj is also possible
% - only one cell division option per state
% - each element is the rate
% 
% Cd    IN: Division matrix (N,N)
% Ct    IN: Transition matrix (N,N)
% A     OUT: Adjacency (N,N)
% J     OUT: Jacobian matrix (N,N)

% Adjacency matrix
A = Cd + Ct;

% Jacobian
D = diag(sum(Ct) + sum(Cd/2));
J = A-D;

end