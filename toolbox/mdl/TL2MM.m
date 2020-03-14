function [MM, eigMM] = TL2MM(T, L, N)
%   [MM, eigMM] = TL2MM(T, L, N)
% This function computes the mean matrix of the process and its
% eigenvalues.
% T     IN: Transition options. The first column contains the rates, the
%           second/third the incoming/outcoming states (M,3)
% L     IN: Division options. The first column contains the rates, the
%           second the incoming states and the third and fourth the outcoming states (L,4)
% N     IN: total number of states
% MM   OUT: Mean matrix
% eigMM OUT: Eigenvalue of the mean matrix

% reproduction generating coefficients
[p0, p1, p2] = repGenFcnCoeff(T, L, N);

% mean matrix
tol = 1e-10;
rgFnc = @(s) repGenFnc(s, p0, p1, p2);
MM = numjac(@(t, x) rgFnc(x), 0, ones(N,1), rgFnc(ones(N,1)), tol, [], 0);

% eigenvalues
eigMM = eig(MM);

return
