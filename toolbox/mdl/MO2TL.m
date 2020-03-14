function [T, L, N] = MO2TL(M, O)
% This function buld the T and L matrices starting from M and O.
% Assumptions:
% - cell division of the type ii --> ii+ii are not possible
% - in O only state transitions
% - in M just one option (transition or division) for each state
% - no multiple transitions in O and between M and O
% 
% M     IN: Deterministic matrix (N,N)
% O     IN: Stochastic matrix (N,N)
% T     OUT: Transition options. The first column contains the rates, the
%           second/third the incoming/outcoming states (M,3)
% L     OUT: Division options. The first column contains the rates, the
%           second the incoming states and the third and fourth the outcoming states (L,4)
% N     OUT: total number of states

% number of states
N = size(M,1); 

% identification of the stochastic transitions
[iO, jO] = find(O);
nO = length(iO);

% identification of the deterministic transitions and division
sumM = sum(M > 0);
[iMt, jMt] = find(M & repmat(sumM == 1, N, 1));
nMt = length(iMt);
[iMd, jMd] = find(M & repmat(sumM == 2, N, 1));
nNd = length(iMd);
if any(sumM > 2)
    error('simCellState:MO2TL:WrongInputs', 'More than one option (transition or division) in input matrix M')
end

% transitions
nT = nO + nMt;
T = zeros(nT, 3);
indxO = sub2ind([N, N], iO, jO);
indxM = sub2ind([N, N], iMt, jMt);
if ~isempty(intersect(indxO, indxM))
    error('simCellState:MO2TL:WrongInputs', 'Repeated transition in input matrix M and O')
end
T(:,1) = [O(indxO); M(indxM)];
T(:,2) = [jO; jMt];
T(:,3) = [iO; iMt];

% division
nL = nNd/2;
L = zeros(nL, 4);
L(:,1) = M(sub2ind([N, N], iMd(1:2:end), jMd(1:2:end)));
L(:,2) = jMd(1:2:end);
L(:,3) = iMd(1:2:end);
L(:,4) = iMd(2:2:end);

% checks

end