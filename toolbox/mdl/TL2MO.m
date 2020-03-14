function [M, O] = TL2MO(T, L, N)
% This function buld the M and O matrices starting from T and L.
% Assumptions:
% - in T only state transitions, no repeated
% - in L division in two states
% - cell division of the type ii --> ii+ii are not possible
% - only 1 possible division for each state
% - If for one state there is a single option, this goes in M
% - If there are more options: in M goes the divisions and in O the
%   transitions.
% 
% T     IN: Transition options. The first column contains the rates, the
%           second/third the incoming/outcoming states (M,3)
% L     IN: Division options. The first column contains the rates, the
%           second the incoming states and the third and fourth the outcoming states (L,4)
% N     IN: total number of states
% M     OUT: Deterministic matrix (N,N)
% O     OUT: Stochastic matrix (N,N)

% initialization
M = zeros(N,N);
O = zeros(N,N);

% Divisions go in M
if ~isempty(L)
    if length(unique(L(:,2))) < size(L, 1)
        error('simCellState:TL2MO:WrongInputs', 'Repeated divisions in input matrix L')
    end
    indxM = sub2ind([N N], [L(:,3); L(:,4)], [L(:,2); L(:,2)]);
    M(indxM) = [L(:,1); L(:,1)];
end

% Transitions
if ~isempty(T)
    for ii = unique(T(:,2))'
        indxTii = find(T(:,2) == ii);
        indx = sub2ind([N N], T(indxTii,3), T(indxTii,2));
        if length(indxTii) > 1 || (~isempty(L) && ~isempty(intersect(L(:,2), T(indxTii,2)))) % goes in T
            O(indx) = T(indxTii,1);
        else % goes in M
            M(indx) = T(indxTii,1);
        end
    end
end

end