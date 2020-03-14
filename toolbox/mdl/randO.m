function O = randO(N, nTmax)
%  O = randM(N)
% This function build the stochastic matrix O based on the following
% assumptions: 
% - only state transitions
% - the sum of each column can be 0 to N with the same probability
% - transition of the type ii --> ii and repeated transitions are neglected
% 
%  N   IN: number of states
%  nTmax IN: maximum number of transition for each state
%  O  OUT: stochastic matrix

% initialization
O = zeros(N);
% sum of column
nnM = randi(nTmax+1, 1, N)-1; % 0 to nTmax
for idi = 1:N
    idoM = setdiff(unique(randi(N, 1, nnM(idi))), idi); % remove idi and repeated transitions
    O(idoM, idi) = 1;
end

end