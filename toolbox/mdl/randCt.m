function Ct = randCt(N, nTmax)
%  Ct = randCt(N)
% This function build the cell transition matrix Ct based on the following
% assumptions: 
% - the sum of each column can be 0 to nTmax with the same probability
% - transition of the type ii --> ii and repeated transitions are neglected
% 
%  N   IN: number of states
%  nTmax IN: maximum number of transition for each state
%  Ct  OUT: cell transition matrix

% initialization
Ct = zeros(N);
% sum of column
nnCt = randi(nTmax+1, 1, N)-1; % 0 to nTmax
for idi = 1:N
    idoCt = setdiff(unique(randi(N, 1, nnCt(idi))), idi); % remove idi and repeated transitions
    Ct(idoCt, idi) = 1;
end

end