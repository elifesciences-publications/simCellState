function M = randM(N)
%  M = randM(N)
% This function build the deterministic matrix M based on the following
% assumptions: 
% - division of the type ii --> jj+jj is not possible
% - the sum of each column can be 0, 1 or 2 with the same probability
% 
%  N   IN: number of states
%  M  OUT: deterministic matrix

% initialization
M = zeros(N);
% sum of column
nnM = randi(3, 1, N)-1; % 0 to 2
for idi = 1:N
    idoM = randi(N, 1, nnM(idi));
    while length(unique(idoM)) < length(idoM) % ii --> jj+jj is not possible
        idoM = randi(N, 1, nnM(idi));
    end
    M(idoM, idi) = 1;
end

end