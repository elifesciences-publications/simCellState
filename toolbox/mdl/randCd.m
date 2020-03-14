function Cd = randCd(N)
%  Cd = randCd(N)
% This function build the matrix with cell division Cd based on the following
% assumptions: 
% - the sum of each column can be 0, or 2 respectively with 2/3 and 1/3 probability
% 
%  N   IN: number of states
%  Cd OUT: cell division matrix

% initialization
Cd = zeros(N);
% sum of column
nnCd = randi(3, 1, N)-1; % 0 to 2
nnCd(nnCd<2) = 0;
for idi = 1:N
    idoCd = randi(N, 1, nnCd(idi));    
    if length(unique(idoCd)) < 2
       Cd(idoCd, idi) = 2;
    else
        Cd(idoCd, idi) = 1;
    end
end

end