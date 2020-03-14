function [S, C, sSCC] = findSCC(R)
%  [S, C, sSCC] = findSCC(R)
% This function identifies the strongly connected components.
% 
% R     IN: sparse of the transpose of the Adjacency (N,N)
% S     OUT: number of strongly connected components
% C     OUT: strongly connected component number associated to each state (1,N)
% sSCC  OUT: size of each strongly connected component, (1,S)

% strongly connected components
[S, C] = graphconncomp(R);

% remove the single nodes that are not ssc
sSCC = zeros(1, S);
for isi = 1:S
    idinode = find(C == isi);
    sizeSCC_i = length(idinode);
    if sizeSCC_i > 1 || ~(sizeSCC_i == 1 && R(idinode, idinode) == 0)
        sSCC(isi) = sizeSCC_i;
    else
        C(idinode) = 0;
    end
end
S = length(setdiff(unique(C), 0));
indxscc = find(sSCC);
for isi = 1:S
    C(C==indxscc(isi)) = isi;
end
sSCC(sSCC == 0) = [];

end