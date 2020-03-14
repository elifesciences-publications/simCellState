function [T, L, N] = CdCt2TL(Cd, Ct)
%   [T, L, N] = CdCt2TL(Cd, Ct)
% This function buld the T and L matrices starting from Cd and Ct.
% Assumptions:
% - in Ct all the state transitions
% - in Cd only cell division, ii --> jj+jj is also possible
% - only one cell division option per state
% - each element is the rate
% 
% Cd    IN: Cell division matrix (N,N)
% Ct    IN: Cell transition matrix (N,N)
% T     OUT: Transition options. The first column contains the rates, the
%           second/third the incoming/outcoming states (M,3)
% L     OUT: Division options. The first column contains the rates, the
%           second the incoming states and the third and fourth the outcoming states (L,4)
% N     OUT: total number of states

% number of states
N = size(Cd,1); 

% identification of the transitions
[iCt, jCt] = find(Ct);
nCt = length(iCt);
T = zeros(nCt , 3);
indxCt = sub2ind([N, N], iCt, jCt);
T(:,1) = Ct(indxCt);
T(:,2) = jCt;
T(:,3) = iCt;

% identification of the division
iCd = find(sum(Cd>0));
nCd = length(iCd);
L = zeros(nCd, 4);
for icd = 1:nCd
    indx = find(Cd(:,iCd(icd)));
    if length(indx) == 1
        r = Cd(indx, iCd(icd))/2;
        indx = indx*ones(1,2);
    else
        r = Cd(indx(1), iCd(icd));
    end
    L(icd,1) = r;
    L(icd,2) = iCd(icd);
    L(icd,3) = indx(1);
    L(icd,4) = indx(2);
end

end