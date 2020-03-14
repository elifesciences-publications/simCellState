function [Cd, Ct] = TL2CdCt(T, L, N)
% This function buld the Cd and Ct matrices starting from T and L.
% Assumptions:
% - in Ct all the state transitions
% - in Cd only cell division, ii --> jj+jj is possible
% - only one cell division possibility per state
% - each element is the rate
% 
% T     IN: Transition options. The first column contains the rates, the
%           second/third the incoming/outcoming states (M,3)
% L     IN: Division options. The first column contains the rates, the
%           second the incoming states and the third and fourth the outcoming states (L,4)
% N     IN: total number of states
% Cd    OUT: Division matrix (N,N)
% Ct    OUT: Transition matrix (N,N)

% initialization
Cd = zeros(N);
Ct = zeros(N);

% add rates
for ii = 1:size(T,1)
    Ct(T(ii,3), T(ii,2)) = T(ii,1);
end
for ii = 1:size(L,1)
    Cd(L(ii,3), L(ii,2)) = Cd(L(ii,3), L(ii,2)) + L(ii,1);
    Cd(L(ii,4), L(ii,2)) = Cd(L(ii,4), L(ii,2)) + L(ii,1);
end

end