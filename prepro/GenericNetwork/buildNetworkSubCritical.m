function [Cd, Ct] = buildNetworkSubCritical(xOut_sub)
%  [Cd, Ct] = buildNetworkSubCritical(xOut_sub)
% This function build the condensed network and fill it with the
% subcritical SCC.
% The diagonal blocks are filled first. Then for each outgoing component
% the incoming one is selected randomly among the nodes of the remaining
% SCC (but only considering the lower triangular part to avoid the
% formation of a macro-SCC).
% 
% xOut_sub  IN: cell containing the structure with the SCC outputs
% Cd    OUT: Division matrix
% Ct    OUT: Transition matrix
% 
% Cristina Parigini, 14/03/2020
% 
% Copyright 2020 Cristina Parigini
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

% number of SCC
Nscc_sub = length(xOut_sub);

% state of subcricitcal SCC
% Nstate_sub = [arrayfun(@(x) size(x.Cd,1)-1, xOut_sub) 1]; % add 1 death at the end
Nstate_sub = [cellfun(@(x) size(x.Cd,1)-1, xOut_sub) 1]; % add 1 death at the end
sumNstate_sub = sum(Nstate_sub);

% initialization
Cd = zeros(sumNstate_sub, sumNstate_sub);
Ct = zeros(sumNstate_sub, sumNstate_sub);

% insert first the diagonal blocks
kk = 0;
for ii = 1:Nscc_sub
    Nii = Nstate_sub(ii);
    kk = kk(end)+(1:Nii);
    Cd(kk, kk) = xOut_sub{ii}.Cd(1:Nii, 1:Nii);
    Ct(kk, kk) = xOut_sub{ii}.Ct(1:Nii, 1:Nii);
end
% now fill (randomly) the lower triangular part
kk = 0; jj = 1;
for ii = 1:Nscc_sub
    Nii = Nstate_sub(ii);
    kk = kk(end)+(1:Nii);
    if ii == Nscc_sub
        jj = jj(1)+Nii:sumNstate_sub;
    else
        jj = jj(1)+Nii:sumNstate_sub-1; % not to death
    end
    % assign randomly the outgoing components for the cell division
    Cd = fillOutgoingCd(Cd, jj, kk, xOut_sub{ii}.Cd);
    % assign randomly the outgoing components for the cell transition
    Ct = fillOutgoingCt(Ct, jj, kk, xOut_sub{ii}.Ct);
end

return