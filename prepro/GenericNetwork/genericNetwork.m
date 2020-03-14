function [Cd_nc, Ct_nc, Nr_nc, Cd_c, Ct_c, Nr_c, ierr] = genericNetwork(outNC, outNCs, outC, propNC, propNCs, propC, NsccMax, eigenMax, eigTol, Rmean, Rmin)
%  [Cd_nc, Ct_nc, Cd_c, Ct_c, ierr] = genericNetwork(outNC, outNCs, outC, NsccMax, eigShift, eigTol)
% This function build a random network starting from given strongly
% connected components (SCC). A random number of downstream SCC are always
% randomply picked from the non conserved SCC. In case these SCC are not
% subcritical, the rates are changed to met the subctritical criterion (max
% eigenvalue of the Jacobian eigenMax). The links between SCC are randomly
% setup by filling the lower diagonal matrix (see buildNetworkSubCritical
% function). A final single node representing death is added.
% Two networks are then build by attacching as upstream SCC the non
% conserved or the conserved SCC. If the non conserved SCC is not critical,
% the rates are changed to met the criteria (absolute value of the max
% eigenvalue of the Jacobian less tham eigTol). The links with remaining SCC
% are randomly setup by filling the lower diagonal matrix (see
% attachNetworkCritical function).
%
% outNC         IN: output structure containing the non conserved SCC
% outNCs        IN: output structure containing the non conserved SCC that
%                   can be stabilized
% outC          IN: output structure containing the conserved SCC
% propNC, propNCs, propC    IN: structure containing the number of transitions,
%                   divsion, size of the SCC and number of networks
% NsccMax       IN: maximum number of subcritical SCC
% eigenMax, eigTol  IN: maximum eigenvalue and tolerance for rates modification
% Rmean, Rmin   IN: Mean and minimum rate
% Cd_nc        OUT: Division matrix non conserved network
% Ct_nc        OUT: Transition matrix non conserved network
% Nr_nc        OUT: Number of states of the critical SCC (non conserved)
% Cd_c         OUT: Division matrix conserved network
% Ct_c         OUT: Transition matrix conserved network
% Nr_c         OUT: Number of states of the critical SCC (conserved)
% ierr         OUT: error flag
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

% total number of subcritical SCC (condensed structure)
Nscc_sub = randi(NsccMax);

% select subcritical SCC
xOut_sub = cell(1,Nscc_sub); ierr_sub = NaN(1,Nscc_sub);
NSCC = arrayfun(@(x) x.NSCC, propNC);
for ii = 1:Nscc_sub
    [out, indx01] = randSelection(outNC, NSCC, []);
    [out, ierr_sub(ii)] = loopModifyRate({out}, propNC(indx01), Rmean, Rmin, eigenMax, eigTol);
    xOut_sub{ii} = struct('Cd', out{1}{1}, 'Ct', out{1}{2});
end
% select critical SCC
NSCC = arrayfun(@(x) x.NSCC, propNCs);
rNSScrit = randi(length(unique(NSCC)), 1); % select SCC of the same size (assuming there are)
[out, indx01] = randSelection(outNCs, NSCC, rNSScrit);
[out, ierr_crit_nc] = loopModifyRate({out}, propNCs(indx01), Rmean, Rmin, 0, eigTol);
xOut_crit_nc = struct('Cd', out{1}{1}, 'Ct', out{1}{2});
NSCC = arrayfun(@(x) x.NSCC, propC);
out = randSelection(outC, NSCC, rNSScrit);
[out, ierr_crit_c] = loopModifyRate({out}, propNCs(indx01), Rmean, Rmin);
xOut_crit_c = struct('Cd', out{1}{1}, 'Ct', out{1}{2});

% build the network with all the subcritical SCC
[Cdsub, Ctsub] = buildNetworkSubCritical(xOut_sub);

% attach upstream the critical SCC
[Cd_nc, Ct_nc] = attachNetworkCritical(xOut_crit_nc, Cdsub, Ctsub);
Nr_nc = length(Cd_nc)-length(Cdsub);
% 2) conserved
[Cd_c, Ct_c] = attachNetworkCritical(xOut_crit_c, Cdsub, Ctsub);
Nr_c = length(Cd_c)-length(Cdsub);

% ierr
ierr = [ierr_sub, ierr_crit_nc, ierr_crit_c];
% % plot Network and check eigenvalues
% sysProp = sysProperties('CdCt', Cd_nc, Ct_nc);
% plotNetwork(sysProp, '')

end

function [out1, ierr] = loopModifyRate(out, propNC, Rmean, Rmin, varargin)

kk = 0; ierr = NaN;
while (ierr ~= 0) && (kk < 10)
    [out1, ierr] = modifyRates(out, propNC, Rmean, Rmin, varargin{:});
    kk = kk+1;
end

end

function [xout, indx01] = randSelection(out, NSCC, rNSCC)

% random selection of a SCC, uniform in size of SCC, and number of
% transition and division
if isempty(rNSCC)
    rNSCC = randi(length(unique(NSCC)), 1); % random size of SCC
end
indx0 = find(NSCC == rNSCC);
rNTD = randi(length(indx0), 1); % random number of Transition/Division
indx01 = indx0(rNTD);
indx1 = randi(size(out{indx01},1), 1); % random network
xout = out{indx01}(indx1,:);

end