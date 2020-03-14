function [Ns, Nl, type, sysProp] = stabilityLogic(Cd, Ct)
%  [Ns, Nl, type, sysProp] = stabilityLogic(Cd, Ct)
% This function implements the logic for determining if a system can be
% stabilized or not. The logic is as follow:
% - Ns = Nl = 0 -> this is just a stable oscillatory case, whatever the parameters are.
% - Ns > 0 and Nl > 0 -> if all the conserved cycles have an outgoing component the system can be stabilized.
% - In the remaining cases nothing can be done.
% System with only one SCC are analyzed. This has to be generalized.
%
%  Cd    IN: Division matrix (N,N)
%  Ct    IN: Transition matrix (N,N)
%  Ns   OUT: number of source nodes
%  Nl   OUT: number of leak node
%  type OUT: Type of system: 0 stable for any rate; 1: can be stabilized,
%           -1 cannot be stabilized.
%  sysProp OUT: system properties structure
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

% system properties
sysProp = sysProperties('CdCt', Cd, Ct);
NSCC = sysProp.NSSC;
% initialization
Ns = -1;
type = NaN;
Nl = -1;
if NSCC == 1
    [Ns, Nl, type] = stabilityLogic_SCCi(sysProp, 1, Cd, Ct);
elseif NSCC > 1 % only if all are decayng, this must be generalized
    idec = cellfun(@(x) all(real(x) < 0), sysProp.eigJ);
    % igrow = cellfun(@(x) all(real(x) > 0), sysProp.eigJ);
    if all(idec) % all are decayng: it is sufficient to stablize one
        for isc = 1:NSCC
            [~, ~, type] = stabilityLogic_SCCi(sysProp, isc, Cd, Ct);
            if type == 1
                type = 2;
                sysProp.SCCstab = isc;
                break
            end
        end
    else % if
    end
end

end

% stability logic for the i-SCC
function [Ns, Nl, type] = stabilityLogic_SCCi(sysProp, isc, Cd, Ct)

Ap = sysProp.Ap;

% indentification of the leak and source nodes in the SCC
inode = find(sysProp.SCC == isc);
[indexS, indexL] = findSourceLeak(Cd, Ct, inode);
Ns = length(indexS);
Nl = length(indexL);

% is there any source/leak?
type = -1;
if Ns == 0 && Nl == 0
    type = 0; % oscillatory stable for any rate
elseif Ns*Nl > 0 % at least one source and one leak
    % find conserved cycles
    type = 1;
    inode1 = setdiff(inode, indexL);
    for ii = 1:length(inode1) % reorder here to try first the shortest and then the longest sequences (1-10, 2-9, ...)
        [flag, icc] = conservedCycle(Ap, inode1, ii, Ct);
        if flag
            type = -1;
            break
        end
    end
end
end

% conserved cycle analysis
function [flag, icc] = conservedCycle(Ap, inode, k, Ct)

% This function find a conserved cycle w/o any outgoing component.
% If the sum of the columns of the Ap matrix for a given combination of
% nodes is greater than one (conserved or growing cycle) then check for any
% outgoing component.

% loop on possible combinations
C = nchoosek(inode, k);
for ic = 1:size(C, 1)
    % node under analysis
    icc = C(ic,:);
    % condition
    flag = all(sum(Ap(icc, icc), 1) >= 1); % conserved cycle
    if flag
        if any(sum(Ct(setdiff(1:1:length(Ap), icc),icc), 1) > 0) % outgoing component
            flag = 0;
        end
    end
    if flag
        break
    end
end

end

% THIS FUNCTION IS NOT USED, but equivalent to the above
function [flag, icc] = conservedCycle_1(Ap, inode, k, Cd, Ct)

% This function find a conserved cycle w/o any outgoing component
% It is based on the expected sum of the columns of the Ap matrix

% sum for conserved cycles
s0 = ones(1, size(Ap, 1));

% loop on possible combinations
C = nchoosek(inode, k);
for ic = 1:size(C, 1)
    % node under analysis
    icc = C(ic,:);
    % in case a source
    ss = s0;
    indexS = findSourceLeak(Cd, Ct, icc);
    ss(indexS) = (2+sum(Ct(:,indexS)~=0))./(sum(Ct(:,indexS)~=0)+1);
    % condition
    flag = all(sum(Ap(icc, icc), 1) >= ss(icc));
    if flag
        break
    end
end

end
