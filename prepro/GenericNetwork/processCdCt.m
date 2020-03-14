function xOut = processCdCt(Cdii, Ctjj)

% This function returns a cell containing:
% - Cd and Ct matrices, xOut{1:2}
% - maximum eigenvalue, xOut{3}
% - type of SCC, xOut{4}:
%      NaN not a SCC
%      0 SCC, stable for any rate
%      1 SCC, can be stabilized,
%     -1 SCC, cannot be stabilized
%     -2 SCC, nothing is going outside the SCC
% If there is nothing going outsude the SCC, then 
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

N = size(Cdii, 1)-1; % number of states
xOut = cell(1, 4); % initialization: Cd, Ct, eigMax, Type
[A, J] = CdCt2AJ(Cdii, Ctjj);
R = sparse(A');
[S, C] = findSCC(R);
xOut{1} = Cdii; xOut{2} = Ctjj;
eg = eig(J(1:N,1:N)); [~, im] = max(real(eg));
xOut{3} = eg(im);
% assinging type of SCC
if any(Cdii(end,:)>0) || any(Ctjj(end,:)>0) % at least something going outside the SCC
    if S == 1 && all(C(1:N) == 1) % 1 SCC with all the nodes
        [~, ~, type] = stabilityLogic(Cdii, Ctjj);
        xOut{4} = type;
    elseif N == 1 % not exactly a SCC (single node)
        xOut{4} = -1;
    else
        xOut{4} = NaN;
    end
else
    xOut{4} = -2; % nothing going outside the SCC
end