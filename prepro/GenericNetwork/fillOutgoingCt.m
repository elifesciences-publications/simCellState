function Ct = fillOutgoingCt(Ct, jj, kk, Ctii)
%  Ct = fillOutgoingCt(Ct, jj, kk, Ctii)
% This function fill randomly the outgoing components of the transition
% matrix Ctii (ii-SCC).
% 
% Ct   IN/OUT: Transition matrix
% jj   IN: index of the Ct rows (incoming nodes) that can be selected (1, M)
% kk   IN: Index of the Ct columns associated to Ctii (outgoing nodes) (1, N)
% Ctii IN: Transition matrix of the SCC under consideration (N+1, N+1)
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

% Number of states
Nii = size(Ctii,1)-1;

% number of outgoing states
Ctoii = Ctii(end,1:Nii);

% loop on number of states
for mm = 1:Nii
    if Ctoii(mm) > 0
        indx = randi(length(jj), 1);
        Ct(jj(indx), kk(mm)) = Ct(jj(indx), kk(mm)) + Ctoii(mm);
    end
end

return