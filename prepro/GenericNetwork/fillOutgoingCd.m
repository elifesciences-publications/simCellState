function Cd = fillOutgoingCd(Cd, jj, kk, Cdii)
%  Cd = fillOutgoingCd(Cd, jj, kk, Cdii)
% This function fill randomly the outgoing components of the division
% matrix Cdii (ii-SCC).
% 
% Cd   IN/OUT: Division matrix
% jj   IN: index of the Cd rows (incoming nodes) that can be selected (1, M)
% kk   IN: Index of the Cd columns associated to Cdii (outgoing nodes) (1, N)
% Cdii IN: Division matrix of the SCC under consideration (N+1, N+1)
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
Nii = size(Cdii,1)-1; 

% number of outgoing cells (cd/total rate)
Cdoii = Cdii(end,1:Nii);
muii = sum(Cdii(:,1:Nii))/2; muii(muii==0) = 1;
Cdoiin = Cdii(end,1:Nii)./muii;

% loop on number of states
for mm = 1:Nii
    if Cdoii(mm) > 0
        indx = randi(length(jj), [1 Cdoiin(mm)]); % one or 2 depending on Cd value
        if length(indx) == 2 && length(unique(indx)) == 2 % not repeated
            Cd(jj(indx), kk(mm)) = Cd(jj(indx), kk(mm))+Cdoii(mm)/2;
        else
            Cd(jj(indx), kk(mm)) = Cd(jj(indx), kk(mm))+Cdoii(mm);
        end
    end
end

return