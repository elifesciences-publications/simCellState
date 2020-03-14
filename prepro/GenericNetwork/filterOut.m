function xOut = filterOut(xOut, icons)

% remove those networks that:
% - do not have any outgoing component (check type = -2)
% - are not a SCC (type = NaN)
% - conserved (optionally)
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

% filter output
type = cellfun(@(x) x, xOut(:,4));
eigM = cellfun(@(x) x, xOut(:,3));
irm = type == -2 | ... % removing SCC w/o outgoing components
    isnan(type) | ... % removing not SCC
    (type == -1 & eigM > 0) | ... % remove SCC that cannot be stabilized and are unstable
    (icons == 0 & type == 0); % remove conserved if icons = 0
xOut(irm,:) = [];

end