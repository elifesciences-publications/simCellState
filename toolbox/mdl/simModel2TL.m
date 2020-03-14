function [T, L, Nstate] = simModel2TL(simModel)
%  [T, L, Nstate] = simModel2TL(simModel)
% This function is a wrapper to generate the sparse matrix (deterministic
% model) compatible with the inputRead function outputs.
% It is assumed that each colum can sum at max 2 (corresponding to cell
% splitting option). All rates are qual to 1 s^-1.
%  simModel     IN: model info
%  inputT      OUT: Transitions
%  inputS      OUT: Splitting
%  Nstate      OUT: Number of state
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

% parameters
Nstate = simModel.Nstate;
indxState = simModel.indxState;
deltaState = simModel.deltaState;
sRate = simModel.sRate;

% initialization
T = []; L = [];
% addDS = zeros(Nstate); addDS(sub2ind([Nstate Nstate], indxState, indxState)) = 1;
% deltaState = deltaState + addDS;
for irate = 1:length(sRate)
    iIn = indxState(irate);
    dsti = deltaState(:,irate);
    dsti(iIn) = dsti(iIn)+1;
    iOut = find(dsti);
    if length(iOut) == 2
        L = cat(1, L, [sRate(irate) iIn iOut']);
    else
        T = cat(1, T, [sRate(irate) iIn iOut]);
    end
end

end