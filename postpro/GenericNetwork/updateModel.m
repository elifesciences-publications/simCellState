function [p, q, inputRaw1] = updateModel(x, inputRaw, indxT, indxS, simOptions)

% inputRaw1 = inputRaw;
% for ii = 1:length(indxR1)
%     irow = find(inputRaw.inputS(:,2) == indxR1(ii));
%     inputRaw1.inputS(irow, 1) = inputRaw.inputS(irow, 1)*x(ii);
% end
% % update param
% simModel = inputSetup(inputRaw1);
% [~, ~, ~, p, q] = findEquivalentParam(simModel, simOptions);
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

inputRaw1 = inputRaw;
for ii = 1:length(indxS)
    inputRaw1.inputS(indxS(ii), 1) = inputRaw.inputS(indxS(ii), 1)*x(ii);
end
for ii = 1:length(indxT)
    inputRaw1.inputT(indxT(ii), 1) = inputRaw.inputT(indxT(ii), 1)*x(ii+length(indxS));
end

% update param
simModel = inputSetup(inputRaw1);
[~, ~, ~, p, q] = findEquivalentParam(simModel, simOptions);