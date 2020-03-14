function [nn1, nnpmf] = getPmf(nn_sim, dn)

% This function estimate the pmf of the nn_sim data.
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

% histcounts
[nnc, nncedges] = histcounts(nn_sim, 'BinLimits', [0-dn/2 ceil(max(nn_sim)/dn)*dn+dn/2], 'BinWidth', dn, 'Normalization', 'count');

% normalize
nn1 = (nncedges(1:end-1)+diff(nncedges)/2);
nnpmf = nnc/length(nn_sim);

end
