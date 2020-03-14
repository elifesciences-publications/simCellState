function [xx1, xxpdf] = getPdf(nn_sim, dx)

% This function estimate the continuous normalized pdf based on the
% cumulative count of the input nn_sim.
% Output is given at steps dx.
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

% cumulative count
[nncc, nnccedges] = histcounts(nn_sim, 'BinLimits', [0 max(nn_sim)], 'BinWidth', 1,'Normalization', 'cumcount');

% normalize
xx0 = (nnccedges(1:end-1)+diff(nnccedges)/2)/mean(nn_sim);
xx1 = xx0:dx:max(xx0); 
nncdf1 = interp1(xx0, nncc/max(nncc), xx1);

% pdf (based on finite difference 2nd order)
xxpdf = varDer2(xx1, nncdf1);
xxpdf(xxpdf<0) = 0; % numerical errors

end
