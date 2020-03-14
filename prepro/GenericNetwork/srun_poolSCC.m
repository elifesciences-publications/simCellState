clear all
close all

% This function generates the pools of SCC: conserved, non-conserved and
% non-conserved subcritical.
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
xNSCC = 1:4;
dOut = 'SCC';
SD0 = 501; rng(SD0)
Rmin = 0.3; Rmean = 1;
eigMaxS = -0.25; eigTgSM = 0; eigTol = 1e-6;

% initialization
xOutNC = []; xOutC = []; xPropNC = []; xPropC = [];
% load data
for ii = xNSCC
    fileName = sprintf('SCC_N%d_C.mat', ii);
    A = load(fullfile(dOut,fileName));
    xOutC = [xOutC; A.xOutF(:)];
    xPropC = [xPropC; A.xProp(:)];
    if ii == 4
        fileName = 'SCC_N4_nC_20190326161724';        
    else
        fileName = sprintf('SCC_N%d_nC.mat', ii);
    end
    A = load(fullfile(dOut,fileName));
    xOutNC = [xOutNC; A.xOutF(:)];
    xPropNC = [xPropNC; A.xProp(:)];
end
% remove empty cases/not used fields
irm = arrayfun(@(x) x.NN == 0, xPropC);
xOutC(irm) = []; xPropC(irm) = [];
xPropC = rmfield(xPropC, 'indxType1');
irm = arrayfun(@(x) x.NN == 0, xPropNC);
xOutNC(irm) = []; xPropNC(irm) = [];
% select the non conserved that can be stabilized
xOutNCS = xOutNC; xPropNCS = xPropNC;
irm = arrayfun(@(x) isempty(x.indxType1), xPropNCS);
xOutNCS(irm) = []; xPropNCS(irm) = [];
xOutNCS = cellfun(@(x,y) x(y,:), xOutNCS, {xPropNCS.indxType1}', 'UniformOutput', false);

% save to file
save(fullfile(dOut,'SCC_All'), 'xOutC', 'xOutNC', 'xOutNCS', 'xPropC', 'xPropNC', 'xPropNCS')

