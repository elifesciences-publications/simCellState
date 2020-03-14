clear all
close all

% This function builds a set of random networks by combining the previously
% generated SCC.
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
NN = 1e3; % total number of networks
NsccMax = 3; % max number of subcritical SCC
SD0 = 2*round(cputime*1e3)+1; rng(SD0)
eigTol = 1e-4; % tolerance for eigenvalues
eigShift = -1; % max eig expected in subcritical SCC
Rmean = 1; Rmin = 0.3; % minimum and mean rates
fileName = 'SCC_All';
A = load(fullfile('SCC', fileName));
outNC = A.xOutNC; propNC = A.xPropNC;
outNCS = A.xOutNCS; propNCS = A.xPropNCS;
outC = A.xOutC; propC = A.xPropC;
clear A

% BUILD NETWORK
% -------------
% initialization
xCdCt_nc = cell(2, NN);
xCdCt_c = cell(2, NN);
NstateR_nc = zeros(1, NN);
NstateR_c = zeros(1, NN);
xierr = cell(1, NN);
% loop
for in = 1:NN
    % build random network
    [Cd_nc, Ct_nc, Nr_nc, Cd_c, Ct_c, Nr_c, ierr] = genericNetwork(outNC, outNCS, outC, propNC, propNCS, propC, ...
        NsccMax, eigShift, eigTol, Rmean, Rmin);
    % output
    xCdCt_nc(:,in) = {Cd_nc; Ct_nc};
    xCdCt_c(:,in) = {Cd_c; Ct_c};
    NstateR_nc(in) = Nr_nc;
    NstateR_c(in) = Nr_c;
    xierr{in} = ierr;
end

% SAVE OUTPUT
outDir = fullfile('../', '../', 'io', 'IN', 'GENERIC');
tag = dateTag(clock);
fileName = 'out';
% not conserved
xCdCt = xCdCt_nc; NstateR = NstateR_nc;
fileOut = sprintf('%s_ncons_%s', fileName, tag);
mkdir(fullfile(outDir, fileOut))
save(fullfile(outDir, fileOut, fileOut), 'xCdCt', 'NstateR', 'xierr', 'SD0')
% conserved
xCdCt = xCdCt_c; NstateR = NstateR_c;
fileOut = sprintf('%s_cons_%s', fileName, tag);
mkdir(fullfile(outDir, fileOut))
save(fullfile(outDir, fileOut, fileOut), 'xCdCt', 'NstateR', 'xierr', 'SD0')