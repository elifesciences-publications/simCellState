function out = ppro_simPaper(simModel, simOptions, simOut, maxExt, xtime)

% This function process the generic network clonal dynamics simulations.
% The output structure contains:
% - indexes of the renewing, committed and death states
% - time and mean number of cells (both all and surviving clones) for all
%   cells, cells in the critical SCC, and in the supercritical SCCs
% - indexes and corresponding extinction rate (matching with maxExt)
% - clone size at the above time points
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

% inputs
Nstate = simModel.inputRaw.Nstate;
sysProp = sysProperties('TL', simModel.inputRaw.inputT, simModel.inputRaw.inputS, Nstate);
indxD = setdiff(1:1:Nstate, simModel.indxActiveState);
indxR = sysProp.sccIndx{cellfun(@(x) any(x == 1), sysProp.sccIndx) == 1};
indxC = setdiff(1:1:Nstate, [indxR indxD]);
xstate = simOut{1}.xstate;
time_sim = simOptions.time;

% check for NaN or extingued case
xxtot = squeeze(sum(xstate([indxR indxC],:,:),1));
nNaN = sum(isnan(xxtot), 1);
nExt = sum(xxtot == 0, 1);
inan = find(nNaN == 0, 1, 'last');
if isempty(maxExt) % use time
    isurv = NaN(1,length(xtime));
    for itime = 1:length(xtime)
        ie = find(time_sim>=xtime(itime), 1, 'first');
        if ~isempty(ie)
            isurv(itime) = ie;
        end
    end
else % extinction levels
    isurv = NaN(1,length(maxExt));
    for iext = 1:length(maxExt)
        ie = find(nExt/simOptions.Nrun>maxExt(iext), 1, 'first');
        if ~isempty(ie)
            isurv(iext) = ie;
        end
    end
end
isurv(isurv>inan | isnan(isurv)) = [];
indxT = unique([isurv inan]);
extRate = nExt(indxT)/simOptions.Nrun;

% get outputs
nn_sim = squeeze(sum(xstate([indxR indxC],:,:),1));
nn1_sim = squeeze(sum(xstate(indxR,:,:),1));
nn2_sim = squeeze(sum(xstate(indxC,:,:),1));

% Subcritical SCCs
% mean all
mean_n_all_sim = mean(nn_sim);
mean_n1_all_sim = mean(nn1_sim);
mean_n2_all_sim = mean(nn2_sim);
% mean (surviving only)
mean_n_sim = NaN(size(time_sim));
mean_n1_sim = NaN(size(time_sim));
mean_n2_sim = NaN(size(time_sim));
for itime = 1:max(indxT)
    iok = find(xxtot(:,itime)>0);
    mean_n_sim(itime) = mean(nn_sim(iok, itime));
    mean_n1_sim(itime) = mean(nn1_sim(iok, itime));
    mean_n2_sim(itime) = mean(nn2_sim(iok, itime));
end

% distribution
nn = cell(1,length(indxT));
for itime = 1:length(indxT)
    nn{itime} = nn_sim(xxtot(:,indxT(itime))>0,indxT(itime));
end
% output
out = struct('indxR', indxR, 'indxC', indxC, 'indxD', indxD, ...
    'time_sim', time_sim, 'mean_n_sim', mean_n_sim, 'mean_n_all_sim', mean_n_all_sim, ...
    'mean_n1_sim', mean_n1_sim, 'mean_n1_all_sim', mean_n1_all_sim, ...
    'mean_n2_sim', mean_n2_sim, 'mean_n2_all_sim', mean_n2_all_sim, ...
    'indxT', indxT, 'iTime', time_sim(indxT), 'extRate', extRate, 'inan', inan);
out.nn_sim = nn;

end