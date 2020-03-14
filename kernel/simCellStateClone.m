function simCellStateClone(icase, scase, simInFolder, simCase, outDir, ii0, jrun)

% This function run multiple stochastic simulations. The external loop is
% on the cases under the input files (xCdCt); the internal on the initial
% conditions 
% if length(ii0) == 1: all the cases starting from the ii0 one are run and
%   for each case, Nstate simulations are run assuming 1 cell in each state.
% if length(ii0) > 1: only the ii0 cases are run ([ii0 ii0] run just this
%   case), starting from initial condition in jrun
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

% INPUT PARAMETERS
% ----------------
if length(ii0) == 1 % starting network, all states
    StoInput = load(fullfile(outDir, 'io', 'IN', simInFolder, icase, icase));
    NN = length(StoInput.xCdCt);
    xrun = ii0:NN;
    if isempty(jrun)
        Nstate = length(StoInput.xCdCt{1});
    elseif length(jrun) == 1
        Nstate = jrun;
        jrun = [];
    else
        Nstate = 1; % 1 sim starting from j for each i
    end
else % selected networks, one selected state
    xrun = unique(ii0);
    Nstate = 1;
end

for ii = 1:length(xrun) % loop on cases
    iirun = xrun(ii);
    for jj = 1:Nstate % loop on initial condition        
        if isempty(jrun)
            jjrun = jj;
        else
            jjrun = jrun(ii);
        end
        % SIMULATION
        [inputRaw, simOptions] = inputRead(simCase, {icase, iirun, jjrun, scase{:}}, outDir);
        if ~isempty(inputRaw)
            simModel = inputSetup(inputRaw);
            % simulation
            tic
            simCellStateLoop(simModel, simOptions);
            toc
        end
    end
end