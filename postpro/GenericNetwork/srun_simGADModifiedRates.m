clear all
close all

% This script modifies the current network parameters to match the desired 
% network parameters and run the stochastic simulations.
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

addpath('C:\Users\cp4u17\OneDrive - University of Southampton\Work\code\MOPSO')

% parameters
d0 = fullfile('..', '..', 'IO', 'OUT', 'GENERIC');
tg = '20190409150503';
outFolder = ['out_cons_' tg];
p_tg = [0.5 1 2 5 10 30]; Rmax = 150;
Nrun = 1e4; indx1 = []; 
% indx = 1:1:1e3; % all sim
% indx1 = [26 33 57]; % only specific sims
indx = 870*ones(size(p_tg)); % [94 57 30 302 870];

% initalization
if length(p_tg) == 1
    simOutFolder = fullfile(d0, sprintf('%s_p%d', outFolder, round(p_tg)));
    p_tg = p_tg*ones(size(indx));
else
    simOutFolder = fullfile(d0, sprintf('%s_indx%d', outFolder, unique(indx)));
end
if ~exist(simOutFolder, 'dir')
    mkdir(simOutFolder)
end
xout = cell(length(indx), 4);
% loop
for ii = 1:length(indx)
    
    % load network and simulation results
    isim = indx(ii); 
    if isempty(indx1) || any(indx1 == isim)
        disp(isim)
        simDir = dir(fullfile(d0, outFolder, ['R' num2str(isim) '_I*']));
        simCase = simDir.name;
        xOut = load(fullfile(d0, outFolder, simCase, simCase));
        
        % SIMULATION
        inputRaw = xOut.simModel.inputRaw;
        simModel = inputSetup(inputRaw);
        simOptions = xOut.simOptions;
        simOptions.isave = 0;
        simOptions.seed0 = 2*round(cputime*1e3)+1;
        simOptions.Nrun = Nrun;
        simOptions.Nrand = 2e5;
        simOptions.simPar = 1;
        simOptions.mex = 1;
        sysProp = sysProperties('TL', inputRaw.inputT, inputRaw.inputS, inputRaw.Nstate);
        indxR1 = sysProp.sccIndx{cellfun(@(x) any(x== 1), sysProp.sccIndx)}; % assuming that critical SCC is upstream
        
        % figure; plotNetwork(sysProp, ['#' num2str(isim)]); saveas(gcf, fullfile(simOutFolder, ['network_' num2str(isim)]))
        
        % update param
        [~, ~, ~, p0, q0, mu0] = findEquivalentParam(simModel, simOptions);
        disp([p0 q0 p0/(1-q0)])
        % update only critical SCC
        inputRaw1 = findModifiedRates(p_tg(ii), inputRaw, indxR1, Rmax, simOptions);
        if ~isempty(inputRaw1)
            simModel1 = inputSetup(inputRaw1);
            [~, ~, ~, p1, q1, mu1] = findEquivalentParam(simModel1, simOptions);
            % run simulation
            tic
            simOut = simCellStateLoop(simModel1, simOptions);
            toc
            out = ppro_simPaper(simModel1, simOptions, simOut, [], [5 10 15 20 25]);
            [lambda1Eq, lambda2Eq, gammaEq, p, q, mu] = findEquivalentParam(simModel1, simOptions);
            out.lambda1Eq = lambda1Eq; out.lambda2Eq = lambda2Eq; out.gammaEq = gammaEq;
            out.p = p; out.q = q; out.mu = mu;
            % out.xsol = xsol;
            xout{ii, 4} = out;
            xout{ii, 3} = simOptions;
            xout{ii, 2} = simModel1;
        else
            xout{ii, 4} = NaN;
            xout{ii, 3} = NaN;
            xout{ii, 2} = NaN;
        end
        xout{ii, 1} = isim;
    end
end

if isempty(indx1)
    % save results
    save(fullfile(simOutFolder, 'out_ppro'), 'xout', '-v7.3')
    
else
     % add/replace simulations
     A = load(fullfile(simOutFolder, 'out_ppro'), 'xout');
     A.xout(indx1,:) = xout(indx1,:);
     xout = A.xout;
     save(fullfile(simOutFolder, 'out_ppro'), 'xout', '-v7.3')
end
