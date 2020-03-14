clear all
% close all

% Sensitivity analysis to p parameter.
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
p_tg = [20];
Nrun = 500; % 5e2;

% indx = [591 122 756 777 33 57]; %1:1:100;
indx = [9]; %1:1:100;
simOutFolder = fullfile(pwd, 'out_sens', sprintf('%s_%d', outFolder, indx));
if ~exist(simOutFolder, 'dir')
    mkdir(simOutFolder)
end

% initalization
xout = cell(length(p_tg), 4);
% loop
ii = 1; % single simulation
for jj = 1:length(p_tg)
    
    % load network and simulation results
    isim = indx(ii); disp(p_tg(jj))
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
    simOptions.Nrand = 1e6;
    sysProp = sysProperties('TL', inputRaw.inputT, inputRaw.inputS, inputRaw.Nstate);
    indxR1 = sysProp.sccIndx{cellfun(@(x) any(x== 1), sysProp.sccIndx)}; % assuming that critical SCC is upstream
    [~, indxS] = intersect(inputRaw.inputS(:,2), indxR1);
    [~, indxT] = intersect(inputRaw.inputT(:,2), indxR1);
    
    % update param
    [~, ~, ~, p0, q0, mu0] = findEquivalentParam(simModel, simOptions);
    % fun = @(x) updateModel(x, inputRaw, indxR1, simOptions);
    fun = @(x) updateModel(x, inputRaw, indxT, indxS, simOptions);
    opt = MOPSOoptions('Nobj', 1, 'plot', 1);
    opt.stop = 15;
    [x,Ygbest] = MOPSO(@(x, iter) abs(fun(x) - p_tg(jj)), 0.1*ones(1, length(indxS)+length(indxT)), ...
        50*ones(1, length(indxS)+length(indxT)), opt);
    [~, ~, inputRaw1] = fun(x);
    simModel1 = inputSetup(inputRaw1);
    [~, ~, ~, p1, q1, mu1] = findEquivalentParam(simModel1, simOptions);
    % run simulation
    tic
    simOut = simCellStateLoop(simModel1, simOptions);
    toc
    xout{jj, 4} = ppro_GAD(simModel1, simOptions, simOut);
    xout{jj, 3} = simOptions;
    xout{jj, 2} = inputRaw1;
    xout{jj, 1} = isim;
end

% save results
save(fullfile(simOutFolder, 'out'), 'xout')

% PLOT
% ----
xp = []; hh = []; xleg = [];
figure; hold on; grid on
for ii = 1:size(xout,1)
    if isstruct(xout{ii, 4})
        % hii = plot(xout{ii, 4}.xx2, xout{ii, 4}.xx2pdf, '.--', 'MarkerSize', 15);
        hii = plot(xout{ii, 4}.xx2, xout{ii, 4}.xx2pdf, '.-', 'MarkerSize', 6);
        hh = [hh hii];
        xleg = [xleg {sprintf('p = %.2f, q = %.2f', xout{ii, 4}.p, xout{ii, 4}.q)}];
        xp = [xp xout{ii, 4}.p];
    end
end
legend(hh, xleg)
xx20 = 0:0.01:4;
% plot(xx20, normpdf(xx20, 1, sqrt(1/p_tg)), 'k', 'LineWidth', 2)
set(gca, 'XLim', [0 3])
saveas(gcf, fullfile(simOutFolder, 'distribution'))

figure; hold on; grid on
plotNetwork(sysProp, ['# ' num2str(indx)], 'EdgeLabel', [])
saveas(gcf, fullfile(simOutFolder, 'network'))