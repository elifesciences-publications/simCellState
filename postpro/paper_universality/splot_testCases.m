clear all
close all

% This script generates the figures for the paper (SI)
% - Clonal dynamics: mean
% - Clonal dynamics: distribution (sim VS master Equation)
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
outDir = '../../';
Tref = 10; % time for plotting the distribution
nAMaxAD = 3; nBMaxAD = 30;
nAMaxPA = 150; nBMaxPA = 160;
figPath = fullfile(pwd, 'fig');
if ~exist(figPath, 'dir')
    mkdir(figPath)
end
itest = [0 0 1]; % AD, PA, MS

% simulation path
simPathAD = {'AsymmetricDivision_20190712130850', 'AsymmetricDivision_20190712130905', 'AsymmetricDivision_20190712130920'}; NsimAD = length(simPathAD);
xlegAD = {'IA#1', 'IA#2', 'IA#3'};
simPathPA = {'PopulationAsymmetryR_20190712133058' 'PopulationAsymmetryR_20190712133259' 'PopulationAsymmetryR_20190712133458'}; NsimPA = length(simPathPA);
xlegPA = {'PA#1', 'PA#2', 'PA#3'};
simPathMS = {'PopulationAsymmetryM_20190715100003' 'PopulationAsymmetryM_20190715100334'}; NsimMS = length(simPathMS);
xlegMS = {'MS#1', 'MS#3'};
xcol = lines(10);

% initialization
lambdaAD = zeros(1, NsimAD); gammaAD = zeros(1, NsimAD);
simOutAD = cell(1, NsimPA); meanODEP = cell(1, NsimPA); QgfAD = cell(1, NsimPA);
lambdaPA = zeros(1, NsimPA); omegaPA = zeros(1, NsimPA); gammaPA = zeros(1, NsimPA);
simOutPA = cell(1, NsimPA); meanODEPA = cell(1, NsimPA); QgfPA = cell(1, NsimPA);
lambdaMS = zeros(1, NsimMS); gammaMS = zeros(1, NsimMS);
simOutMS = cell(1, NsimMS); meanODEMS = cell(1, NsimMS); QgfMS = cell(1, NsimMS);

if itest(1)
    itg = 'test_AD';
    % load simulation outputs
    for ii = 1:NsimAD
        % load outputs
        load(fullfile(outDir, 'io', 'OUT', 'AsymmetricDivision', simPathAD{ii}, simPathAD{ii}))
        lambdaAD(ii) = simModel.inputRaw.inputS(1,1);
        gammaAD(ii) = simModel.inputRaw.inputT(1,1);
        % postpro
        [simOutAD(ii), meanODEAD{ii}, QgfAD{ii}] = ppro_simCellState(simModel, simOptions, simOut);
        simOutAD{ii}.time = simOptions.time; % add time
        simOutAD{ii}.time_s = simOptions.time/Tref; % add time scaled
        % check if at final time I have NaN
        iok = cellfun(@(x) length(x), simOutAD{ii}.ppro.outStat_surv.indxNotNan);
        if any(iok<iok(1))
            disp('check NaN')
        end
    end
else
    NsimAD = 0; xlegAD = [];
end

if itest(2) || itest(3)
    itg = 'test_PA';
    for ii = 1:NsimPA
        % load outputs
        load(fullfile(outDir, 'io', 'OUT', 'PopulationAsymmetryR', simPathPA{ii}, simPathPA{ii}))
        lambdaPA(ii) = sum(simModel.inputRaw.inputS(:,1));
        rPA(ii) = simModel.inputRaw.inputS(1,1)/lambdaPA(ii);
        gammaPA(ii) = simModel.inputRaw.inputT(1,1);
        % postpro
        [simOutPA(ii), meanODEPA{ii}, QgfPA{ii}] = ppro_simCellState(simModel, simOptions, simOut);
        simOutPA{ii}.time = simOptions.time; % add time
        simOutPA{ii}.time_s = simOptions.time/Tref; % add time scaled
        [~, indxT] = min(abs(simOutPA{ii}.time-Tref));
        simOutPA{ii}.indxT = indxT;
        simOutPA{ii}.nn = sum(squeeze(simOutPA{ii}.xstate(simModel.indxActiveState, :, indxT)),1);
        % check if at final time I have NaN
        iok = cellfun(@(x) length(x), simOutPA{ii}.ppro.outStat_surv.indxNotNan);
        if any(iok<iok(1))
            disp('check NaN')
        end
    end
else
    NsimPA = 0; xlegPA = [];
end

if itest(3)
    itg = 'test_MS';
    for ii = 1:NsimMS
        % load outputs
        load(fullfile(outDir, 'io', 'OUT', 'PopulationAsymmetryM', simPathMS{ii}, simPathMS{ii}))
        lambdaMS(ii) = sum(simModel.inputRaw.inputS(:,1));
        rMS(ii) = simModel.inputRaw.inputS(1,1)/lambdaPA(ii);
        gammaMS(ii) = simModel.inputRaw.inputT(1,1);
        % postpro
        [simOutMS(ii), meanODEMS{ii}, QgfMS{ii}] = ppro_simCellState(simModel, simOptions, simOut);
        simOutMS{ii}.time = simOptions.time; % add time
        simOutMS{ii}.time_s = simOptions.time/Tref; % add time scaled
        [~, indxT] = min(abs(simOutMS{ii}.time-Tref));
        simOutMS{ii}.indxT = indxT;
        simOutMS{ii}.nn = sum(squeeze(simOutMS{ii}.xstate(simModel.indxActiveState, :, indxT)),1);
        % check if at final time I have NaN
        iok = cellfun(@(x) length(x), simOutMS{ii}.ppro.outStat_surv.indxNotNan);
        if any(iok<iok(1))
            disp('check NaN')
        end
    end
    NsimPA = 2; xlegPA(2) = [];
    lambdaPA(2) = []; omegaPA(2) = []; gammaPA(2) = [];
    simOutPA(2) = []; meanODEPA(2) = []; QgfPA(2) = [];
else
    NsimMS = 0; xlegMS = [];
end

% MEAN VALUE
% ----------
if itest(3)
    figure; setFigureProp(gcf, 10); subplot(2,1,1); hold on; grid on
else
    figure; setFigureProp(gcf, 10); hold on; grid on
end
for ii = 1:NsimAD
    plot(simOutAD{ii}.time_s, simOutAD{ii}.ppro.outStat_surv.meanT, '.--', 'color', xcol(ii,:), 'MarkerFaceColor', xcol(ii,:), 'MarkerSize', 10, 'LineWidth', 0.5)
end
for ii = 1:NsimPA
    plot(simOutPA{ii}.time_s, simOutPA{ii}.ppro.outStat_surv.meanT, '.--', 'color', xcol(ii+NsimAD,:), 'MarkerFaceColor', xcol(ii+NsimAD,:), 'MarkerSize', 10, 'LineWidth', 0.5)
end
for ii = 1:NsimMS
    plot(simOutMS{ii}.time_s, simOutMS{ii}.ppro.outStat_surv.meanT, '.--', 'color', xcol(ii+NsimAD+NsimPA,:), 'MarkerFaceColor', xcol(ii+NsimAD+NsimPA,:), 'MarkerSize', 10, 'LineWidth', 0.5)
end
legend([xlegAD xlegPA xlegMS], 'Location', 'northwest')
xlabel('t/$\tau$'); ylabel('$\bar{n}$')
set(gca, 'Xlim', [0 1]) % , 'Ylim', [0 50])
if ~itest(3)
    saveas(gcf, fullfile(figPath, [itg '_clone_mean']))
end

% EXTINCTION
% ----------
if itest(3)
    subplot(2,1,2); hold on; grid on
else
    figure; setFigureProp(gcf, 10); hold on; grid on
end
for ii = 1:NsimAD
    plot(simOutAD{ii}.time_s, simOutAD{ii}.ppro.Qsim, '.-', 'color', xcol(ii,:), 'MarkerFaceColor', xcol(ii,:), 'MarkerSize', 10, 'LineWidth', 0.5)
end
for ii = 1:NsimPA
    plot(simOutPA{ii}.time_s, simOutPA{ii}.ppro.Qsim, '.-', 'color', xcol(ii+NsimAD,:), 'MarkerFaceColor', xcol(ii+NsimAD,:), 'MarkerSize', 10, 'LineWidth', 0.5)
end
for ii = 1:NsimMS
    plot(simOutMS{ii}.time_s, simOutMS{ii}.ppro.Qsim, '.-', 'color', xcol(ii+NsimAD+NsimPA,:), 'MarkerFaceColor', xcol(ii+NsimAD+NsimPA,:), 'MarkerSize', 10, 'LineWidth', 0.5)
end
legend([xlegAD xlegPA xlegMS], 'Location', 'northwest')
xlabel('t/$\tau$'); ylabel('$P(n=0)$')
set(gca, 'Xlim', [0 1])
if ~itest(3)
    saveas(gcf, fullfile(figPath, [itg '_clone_extinction']))
else
    saveas(gcf, fullfile(figPath, [itg '_clone_MeanExt']))
end

% DISTRIBUTION
% ------------
if itest(1)
    % Asymmetric division
    figure; setFigureProp(gcf, 10); hold on; grid on
    hh = zeros(1,NsimAD);
    n = 0:1:20;
    for ii = 1:NsimAD
        [~, indxT] = min(abs(simOutAD{ii}.time-Tref));
        % mEqFnc = @(t, p, Na, Nb) mEq_asymmetric(t, p, Na, Nb, lambdaAD(ii), gammaAD(ii));
        % [n, Pn] = integrMEq(nAMaxAD, nBMaxAD, mEqFnc, simOutAD{ii}.time, indxT);
        hk = plot(n+1, poisspdf(n, lambdaAD(ii)/gammaAD(ii)), '.-k', 'MarkerSize', 15, 'LineWidth', 0.5);
        hh(ii) = plot(simOutAD{ii}.ppro.outStat_surv.distrT(indxT).xx, simOutAD{ii}.ppro.outStat_surv.distrT(indxT).pdf, ...
            '.-', 'color', xcol(ii,:), 'MarkerFaceColor', xcol(ii,:), 'MarkerSize', 15, 'LineWidth', 0.5);
    end
    xlabel('$n$'); ylabel('$P(n)$');
    legend([hh hk], [xlegAD {'Poisson(\lambda/\gamma)'}])
    set(gca, 'Xlim', [0 20])
    saveas(gcf, fullfile(figPath, 'clone_distributionAD'))
end

% Population asymmetry
if itest(2)
    figure; setFigureProp(gcf, 10); hold on; grid on
    hh = zeros(1,NsimPA);
    for ii = 1:NsimPA
        indxT = simOutPA{ii}.indxT;
        mEqFnc = @(t, p, Na, Nb) mEq_populationR(t, p, Na, Nb, lambdaPA(ii), rPA(ii), gammaPA(ii));
        [n, Pn] = integrMEq(nAMaxPA, nBMaxPA, mEqFnc, simOutPA{ii}.time, indxT);
        hk = plot(n, Pn, '.-k', 'MarkerSize', 15, 'LineWidth', 0.5);
        nn = simOutPA{ii}.nn; nn(nn==0) = [];
        [hc, edges] = histcounts(nn, 'BinLimits', [0.5 max(nn)+0.5], 'BinWidth', 1, 'Normalization', 'pdf');
        hh(ii) = plot(edges(1:end-1)+diff(edges)/2, hc, ...
            '.-', 'color', xcol(ii+NsimAD,:), 'MarkerFaceColor', xcol(ii+NsimAD,:), 'MarkerSize', 15, 'LineWidth', 0.5);
        % reference solution: plot([0 1 2:20], [1+exp(-Tref)-4/(Tref+2) -exp(-Tref)+8/(Tref+2)^2 8/(Tref+2)^2*(Tref/(Tref+2)).^((2:20)-1)])
        ss = 1:1:100; Pref = [-exp(-Tref)+8/(Tref+2)^2 8/(Tref+2)^2*(Tref/(Tref+2)).^(ss(2:end)-1)]; Pref = Pref/sum(Pref); % rescaled
        if ii == 1
            href = plot(ss, Pref, 'ok', 'MarkerSize', 7, 'LineWidth', 0.5);
        end
    end
    set(gca, 'YScale', 'log');
    set(gca, 'Xlim', [0 25])
    xlabel('$n$'); ylabel('$P(n)$');
    legend([hh hk href], [xlegPA {'Master Equation', 'Analytic Solution PA#1'}])
    saveas(gcf, fullfile(figPath, 'clone_distributionPA'))
end

if itest(3)
    figure; setFigureProp(gcf, 10); hold on; grid on
    hh = zeros(1,NsimPA+NsimMS);
    for ii = 1:NsimPA
        indxT = simOutPA{ii}.indxT;
        nn = simOutPA{ii}.nn; nn(nn==0) = [];
        [hc, edges] = histcounts(nn, 'BinLimits', [0.5 max(nn)+0.5], 'BinWidth', 1, 'Normalization', 'pdf');
        hh(ii) = plot(edges(1:end-1)+diff(edges)/2, hc, ...
            '.-', 'color', xcol(ii+NsimAD,:), 'MarkerFaceColor', xcol(ii+NsimAD,:), 'MarkerSize', 15, 'LineWidth', 0.5);
    end
    for ii = 1:NsimMS
        indxT = simOutMS{ii}.indxT;
        nn = simOutMS{ii}.nn; nn(nn==0) = [];
        [hc, edges] = histcounts(nn, 'BinLimits', [0.5 max(nn)+0.5], 'BinWidth', 1, 'Normalization', 'pdf');
        hh(ii+NsimPA) = plot(edges(1:end-1)+diff(edges)/2, hc, ...
            '.-', 'color', xcol(ii+NsimAD+NsimPA,:), 'MarkerFaceColor', xcol(ii+NsimAD+NsimPA,:), 'MarkerSize', 15, 'LineWidth', 0.5);
    end
    set(gca, 'YScale', 'log');
    set(gca, 'Xlim', [0 25])
    xlabel('$n$'); ylabel('$P(n)$');legend(hh, [xlegPA xlegMS])
    saveas(gcf, fullfile(figPath, 'clone_distributionMS'))
end