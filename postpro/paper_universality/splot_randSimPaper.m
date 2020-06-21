clear all
close all

% This script is run to prepare the figures for the main text and SI
% Output comes from sprocess_genericSimPaper.m and
% srun_simGADModifiedRates.m and it is currently saved under io folder.
% 
% Fig.2/3 (a): set icase = 0
% Fig.2/3 (b): set icase = 2
% Fig.4 (b): set icase = 1
% Fig.14 (SI): set icase = 10
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

% setting parameters
icase = 0; 
iplot = []; isave = 0; xx0Norm = [];
iSimGAD = [94 57 30 302 870];
iSimGPA = [627 733 741 768 897];
outDir = '..\..\io\OUT\GENERIC\';
inormSigma = 0; tgOut = '';
switch icase
    case 0 % GAD: initial networks
        outFolder = 'ppro_out_cons_20190409150503';
        A = load(fullfile(outDir, outFolder, 'out_ppro.mat'));
        dx = 0.2; xaxisDistr = [0 3 0.01 5]; xaxisMean = [0 1 0 25];
        A.xout(:,1) = mat2cell(1:1:size(A.xout,1), 1, ones(size(A.xout,1),1))'; % missing sim newtork
        % filter out those with mean below 2
        mm = cellfun(@(x) x.mean_n_sim(end)-1, A.xout(:,4));
        A.xout(mm<2,:) = [];
        iplot = 'dist_exp'; iextr = 1;
        iTimeNorm = '20rmin';
        indxSim0 = iSimGAD;
    case 1 % GAD: p=30
        outFolder = 'out_cons_20190409150503_p30';
        A = load(fullfile(outDir, outFolder, 'out_ppro.mat'));
        xaxisDistr = [-4 4 0 0.6]; xaxisMean = [0 1 0 600];
        iplot = 'dist_normSigma'; iextr = 1;
        dx = 0.05; inormSigma = 1;
        istr = cellfun(@(x) isstruct(x), A.xout(:,2));
        A.xout(istr == 0,:) = [];
        iTimeNorm = '20rmin';
        indxSim0 = iSimGAD;
    case 10 % GAD: p=30, but not scaled by variance
        outFolder = 'out_cons_20190409150503_p30';
        A = load(fullfile(outDir, outFolder, 'out_ppro.mat'));
        xaxisDistr = [0 2.5 0.01 4]; xaxisMean = [0 1 0 600];
        iplot = 'dist_norm'; iextr = 1;
        dx = 0.05; inormSigma = 0;
        istr = cellfun(@(x) isstruct(x), A.xout(:,2));
        A.xout(istr == 0,:) = [];
        iTimeNorm = '20rmin';
        indxSim0 = iSimGAD; tgOut = '_ns';
    case 2 % GPA: initial networks
        outFolder = 'ppro_out_ncons_20190409150503';
        A = load(fullfile(outDir, outFolder, 'out_ppro.mat'));
        dx = 0.25; xx0Norm = [];
        xaxisDistr = [0 3 0.01 5]; xaxisMean = [0 1 0 150];
        iplot = 'dist_exp'; iextr = 0;
        iTimeNorm = '98ext';
        indxSim0 = iSimGPA;
end
% get mean and distributions
switch iTimeNorm
    case 'final'
        xmean = cellfun(@(x) [x.time_sim(1:x.indxT(end))/x.time_sim(x.indxT(end)); x.mean_n_sim(1:x.indxT(end))], A.xout(:,4), 'UniformOutput', false);
        nn_sim = cellfun(@(x) x.nn_sim{find(x.extRate>=extLev, 1, 'first')}, A.xout(:,4), 'UniformOutput', false); % ?
    case '20rmin'
        iTimeFilt = 20;
        ifilt = cellfun(@(x) ~isempty(find(x.iTime>=iTimeFilt, 1, 'first')), A.xout(:,4));
        A.xout(ifilt == 0,:) = [];
        xmean = cellfun(@(x) [x.time_sim(1:x.indxT(find(x.iTime>=iTimeFilt, 1, 'first')))/x.time_sim(x.indxT(find(x.iTime>=iTimeFilt, 1, 'first'))); ...
            x.mean_n_sim(1:x.indxT(find(x.iTime>=iTimeFilt, 1, 'first')))], A.xout(:,4), 'UniformOutput', false); % scaled by 20/rmin
        nn_sim = cellfun(@(x) x.nn_sim{find(x.iTime>=iTimeFilt, 1, 'first')}, A.xout(:,4), 'UniformOutput', false);
    case '98ext'
        extLev = 98/100;
        iextRate = cellfun(@(x) ~isempty(find(x.extRate>=extLev, 1, 'first')), A.xout(:,4));
        A.xout(iextRate == 0,:) = [];
        xmean = cellfun(@(x) [x.time_sim(1:x.indxT(find(x.extRate>=extLev, 1, 'first')))/x.time_sim(x.indxT(find(x.extRate>=extLev, 1, 'first'))); ...
            x.mean_n_sim(1:x.indxT(find(x.extRate>=extLev, 1, 'first')))], A.xout(:,4), 'UniformOutput', false); % scaled by 20/rmin
        nn_sim = cellfun(@(x) x.nn_sim{find(x.extRate>=extLev, 1, 'first')}, A.xout(:,4), 'UniformOutput', false);
end
xdist = cell(length(nn_sim),1);
for jj = 1:length(nn_sim)
    nn_jj = nn_sim{jj};
    % nn_jj(nn_jj<mean(nn_jj)*0.25/2) = [];
    [xx1, xxpdf] = getPdf(nn_jj, dx);
    if ~isempty(xx0Norm)
        % renormalize pdf
        indx = find(xx1>xx0Norm);
        xxpdf = xxpdf/trapz(xx1(indx), xxpdf(indx));
    end
    if inormSigma
        xmu = mean(nn_jj); xsigma = std(nn_jj);
        xx1 = (xx1-1)*xmu/xsigma;
        xxpdf = xxpdf*xsigma/xmu;
    end
    xdist{jj} =  [xx1; xxpdf];
end
% xdist = cellfun(@(x) [x.xx; x.xxpdf], A.xout(:,4), 'UniformOutput', false);
[~, indxSim] = intersect([A.xout{:,1}], indxSim0); [A.xout{indxSim,1}];

% interpolation and counting
xfilt = [5 95]; % range of percentiles

figDir = fullfile(pwd, 'fig_pproNetwork');
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

% distribution
if inormSigma
    xint = -4:0.01:4; minPdf = 1e-10;
else
    xint = 0:0.01:5; minPdf = 1e-10;
end
[xint, yint, N, ~, indxOut] = getEnvelopeSim(xdist, xint, minPdf, [], xfilt, iextr);
hdist = plotEnvelopeSim(xint, yint, N, xfilt, xaxisDistr, iplot, xdist(indxSim));
if isave
    saveas(hdist, fullfile(figDir, [outFolder, '_dist' tgOut]))
end

% mean
xint = 0:0.01:1;
[xint, yint, N] = getEnvelopeSim(xmean, xint, [], xaxisMean(end)+10, xfilt, iextr);
hmean = plotEnvelopeSim(xint, yint, N, xfilt, xaxisMean, 'mean', xmean(indxSim));
if isave
    saveas(hmean, fullfile(figDir, [outFolder, '_mean' tgOut]))
end