clear all
close all

% This script generates the figure to show how the distribution in the GPA
% models converge to an exponential. 
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
isave = 1;
outDir = 'C:\Users\cp4u17\OneDrive - University of Southampton\Work\code\simCellState\io\OUT\GENERIC\';
xrange = [5 95]; % range of percentiles
x50 = 50; % average
xaxisDistr = [0 5 0.01 5]; iextr = 0;
extLev = [80 90 95 98]/100; dx = 0.25;

figDir = fullfile(pwd, 'fig_pproNetwork');
if ~exist(figDir, 'dir')
    mkdir(figDir)
end
outFolder = 'ppro_out_ncons_20190409150503';
simOut = 'out_ppro.mat';
xcol = lines(7); xcol(4,:) = [];

hfig = figure; setFigureProp(hfig, 10); hold on; grid on
xlabel('$x$'); ylabel('$P(x)$');
xleg = [];
for ii = 1:length(extLev)
    
    % load data and filter simulations that do not reach the required level of extinction
    A = load(fullfile(outDir, outFolder, simOut));
    iextRate = cellfun(@(x) ~isempty(find(x.extRate>extLev(ii), 1, 'first')), A.xout(:,4));
    A.xout(iextRate == 0,:) = [];
    ext = cellfun(@(x) x.extRate(find(x.extRate>extLev(ii), 1, 'first')), A.xout(:,4), 'UniformOutput', false);
    
    % compute distributions
    nn_sim = cellfun(@(x) x.nn_sim{find(x.extRate>extLev(ii), 1, 'first')}, A.xout(:,4), 'UniformOutput', false);
    xdist = cell(length(nn_sim),1);
    for jj = 1:length(nn_sim)
        [xx1, xxpdf] = getPdf(nn_sim{jj}, dx);
        xdist{jj} =  [xx1; xxpdf];
    end
    
    % compute envelope
    xint = 0.05:0.05:10; minPdf = 1e-10;
    [xint, yint, N, ~, indxOut] = getEnvelopeSim(xdist, xint, minPdf, [], xrange, iextr);
        
%     [cminMax, hminMax] = contour(xint, yint, N, xrange);
%     hminMax.LineColor = xcol(ii,:);
    [c50, h50] = contour(xint, yint, N, x50*[1 1]);
    h50.LineColor = xcol(ii,:);
    h50.LineWidth = 2;
    xleg = [xleg {sprintf('Extinction Fraction = %.d%%',extLev(ii)*100)}];
end
set(gca, 'Layer', 'top', 'YScale', 'log')
axis(xaxisDistr)
plot(linspace(xaxisDistr(1), xaxisDistr(2), 100), exppdf(linspace(xaxisDistr(1), xaxisDistr(2), 100)), 'k', 'LineWidth', 2)
legend([xleg 'Exp(1)'])

if isave
    saveas(gcf, fullfile(figDir, [outFolder, '_conv']))
end