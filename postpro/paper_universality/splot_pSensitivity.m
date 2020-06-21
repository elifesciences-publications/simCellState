clear all
close all

% This script is run to prepare the figures for the main text.
% Output comes from sprocess_genericSimPaper.m and srun_simGADModifiedRates.m
% Fig.4 (a)
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
icase = 0; isave = 1;
iSimGAD = 870;
outDir = '..\..\io\OUT\GENERIC\';
p_tg = [0.5 1 5 10 30];

% figure
hfig = figure; setFigureProp(hfig, 20); hold on; grid on
xlabel('$x$'); ylabel('$P(x)$');
figDir = fullfile(pwd, 'fig_pproNetwork');
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

% other parameters
dx = 0.2; xaxisDistr = [0 3 0.01 4]; xaxisMean = [0 1 0 25];
iTimeNorm = '20rmin';
% cmap = colormap('gray'); xcmap = linspace(35, -10, size(cmap, 1)); % to interpolate
cmap = [0 0 0; 1 0 0]; xcmap = linspace(log10(0.4), log10(35), size(cmap, 1)); % to interpolate

% loop on case
for isim = 1:length(iSimGAD)
    
    % load simulation
    outFolder = ['out_cons_20190409150503_indx' num2str(iSimGAD(isim))];
    A = load(fullfile(outDir, outFolder, 'out_ppro.mat'));
    
    % get mean and distributions
    switch iTimeNorm
        case '20rmin'
            iTimeFilt = 20;
            ifilt = cellfun(@(x) ~isempty(find(x.iTime>=iTimeFilt, 1, 'first')), A.xout(:,4));
            A.xout(ifilt == 0,:) = [];
            xp = cellfun(@(x) x.p, A.xout(:,4)); % xp(xp<1) = round(xp(xp<1)*10)/10; xp(xp>=1) = round(xp(xp>=1));
            % ifilt1 = ismember(xp, p_tg); 
            ifilt1 = zeros(size(xp));
            for ii = 1:length(xp)
                if any(abs(p_tg-xp(ii)) < 0.5)
                    ifilt1(ii) = 1;
                end
            end
            A.xout(ifilt1 == 0,:) = []; xp(ifilt1 == 0) = [];
            xmean = cellfun(@(x) [x.time_sim(1:x.indxT(find(x.iTime>=iTimeFilt, 1, 'first')))/x.time_sim(x.indxT(find(x.iTime>=iTimeFilt, 1, 'first'))); ...
                x.mean_n_sim(1:x.indxT(find(x.iTime>=iTimeFilt, 1, 'first')))], A.xout(:,4), 'UniformOutput', false); % scaled by 20/rmin
            nn_sim = cellfun(@(x) x.nn_sim{find(x.iTime>=iTimeFilt, 1, 'first')}, A.xout(:,4), 'UniformOutput', false);
    end
    xdist = cell(length(nn_sim),1); hprof = zeros(1,length(nn_sim)); xleg = cell(1,length(nn_sim));
    for jj = 1:length(nn_sim)
        nn_jj = nn_sim{jj};
        % nn_jj(nn_jj<mean(nn_jj)*0.25/2) = [];
        if xp(jj) <=1 
            [xx1, xxpdf] = getPdf(nn_jj, dx);
        elseif xp(jj) <=10
            [xx1, xxpdf] = getPdf(nn_jj, dx/2);
        else
            [xx1, xxpdf] = getPdf(nn_jj, dx/4);
        end
        xdist{jj} =  [xx1; xxpdf];
        % plot
        hprof(jj) = plot(xdist{jj}(1,:), xdist{jj}(2,:), 'color', interp1(xcmap, cmap, log10(xp(jj))), 'LineWidth', 1);
        xleg{jj} = sprintf('$%s_R = %.1f$', '\hat \lambda', p_tg(jj)); % patch here
    end

end
axis(xaxisDistr)
legend(hprof, xleg, 'Interpreter', 'LaTex')
% colormap(interp1(xcmap, cmap, log10(0.5:0.5:30)))
% hcb = colorbar; caxis([0.5 30])

if isave
    saveas(hfig, fullfile(figDir, 'p_sensitivity'))
end