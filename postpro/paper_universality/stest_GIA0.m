clear all
close all

% Script to analyze the GAD equivalent model and its approximations.
% Fig.5, 6, 7, 8 in SI
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

addpath('../GenericNetwork')

if ~exist(fullfile(pwd, 'fig'), 'dir')
    mkdir(fullfile(pwd, 'fig'))
end

% domain
p = 0.01:0.01:100;
q = 0:0.01:0.9999;

% solution
% Pn = @(n, p, q) p/q*(1-q)^(p/q)/gamma(p/q + 1)*q.^(n).*gamma(p/q + n)./gamma(n+1);
% Pn1 = @(n, p, q) p/q*(1-q)^(p/q)*q.^(n).*exp(gammaln(p/q + n) - gammaln(p/q + 1))./gamma(n+1); % inf for large p
Pn1 = @(n, p, q) (1-q)^(p/q)*q.^(n).*exp(gammaln(p/q + n) - gammaln(p/q))./gamma(n+1); % inf for large p
% Px = @(x, p, q) (1-q)^(p/q)/gamma(p/q)./x.*q.^(x*p/(1-q)).*gamma(p/q + p/(1-q)*x)./gamma(p/(1-q)*x);
% Px1 = @(x, p, q) (1-q)^(p/q)./x.*q.^(x*p/(1-q)).*exp(gammaln(p/q + p/(1-q)*x)-gammaln(p/q)-gammaln(p/(1-q)*x));
Px1 = @(x, p, q) (1-q)^(p/q)./(q*(x - 1) + 1).*q.^(x*p/(1-q)+1).*exp(gammaln(p/q + p/(1-q)*x +1)-gammaln(p/q)-gammaln(p/(1-q)*x +1));
% simulation parameters
Nrun = 1e4; xTime = 0:1:20; G = 1;

% contour
[Q, P] = meshgrid(q, p);
MU = P./(1-Q);

xlev = 10.^(-2:0.5:4); % [0.01 0.05 0.1 0.5 1 5 10 50 100 500 1e3 1e4];
hfig = figure; setFigureProp(hfig, 10); hold on; grid on
[c0, h] = contourf(Q, P, log10(MU), log10(xlev), 'LineColor', 'w');
colormap(gray); % caxis(log10([0 1e4]))
% c = contourc(q, p, MU, xlev); clabel(c, h);
xlabel('$\hat \lambda_2$'); ylabel('$\hat \lambda_1$')
set(gca, 'YScale', 'log')
hcb=colorbar;
title(hcb,'$\log_{10}\bar{n}^*_2$', 'Interpreter', 'Latex')

xcol = lines(7);

% 1) q -> 0 Poisson
p1 = [0.5 1 5 10]; q1 = 0.01*ones(size(p1)); 
hp = plot(q1, p1, 'ok', 'linewidth', 0.5, 'MarkerFaceColor', xcol(3,:));
lab1 = sprintf('$%s_2=%.2g$', '\hat \lambda', q1(1));

% 2) q -> 1 Gamma
p2 = [0.5 1 5 10]; q2 = 0.99*ones(size(p2)); 
hg = plot(q2, p2, 'vk', 'linewidth', 0.5, 'MarkerFaceColor', xcol(2,:));
lab2 = sprintf('$%s_2=%.2g$', '\hat \lambda', q2(1));

% 3) p large (norm)
q3 = [0.1 0.4 0.6 0.9]; p3 = 60*ones(size(q3)); 
hn = plot(q3, p3, 'sk', 'linewidth', 0.5, 'MarkerFaceColor', xcol(1,:));
lab3 = sprintf('$%s_1=%.2g$', '\hat \lambda', p3(1));

hl = legend([hp, hg, hn], {lab1, lab2, lab3}, 'Interpreter', 'Latex');
set(hl, 'Location', 'best')
set(gca,'layer','top');
saveas(hfig, fullfile(pwd, 'fig', 'map_qVSp'))

% run simulations
% ---------------
% 1) Poisson
hh = []; xleg = [];
hfig1 = figure; setFigureProp(hfig1, 10); hold on; grid on
for ii = 1:length(p1)
   simOut = run_simGAD(p1(ii), q1(ii), G, xTime, Nrun);
   nn_sim = simOut.nn_sim{1}-1;
   [nnpdf, nnEdges] = histcounts(nn_sim, 'BinLimits', [0.5 max(nn_sim)+0.5], 'BinWidth', 1,'Normalization', 'pdf');
   nn = nnEdges(1:end-1) + diff(nnEdges)/2;
   ha = plot(nn, Pn1(nn, p1(ii), q1(ii)), '-', 'linewidth', 0.5, 'Color', 'k'); % xcol(ii,:));
   hp = plot(nn, poisspdf(nn, p1(ii)), '--', 'linewidth', 2, 'Color', 'k'); % , xcol(ii,:));
   hs = plot(nn, nnpdf, '.', 'MarkerSize', 15, 'Color', xcol(ii,:));
   hh = [hh, hs]; 
   xleg = [xleg, {sprintf('Simulation, $%s_1 = %.2g$', '\hat \lambda', p1(ii))}];
end
xlabel('$n_2$'); ylabel('$P^*(n_2)$'); 
title(sprintf('%s', lab1), 'Interpreter', 'Latex')
legend([ha hp hh], ['Analytic solution', 'Approximation', xleg], 'Interpreter', 'Latex')
saveas(hfig1, fullfile(pwd, 'fig', 'seq_pois'))

% 2) Gamma
hh = []; xleg = []; dx = 0.05;
hfig2 = figure; setFigureProp(hfig2, 10); hold on; grid on
for ii = 1:length(p2)
   simOut = run_simGAD(p2(ii), q2(ii), G, xTime, Nrun);
   nn_sim = simOut.nn_sim{1}-1;
   xx0 = 1/mean(nn_sim)/2:0.01:max(nn_sim)/mean(nn_sim)+1/mean(nn_sim)/2; 
   [xx1, xxpdf] = getPdf(nn_sim, dx);
   ha = plot(xx0, Px1(xx0, p2(ii), q2(ii)), '-', 'linewidth', 0.5, 'Color', 'k'); % , xcol(ii,:));
   hp = plot(xx0, gampdf(xx0, p2(ii), 1/p2(ii)), '--', 'linewidth', 2, 'Color', 'k'); % , xcol(ii,:));
   hs = plot(xx1, xxpdf, '.', 'MarkerSize', 15, 'Color', xcol(ii,:));
   hh = [hh, hs]; 
   xleg = [xleg, {sprintf('Simulation, $%s_1 = %.2g$', '\hat \lambda', p2(ii))}];
end
xlabel('$x_2$'); ylabel('$P^*(x_2)$'); 
title(sprintf('%s', lab2), 'Interpreter', 'Latex')
legend([ha hp hh], ['Analytic solution', 'Approximation', xleg], 'Interpreter', 'Latex')
set(gca, 'XLim', [0 2.5])
saveas(hfig2, fullfile(pwd, 'fig', 'seq_gamma'))

% 3) Normal
hh = []; xleg = []; xx0 = 0:0.01:4; dx = 0.05;
hfig3 = figure; setFigureProp(hfig3, 10); hold on; grid on
for ii = 1:length(p3)
   simOut = run_simGAD(p3(ii), q3(ii), G, xTime, Nrun);
   nn_sim = simOut.nn_sim{1}-1;
   [xx1, xxpdf] = getPdf(nn_sim, dx);
   ha = plot(xx0, Px1(xx0, p3(ii), q3(ii)), '-', 'linewidth', 0.5, 'Color', 'k'); % , xcol(ii,:));
   hp = plot(xx0, normpdf(xx0, 1, 1/sqrt(p3(ii))), '--', 'linewidth', 2, 'Color', 'k');
   hs = plot(xx1, xxpdf, '.', 'MarkerSize', 15, 'Color', xcol(ii,:));
   hh = [hh, hs]; 
   xleg = [xleg, {sprintf('Simulation, $%s_2 = %.2g$', '\hat \lambda', q3(ii))}];
end
xlabel('$x_2$'); ylabel('$P^*(x_2)$'); 
title(sprintf('%s', lab3), 'Interpreter', 'Latex')
legend([ha hp hh], ['Analytic solution', 'Approximation', xleg], 'Interpreter', 'Latex')
set(gca, 'XLim', [0 2.5])
saveas(hfig3, fullfile(pwd, 'fig', 'seq_norm'))

rmpath('../GenericNetwork')
