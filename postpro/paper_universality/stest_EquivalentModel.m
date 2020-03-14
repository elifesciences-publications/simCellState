clear all
close all

% Script to check the applicability of the equivalent model as an
% approximation of a generic complex model.
% Figures of section Equivalent model for complex GAD networks.
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
d0 = fullfile('..', '..', 'IO', 'OUT', 'GENERIC');
tg = '20190409150503';
outFolder = ['ppro_out_cons_' tg];
xOut = load(fullfile(d0, outFolder, 'out_ppro'));
Pn = @(n, p, q) (1-q)^(p/q)*q.^(n).*exp(gammaln(p/q + n) - gammaln(p/q))./gamma(n+1); % inf for large p
Px = @(x, p, q) (1-q)^(p/q)./(q*(x - 1) + 1).*q.^(x*p/(1-q)+1).*exp(gammaln(p/q + p/(1-q)*x +1)-gammaln(p/q)-gammaln(p/(1-q)*x +1));
xcol = lines(7);

tgPoi = 0.2; tgNorm = 20; tgGamma = 0.8; tgMu = 10;

outDir = fullfile(pwd, 'figGAD');
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

indx = 1:1:size(xOut.xout,1);
xp = zeros(1, length(indx)); xq = zeros(1, length(indx)); xtype = zeros(1, length(indx));
xdata = cell(length(indx), 9);
for ii = 1:length(indx)
    
    % load and postprocess simulation
    outii = xOut.xout{ii, 4};
    isim = xOut.xout{ii,1};
    p = outii.p; q = outii.q;
    q(q<1e-10) = 0;
    xp(ii) = p; xq(ii) = q; mu = p/(1-q);
    nn_simjj = outii.nn_sim{1}-1; % shift by one (GAD model, 1 renewing cell)
    simtg = sprintf('$%s = %.2f, %s = %.2f, %s{n}^*_C = %.2f$', '\hat\lambda_1', p, '\hat\lambda_2', q, '\bar', mu);
    % Plot Network
    % figure; hold on; grid on
    % plotNetwork(sysProp, '');
    % saveas(gcf, fullfile(pwd, 'figGAD', ['network_' num2str(isim)]))
    
    % PLOT
    %     % mean
    %     figure; hold on; grid on
    %     hEq = plot(out.xsol.timeEq, out.xsol.n2Eq);
    %     hSim = plot(out.time_sim, out.mean_n2_sim, '.--');
    %     xlabel('time [t]'); ylabel('n_2^*')
    %     set(gca, 'XLim', [0 max(out.time_sim)])
    %     legend([hEq, hSim], {'Eq. model',  'GIA random model'})
    %     saveas(gcf, fullfile(outDir, ['mean_' num2str(isim)]), 'fig')
    
    % distribution
    dn = max([round(mean(nn_simjj)/5) 1]); dx = 0.1;
    xx20 = 0:0.01:4;
    nn20 = 0:1:max(nn_simjj);
    if length(indx) == 2
        figure(2); subplot(1,2,ii); hold on; grid on
    else
        hfig0 = figure; setFigureProp(hfig0, 10); hold on; grid on
    end
    hSim = []; ntg = 'oth_';
    if mu > tgMu % cont
        [xx, xxpdf] = getPdf(nn_simjj, dx); % continuous distribution
        hEqA = plot(xx20, Px(xx20, p, q), '-', 'linewidth', 0.5, 'Color', 'k');
        hSim = plot(xx, xxpdf, '.', 'MarkerSize', 15, 'Color', xcol(2,:));
        if p > tgNorm % normal
            xtype(ii) = 1;
            hEq = plot(xx20, normpdf(xx20, 1, sqrt(1/p)), '--', 'linewidth', 2, 'Color', 'k');
            % hSim = plot(xx, xxpdf, '.', 'MarkerSize', 15, 'Color', xcol(4,:));
            set(gca, 'XLim', [0 2.5])
            hh = [hEqA, hEq, hSim];
            % xleg = {'Analytic solution', sprintf('Eq. model, Approximation (%s)', 'Norm(1, $1/\tilde p$)'), 'GIA random model'};
            xleg = {'Analytic solution', sprintf('Eq. model, Approximation'), 'GIA random model'};
            ntg = 'norm_';
        elseif q >= tgGamma % gamma
            xtype(ii) = 3;
            hEq = plot(xx20, gampdf(xx20, p, 1/p), '--', 'linewidth', 2, 'Color', 'k');
            % hSim = plot(xx, xxpdf, '.', 'MarkerSize', 15, 'Color', xcol(1,:));
            set(gca, 'XLim', [0 4])
            ntg = 'gamma_';
            hh = [hEqA, hEq, hSim];
            % xleg = {'Eq. model, Analytic solution', sprintf('Eq. model, Approximation (%s)', '$\Gamma(\tilde p, 1/\tildep)$'), 'GIA random model'};
            xleg = {'Eq. model, Analytic solution', sprintf('Eq. model, Approximation'), 'GIA random model'};
        else % other
            % hSim = plot(xx, xxpdf, '.', 'MarkerSize', 15, 'Color', xcol(5,:));
            hh = [hEqA, hSim];
            xleg = {'Eq. model, Analytic solution', 'GIA random model'};
            set(gca, 'XLim', [0 4])
        end
        xlabel('$x_C$'); ylabel('$P^*(x_C)$');
        mse = immse(xxpdf,Px(xx, p, q));
        errMax = max(abs(xxpdf-Px(xx, p, q)));
        maxpdfref = max(Px(xx, p, q));
        maxpdf = max(xxpdf);
    else
        [nn, nnpmf] = getPmf(nn_simjj, dn);
        if q == 0
            hEqA = plot(nn20, poisspdf(nn20, p), '-', 'linewidth', 0.5, 'Color', 'k');
        else
            hEqA = plot(nn20, Pn(nn20, p, q), '-', 'linewidth', 0.5, 'Color', 'k');
        end
        hSim = plot(nn, nnpmf, '.', 'MarkerSize', 15, 'Color', xcol(2,:));
        if q <= tgPoi
            xtype(ii) = 2;
            hEq = plot(nn20, poisspdf(nn20, p), '--', 'linewidth', 2, 'Color', 'k');
            % hSim = plot(nn, nnpmf, '.', 'MarkerSize', 15, 'Color', xcol(2,:));
            ntg = 'pois_';
            hh = [hEqA, hEq, hSim];
            % xleg = {'Eq. model, Analytic solution', sprintf('Eq. model, Approximation (%s)', 'Poisson($\tilde p$)'), 'GIA random model'};
            xleg = {'Eq. model, Analytic solution', sprintf('Eq. model, Approximation'), 'GIA random model'};
        else
            % hSim = plot(nn, nnpmf, '.', 'MarkerSize', 15, 'Color', xcol(5,:));
            hh = [hEqA, hSim];
            xleg = {'Eq. model, Analytic solution', 'GIA random model'};
        end
        xlabel('$n_C$'); ylabel('$P^*(n_C)$'); 
        if q == 0
            mse = immse(nnpmf, poisspdf(nn, p));
            errMax = max(abs(nnpmf-poisspdf(nn, p)));
            maxpdfref = max(poisspdf(nn, p));
        else
            mse = immse(nnpmf, Pn(nn, p, q));
            errMax = max(abs(nnpmf-Pn(nn, p, q)));
            maxpdfref = max(Pn(nn, p, q));
        end
        maxpdf = max(nnpmf);
    end
    if ~isempty(hSim)
        title(sprintf('%s', simtg), 'Interpreter', 'latex')
        legend(hh, xleg, 'Interpreter', 'latex')
        saveas(hfig0, fullfile(outDir, [ntg 'distr_' num2str(isim)]), 'fig')
    end
    % saveas(gcf, fullfile(pwd, 'figGAD', ['distr_' num2str(isim)]), 'jpg')
    % close when running multiple cases
    if length(indx) > 2
        close(hfig0)
    end
    
    xdata(ii,:) = {isim p q mu xtype(ii) mse errMax maxpdf maxpdfref};
end

% selected points
indxs = [591 122 ... % Poisson
    33 57 ... % gamma 19 
    19 777]; % other

% plot summary
hfig = figure; setFigureProp(hfig, 10); hold on; grid on
plot([xdata{:,3}], [xdata{:,7}]./[xdata{:,8}], '.', 'MarkerSize', 10, 'color', xcol(1,:))
plot([xdata{indxs,3}], [xdata{indxs,7}]./[xdata{indxs,8}], 'o', 'LineWidth', 2)
legend('GIA random model', 'Selected model', 'Location', 'best')
xlabel('$\hat \lambda_2$'); ylabel('$\epsilon$')
if length(indx) > 2
    saveas(hfig, fullfile(outDir, 'err_qVSp_sim'))
end

% domain
p = 0.01:0.01:100;
q = 0:0.01:0.9999;

% contour
[Q, P] = meshgrid(q, p);
MU = P./(1-Q);

xlev = 10.^(-2:0.5:4);

hfig1 = figure; setFigureProp(hfig1, 10); hold on; grid on
[c0, h] = contourf(Q, P, log10(MU), log10(xlev), 'LineColor', 'w');
colormap(gray); % caxis(log10([0 1e4]))
% c = contourc(q, p, MU, xlev); clabel(c, h);
xlabel('$\hat \lambda_2$'); ylabel('$\hat \lambda_1$')
set(gca, 'YScale', 'log')
hcb=colorbar;
title(hcb,'$\log_{10}\bar{n}^*_C$', 'Interpreter', 'Latex')

% plot simulation points
ha = plot(xq, xp, '.', 'MarkerSize', 10, 'color', xcol(1,:));
% ha = plot(xq(xtype==0), xp(xtype==0), '.', 'color', xcol(5,:));
% hp = plot(xq(xtype==2), xp(xtype==2), '+', 'linewidth', 2, 'color', xcol(2,:));
% hg = plot(xq(xtype==3), xp(xtype==3), 'x', 'linewidth', 2, 'color', xcol(1,:));
% hn = plot(xq(xtype==1), xp(xtype==1), '*', 'linewidth', 1, 'color', xcol(4,:));

% selected points
hs = plot(xq(indxs), xp(indxs), 'o', 'LineWidth', 2, 'color', xcol(2,:));
legend([ha hs], {'GIA random model', 'Selected model'})
% legend([hp hg ha hs], {sprintf('q<%.2g, Poisson(p)', tgPoi), ...
%     sprintf('q>%.2g, %s > %d, %s(p, 1/p)', tgGamma, '\mu', tgMu, '\Gamma'), ...
%     'all other simulations', 'selected simulations'})
if length(indx) > 2
    saveas(hfig1, fullfile(outDir, 'map_qVSp_sim'))
end