clear all
close all

% Script to run the test case for the bimodal distribution.
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

% simulation parameters
Nrun = 1e4; xTime = 0:1:20; G = 1;
p_tg = 30; % Lambda_R hat
dx = 0.05;

% reference curve
xx0 = -2.5:0.01:2.5; 
hfig = figure; setFigureProp(hfig, 10); hold on; grid on
% hn = plot(xx0, normpdf(xx0, 1, sqrt(1/p_tg)), 'k', 'LineWidth', 2);
hn = plot(xx0, normpdf(xx0, 0, 1), 'k', 'LineWidth', 2);

dataFnc = @(x) [x.outParam.a, x.outParam.b, x.outParam.c, ...
    x.outParam.lambda1Eq, x.outParam.lambda2Eq, x.outParam.gammaEq, x.outParam.p, x.outParam.q, ...
    x.outParam.L1, x.outParam.L2, x.outParam.L3, x.outParam.O12, x.outParam.O21, x.outParam.G];
xdata = [];

% TEST PARAMETERS
xa = [1 1 10 10 10 10 10];
xomegaHat = p_tg*[1 1/1e3 10 1 1/10 1/1e2 1/1e3];
xq_tg = zeros(size(xa));
xb = ones(size(xa));
for ia = 1:length(xa)
    a = xa(ia); b = xb(ia); c = xomegaHat(ia); q_tg = xq_tg(ia);
    simOut(ia) = run_simGAD1(p_tg, q_tg, a, b, c, G, xTime, Nrun);
    xlegtc{ia} = sprintf('GIA^B#%d', ia);
    
    % plot
    mu_tgii = p_tg/(1-q_tg);
    if ia == length(xa) % bimodal is the last one
        hb = plot(simOut(ia).outParam.xxtilde, simOut(ia).outParam.gxxtilde, 'k-.', 'LineWidth', 2);
    end
    % [nn1, nnpmf] = getPmf(simOut(ia).nn_sim{1}-1, 1);
    nn_jj = simOut(ia).nn_sim{1}-1;
    [xx1, xxpdf] = getPdf(nn_jj, dx);
    % rescaling for the variance
    xmu = mean(nn_jj); xsigma = std(nn_jj);
    xx1 = (xx1-1)*xmu/xsigma;
    xxpdf = xxpdf*xsigma/xmu;
    
    htc(ia) = plot(xx1, xxpdf, '.', 'MarkerSize', 15);    
    xdata = [xdata; dataFnc(simOut(ia))];
end
xlabel('$\tilde{x}_C$'); ylabel('$P(\tilde{x}_C)$');
legend([hn hb htc], ['Normal', 'Bimodal', xlegtc])
set(gca, 'XLim', [-2.5 2.5])
saveas(hfig, fullfile(pwd, 'fig', 'testCase_bimodal'))

rmpath('../GenericNetwork')
