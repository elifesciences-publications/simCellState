function [maxMSE, dev_exp, dev_geo, dev_poi, ext_f] = plot_meanAll(xData, plotOptions, figName)
%  [maxMSE, dev_exp, dev_geo, dev_poi, ext_f] = plot_meanAll(xData, plotOptions, figName)
% This function plots the mean of multiple cases in the same figure. It
% also returns some summary numerical outputs.
% 
% xData     IN: simulation input
% plotOptions IN: plot options structure
% figName   IN: figure name
% maxMSE    OUT: maximum mean squared error of the simulations wrt the mean value
% dev_<xx>  OUT: maximum deviation of the variance from epo, geo and poisson distributions
% ext_f     OUT: minimum extinction

% options
simOptions = xData{1}.simOptions;
simOptions.plot = plotOptions;
simOptions.iode = 1;
simOptions.iext = 0;
simOptions.plot.title = 'Total number of cells (except death states)';
simOptions.fileName = 'mean_all';
simOptions.plot.style = '*-';
xcol = hsv(length(xData));
tol = 1e-5;

% initialization
MSE = zeros(1, length(xData));
mean_surv_f = zeros(1, length(xData));
var_surv_f = zeros(1, length(xData));
ext_f = zeros(1, length(xData));
% loop on cases
for ifile = 1:length(xData)
    % postpro
    [simOut, meanODE] = ppro_simCellState(xData{ifile}.simModel, simOptions, xData{ifile}.simOut);
    meanSim = simOut{simOptions.plot.irs}.ppro.outStat_all.meanT;
    meanOdeInt = interp1(meanODE.time, sum(meanODE.mean(xData{ifile}.simModel.indxActiveState,:),1), simOptions.time);
    if meanOdeInt(end) >= tol % relative error
        meanSim = meanSim./meanOdeInt;
        meanOdeInt = ones(size(meanOdeInt));
    end
    MSE(ifile) = immse(meanSim, meanOdeInt);
    surv_i = cellfun(@(x) length(x),simOut{simOptions.plot.irs}.ppro.outStat_surv.indxSurv); 
    ext_f(ifile) = 1-surv_i(end)/simOptions.Nrun;
    if ext_f(ifile) < 0.99
        mean_surv_f(ifile) = simOut{simOptions.plot.irs}.ppro.outStat_surv.meanT(end);
        var_surv_f(ifile) = (simOut{simOptions.plot.irs}.ppro.outStat_surv.stdT(end))^2;
    else
        mean_surv_f(ifile) = 0; var_surv_f(ifile) = 0;
    end
    if ifile == 1
        simOptions.plot.isave = 0;
        simOptions.plot.figHandle = figure;
    else
        simOptions.plot.figHandle = openfig(figName);
        simOptions.time = xData{ifile}.simOptions.time; % in case they have different time step
        simOptions.plot.imeanTMultiple = 1;
    end
    simOptions.plot.color = xcol(ifile,:);
    indx_X0 = find(xData{ifile}.simOptions.iniCondition.X0);
    X0_txt = sprintf('n_{%d}(0) = %d', indx_X0, xData{ifile}.simOptions.iniCondition.X0(indx_X0));    
    simOptions.plot.tag = X0_txt;
    plot_mean(xData{ifile}.simModel, simOptions, simOut, meanODE)
    % save
    saveas(gcf, figName)
    close(gcf)
end

% return the maximum
maxMSE = max(MSE);

% mean wrt variance (for different distributions)
var_exp = mean_surv_f.^2;
var_geo = (1-1./mean_surv_f).*mean_surv_f.^2;
var_poi = mean_surv_f;
% return the maximum deviation of the variance normalized
mean_surv_f(mean_surv_f<tol) = 1;
dev_exp = max(abs(var_exp-var_surv_f)./mean_surv_f);
dev_geo = max(abs(var_geo-var_surv_f)./mean_surv_f);
dev_poi = max(abs(var_poi-var_surv_f)./mean_surv_f);
ext_f = min(ext_f); % minimum extinction

return