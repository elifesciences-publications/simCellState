function [meanF, dmean, ext, mui0, MSE0, MSE2ap, MSE2bp] = plot_distributionCheck(xData, plotOptions, figName)

% check the type of distribution

% options
simOptions = xData{1}.simOptions;
simOptions.iode = 1;
simOptions.iext = 0;

% loop on cases: find the one with maximum mean
for ifile = 1:length(xData)
    % postpro
    [simOut, meanODE] = ppro_simCellState(xData{ifile}.simModel, simOptions, xData{ifile}.simOut);
    % max mean: save pdf and cdf
    irs = xData{ifile}.simOptions.plot.irs(1);
    mui = sum(simOut{irs}.ppro.outStat_surv.mean(xData{ifile}.simModel.indxActiveState,:));
    itime = find(~isnan(mui), 1, 'last'); % take the last point before extinction (or set a value for the extinction ex. 90%?)
    mui = mui(itime);
    if isempty(mui); mui = 0; end
    if ifile == 1
        xx = []; indx_ii = ifile;
        mui0 = mui; itime0 = itime;
        if ~isnan(mui) && mui ~= 0
            [xx, cdf, pdf, meanF, dmean, ext] = getData(simOut, meanODE, xData{ifile}, irs, itime);
        end
    else
        if mui > mui0 | (mui == mui0 & itime > itime0)
            indx_ii = ifile;
            mui0 = mui;
            [xx, cdf, pdf, meanF, dmean, ext] = getData(simOut, meanODE, xData{ifile}, irs, itime);
        end
    end
end

% distribution test: expo
if ~isempty(xx)
    
    % distribution
    pdf_ref = exppdf(xx, mui0);
    cdf_ref = expcdf(xx, mui0);
    
    % plot pdf and cdf
    MSE0 = immse(pdf, pdf_ref);
    figure; subplot(2,1,1); hold on; grid on
    plot(xx/mui0, pdf, '*:b');
    plot(xx/mui0, pdf_ref, 'k');
    xlabel('n/\mu'); ylabel('pdf(n)')
    legend(sprintf('simulation (MSE = %.2e)', MSE0), 'exponential distribution')
    % plot cdf
    MSE0 = immse(cdf, cdf_ref);
    subplot(2,1,2); hold on; grid on
    plot(xx/mui0, cdf, '*:b');
    plot(xx/mui0, cdf_ref, 'k');
    xlabel('n/\mu'); ylabel('cdf(n)')
    legend(sprintf('simulation (MSE = %.2e)', MSE0), 'exponential distribution')
    if plotOptions.isave
        saveas(gcf, [figName '_cdfpdf_' num2str(indx_ii)]); close(gcf)
    end
    
    % PDF/CDF fit
    thrs = 1/4;
    indx = find(pdf > 0); indxT = find(xx/mui0 > thrs & pdf > 0);
    figure; hold on; grid on;
    % pdf
    isp = {2, 2, 1};
    [MSE1ap, MSE1bp] = testExp(@(xx,pdf,mu) log(pdf*mu)./-xx*mu, xx, pdf, pdf_ref, mui0, indx, indxT, 'ln(pdf(n)\mu)/-n\mu', isp);
    isp = {2, 2, 3};
    [MSE2ap, MSE2bp] = testExp(@(xx,pdf,mu) (log(pdf*mu)./-xx*mu-1).*pdf, xx, pdf, pdf_ref, mui0, indx, indxT, 'ln(pdf(n)\mu)/-n\mupdf(n)', isp);
    % cdf
    isp = {2, 2, 2};
    % indx = find(cdf < 1); indxT = find(xx/mui0 > thrs & cdf < 1);
    indx = 1:length(xx)-1; indxT = find(xx(indx)/mui0 > thrs);
    [MSE1ac, MSE1bc] = testExp(@(xx,cdf,mu) log(1-cdf)./-xx*mu, xx, cdf, cdf_ref, mui0, indx, indxT, 'ln(1-cdf(n))/-n\mu-1', isp);
    isp = {2, 2, 4};
    [MSE2ac, MSE2bc] = testExp(@(xx,cdf,mu) (log(1-cdf)./-xx*mu-1).*(1-cdf), xx, cdf, cdf_ref, mui0, indx, indxT, '(ln(1-cdf(n))/-n\mu-1)(1-cdf(n))', isp);
    if plotOptions.isave
        saveas(gcf, [figName '_fit_' num2str(indx_ii)]); close(gcf)
    end
    
    % %     fnc1 =
    % %     MSE1a = immse(fnc1(xx(indx), cdf(indx),mui0), fnc1(xx(indx), cdf_ref(indx),mui0));
    % %     MSE1b = immse(fnc1(xx(indxT), cdf(indxT),mui0), fnc1(xx(indxT), cdf_ref(indxT), mui0));
    % %     figure; hold on; grid on
    % %     plot(xx/mui0, fnc1(xx, pdf, mui0), '*:b');
    % %     plot(xx/mui0, fnc1(xx, pdf_ref, mui0), 'k');
    % %     xlabel('n/Mean'); ylabel('ln(1-cdf(n))/-n')
    % %     legend(sprintf('simulation (MSE = %.2e/%.2e)', MSE1a, MSE1b), 'exponential distribution')
    %
    %     % CDF fit
    %
%         indx = 1:length(xx)-1; % remove last point (log(1) = 0)
% %         if length(indx) > 1
%             % plot cdf for fitting (I)
%             indxT = find(xx(indx)/mui0 > thrs);
%             fnc1 = @(xx,cdf,mu) log(1-cdf)./-xx*mu;
%             MSE1a = immse(fnc1(xx(indx), cdf(indx),mui0), fnc1(xx(indx), cdf_ref(indx),mui0));
%             MSE1b = immse(fnc1(xx(indxT), cdf(indxT),mui0), fnc1(xx(indxT), cdf_ref(indxT), mui0));
%             figure; hold on; grid on
%             plot(xx(indx)/mui0, fnc1(xx(indx), cdf(indx),mui0), '*:b');
%             plot(xx/mui0, fnc1(xx, cdf_ref, mui0), 'k');
%             xlabel('n/Mean'); ylabel('ln(1-cdf(n))/-n')
%             legend(sprintf('simulation (MSE = %.2e/%.2e)', MSE1a, MSE1b), 'exponential distribution')
%             if plotOptions.isave
%                 saveas(gcf, [figName '_cdfFit1_' num2str(indx_ii)]); close(gcf)
%             end
%             % plot cdf for fitting (II)
%             fnc1 = @(xx,cdf,mu) (log(1-cdf)./-xx*mu-1).*(1-cdf);
%             MSE2a = immse(fnc1(xx(indx), cdf(indx),mui0), fnc1(xx(indx), cdf_ref(indx),mui0));
%             MSE2b = immse(fnc1(xx(indxT), cdf(indxT),mui0), fnc1(xx(indxT), cdf_ref(indxT), mui0));
%             figure; hold on; grid on
%             plot(xx(indx)/mui0, fnc1(xx(indx), cdf(indx),mui0), '*:b');
%             plot(xx/mui0, fnc1(xx, cdf_ref, mui0), 'k');
%             xlabel('n/Mean'); ylabel('(ln(1-cdf(n))/-n - 1)*(1-cdf(n))')
%             legend(sprintf('simulation (MSE = %.2e/%.2e)', MSE2a, MSE2b), 'exponential distribution')
%             if plotOptions.isave
%                 saveas(gcf, [figName '_cdfFit2_' num2str(indx_ii)]); close(gcf)
%             end
    %     else
    %         MSE1a = NaN; MSE1b = NaN; MSE2a = NaN; MSE2b = NaN;
    %     end
else
    MSE0 = NaN; meanF = NaN; dmean = NaN; ext = NaN;
    MSE1ac = NaN; MSE1bc = NaN; MSE2ac = NaN; MSE2bc = NaN;
    MSE1ap = NaN; MSE1bp = NaN; MSE2ap = NaN; MSE2bp = NaN;
end

end

% BUILDIN FUNCTIONS

% get data
function [xx, cdf, pdf, meanF, dmean, ext] = getData(simOut, meanODE, xData, irs, itime)

xx = simOut{irs}.ppro.outStat_surv.distrT(itime).xx;
cdf = simOut{irs}.ppro.outStat_surv.distrT(itime).cdf;
pdf = simOut{irs}.ppro.outStat_surv.distrT(itime).pdf;
meanF = simOut{irs}.ppro.outStat_all.meanT(itime);
meanFode = interp1(meanODE.time, sum(meanODE.mean(xData.simModel.indxActiveState,:),1), xData.simOptions.time(itime));
dmean = abs(meanF - meanFode);
ext = 1-length(simOut{irs}.ppro.outStat_all.indxSurv{itime})/xData.simOptions.Nrun;

end

% plot fit
function [MSEa, MSEb] = testExp(fnc, xx, pcdf, pcdf_ref, mui0, indx, indxT, ylab, isbp)

if length(indx) > 1
    % compute errors
    MSEa = immse(fnc(xx(indx), pcdf(indx), mui0), fnc(xx(indx), pcdf_ref(indx), mui0));
    MSEb = immse(fnc(xx(indxT), pcdf(indxT), mui0), fnc(xx(indxT), pcdf_ref(indxT), mui0));
    % plot
    subplot(isbp{:}); hold on; grid on
    plot(xx(indx)/mui0, fnc(xx(indx), pcdf(indx), mui0), '*:b');
    plot(xx/mui0, fnc(xx, pcdf_ref, mui0), 'k');
    xlabel('n/\mu'); ylabel(ylab)
    legend(sprintf('simulation (MSE = %.2e/%.2e)', MSEa, MSEb), 'exponential distribution')
else
    MSEa = NaN;
    MSEb = NaN;
end
end