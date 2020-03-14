function plot_distribution(simModel, simOptions, simOut)

% This function plot the distribution of each state and the total at a given time.
% default options
if ~isfield(simOptions.plot, 'titleInterpreter')
    simOptions.plot.titleInterpreter = 'none';
end
if ~isfield(simOptions.plot, 'figHandle')
    figHandle = [];
else
    figHandle = simOptions.plot.figHandle;
end
[plotOpt, plotOpt1] = plotProperties(simOptions.plot);
if ~isfield(simOptions.plot, 'tag')
    tag = [];
else
    tag = simOptions.plot.tag;
end
if ~isfield(simOptions, 'timeUnit')
    timeUnit = 's';
else
    timeUnit = simOptions.timeUnit;
end

% distType = {'pdf', 'cdf'}; % plot pdf and cdf
distType = {'pdf'}; % plot pdf and cdf

for idist = 1:length(distType)
    
    for itime = simOptions.plot.indxtime
        
        % DISTRIBUTION OF EACH STATE
        %---------------------------
        if simOptions.plot.idist
            if isempty(figHandle)
                figure;
            else
                figure(figHandle)
            end
            hold on; grid on
            for istate = 1:simModel.Nstate
                subplot(simModel.Nstate,1,istate); hold on; grid on
                plotDistr(simModel, simOut, simOptions, itime, istate, simOptions.plot.idistScaled, plotOpt, plotOpt1, tag, distType{idist})
                % set(gca, 'YScale', 'log')
            end
            subplot(simModel.Nstate,1,1); title(sprintf('%s distribution at t = %.2g %s', simOptions.plot.title, simOptions.time(itime), timeUnit), 'Interpreter',simOptions.plot.titleInterpreter)
            % save
            if simOptions.plot.isave
                if simOptions.plot.idistScaled
                    saveas(gcf, fullfile(simOptions.outFolder, ['dist_' distType 'scl_t' num2str(itime) '_' simOptions.fileName]))
                else
                    saveas(gcf, fullfile(simOptions.outFolder, ['dist_' distType '_t' num2str(itime) '_' simOptions.fileName]))
                end
            end
        end
        
        % DISTRIBUTION OF THE SUM OF ALL THE ACTIVE STATES
        %-------------------------------------------------
        if simOptions.plot.idistT
            if isempty(figHandle)
                figure;
            else
                figure(figHandle)
            end
            hold on; grid on
            plotDistr(simModel, simOut, simOptions, itime, [], simOptions.plot.idistScaled, plotOpt, plotOpt1, tag, distType{idist})
            % set(gca, 'YScale', 'log')
            title(sprintf('%s, distribution at t = %.2g %s', simOptions.plot.title, simOptions.time(itime), timeUnit), 'Interpreter',simOptions.plot.titleInterpreter)
            % save
            if simOptions.plot.isave
                if simOptions.plot.idistScaled
                    saveas(gcf, fullfile(simOptions.outFolder, ['distT_' distType 'scl_t' num2str(itime) '_' simOptions.fileName]))
                else
                    saveas(gcf, fullfile(simOptions.outFolder, ['distT_' distType '_t' num2str(itime) '_' simOptions.fileName]))
                end
            end
        end
        
    end
    
end
end

% plot pdf/cdf
function plotDistr(simModel, simOut, simOptions, itime, istate, iscale, plotOpt, plotOpt1, tag, distType)

% plot distribution
for irs = simOptions.plot.irs
    [xx, yy, Nsmean, Nsstd] = getData(simModel, simOut, irs, istate, itime, distType);
    switch iscale
        case 1 % Mean
            xscl = 1/Nsmean; yscl = 100*Nsmean;
        otherwise
            xscl = 1; yscl = 1;
    end
    if ~isnan(xx)
        hs = plot(xx*xscl, yy*yscl, plotOpt{:});
        switch simOptions.plot.idistRef
            case 'exp'
                if strcmp(distType, 'pdf')
                    href = plot(xx*xscl, exppdf(xx, Nsmean)*yscl, plotOpt1{:});
                else
                    href = plot(xx*xscl, expcdf(xx, Nsmean)*yscl, plotOpt1{:});
                end
                legend([hs href], ['sim' tag], sprintf('exp(%.2g)', Nsmean)')
            case 'geo'
                if strcmp(distType, 'pdf')
                    ref_geo = @(x, mu) (1-(1/mu)).^(x-1)*(1/mu); % only for surviving clones here (mean >= 1)
                else
                    ref_geo = @(x, mu) 1-(1-(1/mu)).^x; % only for surviving clones here (mean >= 1)
                end
                href = plot(xx*xscl, ref_geo(xx, Nsmean)*yscl, plotOpt1{:});
                legend([hs href], ['sim' tag], sprintf('geo(%.2g)', Nsmean))
            case 'poi'
                if strcmp(distType, 'pdf')
                    href = plot(round(min(xx):1:max(xx))*xscl, poisspdf(round(min(xx):1:max(xx)), Nsmean)*yscl, plotOpt1{:}); % valid for integer only
                else
                    href = plot(round(min(xx):1:max(xx))*xscl, poisscdf(round(min(xx):1:max(xx)), Nsmean)*yscl, plotOpt1{:}); % valid for integer only
                end
                legend([hs href], ['sim' tag], sprintf('poi(%.2g)', Nsmean))
            case 'nor'
                if strcmp(distType, 'pdf')
                    href = plot(xx*xscl, normpdf(xx, Nsmean, Nsstd)*yscl, plotOpt1{:});
                else
                    href = plot(xx*xscl, normcdf(xx, Nsmean, Nsstd)*yscl, plotOpt1{:});
                end
                legend([hs href], ['sim' tag], sprintf('norm(%.2g, %.2g)', Nsmean, Nsstd))
            otherwise
                if ~isempty(tag)
                    legend(hs, ['sim' tag])
                end
        end
        legend('off'); legend
    end
end
if isempty(istate) % sum
    xlab = 'N'; ylab = 'N';
else % each state
    xlab = ['n_' num2str(istate)]; ylab = ['n_' num2str(istate)];
end
switch iscale
    case 1 % Mean
        xlabel([xlab '/Mean'])
        ylabel(['frequency(' ylab ') (%)'])
    otherwise
        xlabel(xlab)
        ylabel([distType '(' ylab ')'])
end

end

% get data
function [xx, yy, Nsmean, Nsstd] = getData(simModel, simOut, irs, istate, itime, distType)
% simOut{irs}.ppro.outStat_surv = simOut{irs}.ppro.outStat_all;
if isempty(istate) % sum
    Nsmean = sum(simOut{irs}.ppro.outStat_surv.mean(simModel.indxActiveState,itime));
    C = cov(simOut{irs}.xstate(:,:,itime)');
    Nsstd = sqrt(sum(sum(C(simModel.indxActiveState, simModel.indxActiveState))));
    if ~isnan(Nsmean) && Nsmean ~= 0
        xx = simOut{irs}.ppro.outStat_surv.distrT(itime).xx;
        yy = simOut{irs}.ppro.outStat_surv.distrT(itime).(distType);
    else
        xx = NaN; yy = NaN;
    end
else % each state
    Nsmean = simOut{irs}.ppro.outStat_surv.mean(istate,itime);
    Nsstd = simOut{irs}.ppro.outStat_surv.std(istate,itime);
    if ~isnan(Nsmean)
        xx = simOut{irs}.ppro.outStat_surv.distr(istate,itime).xx;
        yy = simOut{irs}.ppro.outStat_surv.distr(istate,itime).(distType);
    else
        xx = NaN; yy = NaN;
    end
end
end


% plot properties
function [plotOpt, plotOpt1] = plotProperties(plotOptStr)

% default plot properties
plotOpt = {'color', 'b', 'LineStyle', ':', 'Marker', '*', 'MarkerFace', 'b', 'MarkerEdge', 'b'};
for ii = 1:2:length(plotOpt)
    if isfield(plotOptStr, plotOpt{ii})
        plotOpt{ii+1} = plotOptStr.(plotOpt{ii});
    end
end
plotOpt1 = {'color', 'k', 'LineStyle', '-', 'Marker', 'none', 'MarkerFace', 'k', 'MarkerEdge', 'k'};
for ii = 1:2:length(plotOpt)
    if isfield(plotOptStr, plotOpt{ii})
        plotOpt{ii+1} = plotOptStr.(plotOpt{ii});
    end
end

end