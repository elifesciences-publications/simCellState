function plot_mean(simModel, simOptions, simOut, meanODE)

% This function plot the mean of each state and the total. The results of
% the integration of the differential equation for the mean are also shown.

% default options
if ~isfield(simOptions.plot, 'titleInterpreter')
    simOptions.plot.titleInterpreter = 'none';
end
if ~isfield(simOptions.plot, 'figHandle')
    figHandle = [];
else
    figHandle = simOptions.plot.figHandle;
end
if ~isfield(simOptions.plot, 'color')
    color = 'b'; color1 = 'r';
else
    color = simOptions.plot.color;
    color1 = color;
end
if ~isfield(simOptions.plot, 'style')
    style = '*:'; style1 = 'x:';
else
    style = simOptions.plot.style;
    style1 = style;
end
if ~isfield(simOptions.plot, 'tag')
    tag = '';
else
    tag = simOptions.plot.tag;
end
if ~isfield(simOptions.plot, 'imeanTMultiple')
    simOptions.plot.imeanTMultiple = 0;
end
if ~isfield(simOptions.plot, 'imeanTsubplot')
    simOptions.plot.imeanTsubplot = 0;
end
if ~isfield(simOptions, 'timeUnit')
    timeUnit = 's';
else
    timeUnit = simOptions.timeUnit;
end

% MEAN OF EACH STATE
%-------------------
if simOptions.plot.imean
    if isempty(figHandle)
        figure;
    else
        figure(figHandle)
    end
    hold on; grid on
    for istate = 1:simModel.Nstate
        subplot(simModel.Nstate,1,istate); hold on; grid on
        xlabel(['time (' timeUnit ')']); ylabel(['mean(n_' num2str(istate) ')']);
        for iX0 = 1:length(meanODE)
            h0 = plot(meanODE(iX0).time, meanODE(iX0).mean(istate,:), 'k', 'linewidth', 2); % ode
        end
        for irs = simOptions.plot.irs
            h1 = plot(simOptions.time, simOut{irs}.ppro.outStat_all.mean(istate,:), style, 'color', color); % sim
            h2 = plot(simOptions.time, simOut{irs}.ppro.outStat_surv.mean(istate,:), style1, 'color', color1); % sim
        end
        if ~isempty(meanODE)
            h = [h0(1) h1(1) h2(1)]; xleg = {'ode', 'sim (all)', 'sim (surv)'};
        else
            h = [h1(1) h2(1)]; xleg = {'sim (all)', 'sim (surv)'};
        end
        legend(h, xleg, 'Location', 'best')
    end
    % add title
    subplot(simModel.Nstate,1,1); title(simOptions.plot.title, 'Interpreter','none')
    % save figure
    if simOptions.plot.isave
        saveas(gcf, fullfile(simOptions.outFolder, ['mean_' simOptions.fileName]))
    end
end

% MEAN OF THE SUM OF ALL THE ACTIVE STATES
%-----------------------------------------
if simOptions.plot.imeanT
    if isempty(figHandle)
        figure;
    else
        figure(figHandle)
    end
    if simOptions.plot.imeanTsubplot == 1
        subplot(2,1,1); hold on; grid on
        xlabel(['time (' timeUnit ')']); ylabel('mean all (n)');
        subplot(2,1,2); hold on; grid on
        xlabel(['time (' timeUnit ')']); ylabel('mean surviving (n)');
    else
        hold on; grid on
        xlabel(['time (' timeUnit ')']); ylabel('mean(n)');
    end
    if simOptions.plot.imeanTsubplot == 1
        subplot(2,1,1)
    end
    for iX0 = 1:length(meanODE)
        h0 = plot(meanODE(iX0).time, sum(meanODE(iX0).mean(simModel.indxActiveState,:), 1), 'k', 'linewidth', 2); % ode
    end
    for irs = simOptions.plot.irs
        if simOptions.plot.imeanTsubplot == 1
            subplot(2,1,1)
        end
        h1 = plot(simOptions.time, sum(simOut{irs}.ppro.outStat_all.mean(simModel.indxActiveState,:), 1), style, 'color', color); % sim
        if simOptions.plot.imeanTsubplot == 1
            subplot(2,1,2)
        end
        h2 = plot(simOptions.time, sum(simOut{irs}.ppro.outStat_surv.mean(simModel.indxActiveState,:), 1), style1, 'color', color1); % sim
    end
    if ~isempty(meanODE)
        if simOptions.plot.imeanTsubplot == 1
            subplot(2,1,1)
            h = [h0(1) h1(1)]; xleg = {'ode', ['sim ' tag]};
            if simOptions.plot.imeanTMultiple == 1 % turn off ode
                set(get(get(h0(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        else
            h = [h0(1) h1(1) h2(1)]; xleg = {'ode', 'sim (all)', 'sim (surv)'};
        end
    else
        if simOptions.plot.imeanTsubplot == 1
            subplot(2,1,1)
            h = []; xleg = {};
        else
            h = [h1(1) h2(1)]; xleg = {'sim (all)', 'sim (surv)'};
        end
    end
    legend(h, xleg, 'Location', 'best'); legend('off'); legend;
    % add title
    title(simOptions.plot.title, 'Interpreter',simOptions.plot.titleInterpreter)
    % save figure
    if simOptions.plot.isave
        saveas(gcf, fullfile(simOptions.outFolder, ['meanT_' simOptions.fileName]))
    end
    
end