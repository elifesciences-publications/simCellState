function plot_extinction(~, simOptions, simOut, Qgf)

% This function plot the extinction or the survival probability. The
% results of the integration of the differential equation are also shown.

% EXTNCTION/SURVIVABILITY
%------------------------
if simOptions.plot.iext && ~isempty(Qgf)
    figure; hold on; grid on
    for iX0 = 1:length(Qgf)
        h0 = plot(Qgf(iX0).time, Qgf(iX0).Q, 'k', 'linewidth', 2);
        if simOptions.plot.iext == -1 % survivability
            set(h0, 'Ydata', 1-get(h0, 'Ydata'))
        end
    end
    for irs = 1:simOptions.Nrs
        h1 = plot(simOptions.time, simOut{irs}.ppro.Qsim, '*-b');
        if simOptions.plot.iext == -1 % survivability
            set(h1, 'Ydata', 1-get(h1, 'Ydata'))
        end
   end
    h = [h0(1) h1(1)]; xleg = {'analytical solution', 'simulation'};
    xlabel('time (s)'); 
    if simOptions.plot.iext == -1 % survivability
    	ylabel('survival probability ()')
    else
        ylabel('extinction probability ()')
    end
    legend(h, xleg, 'Location', 'best')
    title(simOptions.plot.title, 'Interpreter','none')
    if simOptions.plot.isave
        if simOptions.plot.iext == -1 % survivability
            saveas(gcf, fullfile(simOptions.outFolder, ['ext_' simOptions.fileName]))
        else
            saveas(gcf, fullfile(simOptions.outFolder, ['surv_' simOptions.fileName]))
        end
    end    
end
