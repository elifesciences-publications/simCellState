function plot_distributionAll(xData, plotOptions, figName)

% plot distribution of multiple cases over the same graph

% options
simOptions = xData{1}.simOptions;
simOptions.plot = plotOptions;
simOptions.iode = 0;
simOptions.iext = 0;
simOptions.plot.title = 'Total number of cells (except death states), surviving clones';
simOptions.plot.style = '*-';
xcol = hsv(length(xData));

% loop on cases
for ifile = 1:length(xData)
    % postpro
    simOut = ppro_simCellState(xData{ifile}.simModel, simOptions, xData{ifile}.simOut);
    if ifile == 1
        simOptions.plot.isave = 0;
        simOptions.plot.figHandle = figure;
    else
        simOptions.plot.figHandle = openfig(figName);
        simOptions.time = xData{ifile}.simOptions.time; % in case they have different time step
        simOptions.plot.imeanTMultiple = 1;
    end
    simOptions.plot.indxtime = plotOptions.indxtime(simOptions.time);
    simOptions.plot.color = xcol(ifile,:);
    simOptions.plot.MarkerFace = xcol(ifile,:);
    simOptions.plot.MarkerEdge = xcol(ifile,:);
    indx_X0 = find(xData{ifile}.simOptions.iniCondition.X0);
    X0_txt = sprintf('n_{%d}(0) = %d', indx_X0, xData{ifile}.simOptions.iniCondition.X0(indx_X0));
    simOptions.plot.tag = X0_txt;
    plot_distribution(xData{ifile}.simModel, simOptions, simOut)
    
    % save
    saveas(gcf, figName)
    close(gcf)
    
end

return