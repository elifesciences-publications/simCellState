clear all
close all

% This script run and process the simulations for the asymmetric division
% case.
% 1) Tissue dynamics --> only mean VS ODE
% 2) Clone dynamics --> mean VS ODE and distributio VS Master Equation

% parameters
outDir = '../../';
irun = 1;
simCase = 'AsymmetricDivision';

if irun
    gamma = [1 1 1]; lambda = [1 2 5]; Nsim = length(lambda); 
    X0 = 1; Nrun = 1e4;
    % X0 = 1000; Nrun = 1e4;
else
    simOutDir = {'AsymmetricDivision_20190712130850', 'AsymmetricDivision_20190712130905', 'AsymmetricDivision_20190712130920'}; Nsim = length(simOutDir);
end
xleg = cell(1, Nsim);
xcol = lines(Nsim);

% Run and plot results
for ii = 1:Nsim
    
    if irun
        % SIMULATION
        % ----------
        lambdaii = lambda(ii); gammaii = gamma(ii);
        [inputRaw, simOptions] = inputRead(simCase, {gammaii lambdaii X0 Nrun}, outDir);
        simModel = inputSetup(inputRaw);
        % simulation
        tic
        simOut = simCellStateLoop(simModel, simOptions);
        toc
    else
        load(fullfile(outDir, 'io', 'OUT', simCase, simOutDir{ii}, simOutDir{ii}))
        lambdaii = simModel.inputRaw.inputS(1,1);
        gammaii = simModel.inputRaw.inputT(1,1);
        X0 = simOptions.iniCondition.X0(1);
    end
    % POSTPRO
    % -------
    [simOut, meanODE, Qgf] = ppro_simCellState(simModel, simOptions, simOut);
    % check if at final time I have NaN
    iok = cellfun(@(x) length(x), simOut{1}.ppro.outStat_surv.indxNotNan);
    indxT = find(iok == iok(1), 1, 'last');
    if any(iok<iok(1))
        disp('check NaN')
    end
    
    % PLOT
    % ----
    tscale = min([simModel.inputRaw.inputS(:,1); simModel.inputRaw.inputT(:,1)]);
    if gammaii<1
        xleg{ii} = sprintf('%s = %d s^{-1}, %s = 1/%d s^{-1}', '\lambda', lambdaii, '\gamma', 1/gammaii);
    else
        xleg{ii} = sprintf('%s = %d s^{-1}, %s = %d s^{-1}', '\lambda', lambdaii, '\gamma', gammaii);
    end
    
    % plot network
    sysProp = sysProperties('TL', simModel.inputRaw.inputT, simModel.inputRaw.inputS, simModel.inputRaw.Nstate);
    figure; plotNetwork(sysProp, '');
    
    % mean of surviving
    figure(10); hold on; grid on
    plot(meanODE.time*tscale, sum(meanODE.mean(1:2,:))/X0, 'color', xcol(ii,:)) 
    plot(simOptions.time*tscale, simOut{1}.ppro.outStat_surv.meanT/X0, '*', 'color', xcol(ii,:))
    
    % distribution at Tf
    if X0 == 1
        
        % integration of the master equation
        na = 0:1:3; nb = 0:1:25;
        [Na, Nb] = meshgrid(na, nb);
        Pv0 = zeros(size(Na)); Pv0(nb==0, na==1) = 1; Pv0 = Pv0(:);
        [time, Pv] = ode45(@(t, p) mEq_asymmetric(t, p, Na, Nb, lambdaii, gammaii), simOptions.time, Pv0);
        indx = length(time); 
        P = reshape(Pv(indx,:), size(Na));
        n = (na(1)+nb(1)):1:(na(end)+nb(end));
        Pn = zeros(size(n));
        N = Na + Nb;
        for jj = 1:length(n)
            Pn(jj) = sum(P(N==n(jj)));
        end
        
        % plot
        figure(20); hold on; grid on
        plot(n, Pn, 'color', xcol(ii,:))
        plot(simOut{1}.ppro.outStat_surv.distrT(indxT).xx, simOut{1}.ppro.outStat_surv.distrT(indxT).pdf, '*', 'color', xcol(ii,:))
    end
end

figPath = fullfile(outDir, 'io', 'OUT', simCase, 'fig');
if ~exist(figPath, 'dir')
    mkdir(figPath)
end

% add legend/labels and save
tag = dateTag(clock);
figure(1)
saveas(gcf, fullfile(figPath, 'Network'))
figure(10); % legend(xleg)
xlabel('\tau [-]'); ylabel('\mu/\mu_0')
if X0 == 1
    saveas(gcf, fullfile(figPath, ['cloneMean_' tag]))
else
    saveas(gcf, fullfile(figPath, ['tissueMean' tag]))
end
if X0 == 1
    figure(20); % legend(xleg)
    xlabel('n'); ylabel('pdf(n)')
    saveas(gcf, fullfile(figPath, ['cloneDistr' tag]))
end