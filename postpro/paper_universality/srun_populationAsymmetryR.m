clear all
close all

% This script run and process the simulations for the population asymmetry
% case.
% 1) Tissue dynamics --> only mean VS ODE
% 2) Clone dynamics --> mean VS ODE and distributio VS Master Equation

% parameters
outDir = '../../';
irun = 1;
simCase = 'PopulationAsymmetryR';

if irun
    gamma = [1 1 1]; lambda = [1 2 2]; r = [1/4 1/4 1/6] ; Nsim = length(lambda); 
    X0 = 1; Nrun = 5e4;
    % X0 = 1000; Nrun = 1e4; % tissue dynamics
else
    simOutDir = {''}; Nsim = length(simOutDir);
end
xleg = cell(1, Nsim);
xcol = lines(Nsim);

% Run and plot results
for ii = 1:Nsim
    
    if irun
        % SIMULATION
        % ----------
        lambdaii = lambda(ii); rii = r(ii); gammaii = gamma(ii);
        [inputRaw, simOptions] = inputRead(simCase, {gammaii lambdaii rii X0 Nrun}, outDir);
        simModel = inputSetup(inputRaw);
        % simulation
        tic
        simOut = simCellStateLoop(simModel, simOptions);
        toc
    else
        load(fullfile(outDir, 'io', 'OUT', simCase, simOutDir{ii}, simOutDir{ii}))
        lambdaii = simModel.inputRaw.inputS(1,1);
        omegaii = simModel.inputRaw.inputT(2,1);
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
    tscale = simModel.inputRaw.addStr.tscale;
    xleg{ii} = sprintf('%s = %d s^{-1}, %s = %d s^{-1}, %s = 1/%d', '\lambda', lambdaii, '\gamma', gammaii, 'r', 1/rii);
    
    % plot network
    sysProp = sysProperties('TL', simModel.inputRaw.inputT, simModel.inputRaw.inputS, simModel.inputRaw.Nstate);
    figure; plotNetwork(sysProp, '');
    
    % check extinction
    figure; hold on; grid on
    plot(Qgf.time*tscale, Qgf.Q, 'color', xcol(ii,:)) 
    plot(simOptions.time*tscale, simOut{1}.ppro.Qsim, '*', 'color', xcol(ii,:))
    
    xlabel('\tau [-]'); ylabel('Q [-]')
    % mean of surviving
    Q = interp1(Qgf.time, Qgf.Q, meanODE.time);
    figure(10); hold on; grid on
    if X0 > 1
        plot(meanODE.time*tscale, sum(meanODE.mean(1:2,:))./(1-Q)/X0, 'color', xcol(ii,:)) % this takes into account the extinction
    end
    plot(simOptions.time*tscale, simOut{1}.ppro.outStat_surv.meanT/X0, '*', 'color', xcol(ii,:))
    
    % distribution at Tf
    if X0 == 1
        
        % integration of the master equation
        % na = 0:1:100; nb = 0:1:90;
        na = 0:1:50; nb = 0:1:50;
        [Na, Nb] = meshgrid(na, nb);
        Pv0 = zeros(size(Na)); Pv0(nb==0, na==1) = 1; Pv0 = Pv0(:);
        [time, Pv] = ode45(@(t, p) mEq_populationR(t, p, Na, Nb, lambdaii, rii, gammaii), simOptions.time, Pv0);
        
        indxTime = [3 7 10];
        xcold = lines(length(indxTime));
        for kk = 1:length(indxTime)
            indxT = indxTime(kk);
            indx = indxT;
            P = reshape(Pv(indx,:), size(Na));
            n = (na(1)+nb(1)):1:(na(end)+nb(end));
            Pn = zeros(size(n));
            N = Na + Nb;
            for jj = 1:length(n)
                Pn(jj) = sum(P(N==n(jj)));
            end            
            figure(20); hold on; grid on
            plot(n(2:end), Pn(2:end)/sum(Pn(2:end)), 'color', xcold(kk,:))
            plot(simOut{1}.ppro.outStat_surv.distrT(indxT).xx, simOut{1}.ppro.outStat_surv.distrT(indxT).pdf, '*', 'color', xcold(kk,:))
            mu = sum(n(2:end).*Pn(2:end)/sum(Pn(2:end)));
            plot(n(2:end), 1/mu*exp(-1/mu*n(2:end)), ':', 'color', xcold(kk,:)) % expo ref
            set(gca,'YScale','log', 'YLim', [1e-4 1])
            
        end
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
    saveas(gcf, fullfile(figPath, ['tissueMean_' tag]))
end
if X0 == 1
    figure(20); % legend(xleg)
    xlabel('n'); ylabel('pdf(n)')
    saveas(gcf, fullfile(figPath, ['cloneDistr_' tag]))
end