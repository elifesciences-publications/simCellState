clear all
close all

% This script run and process the simulations for the population asymmetry
% and compare the results in case metastates are modelled.

% METASTATES
% original model
% X1 -> X1+X1 (L1)
% X1 -> X1+X2 (L2)
% X1 -> X2+X2 (L3)
% metastate
% X1 -> X4 (Om1), X4 -> X1+X1 (Lm1)
% X1 -> X5 (Om2), X5 -> X1+X2 (Lm2)
% X1 -> X6 (Om3), X6 -> X2+X2 (Lm3)
%  6 unknown variables: Om1, Om2, Om3, Lm1, Lm2, Lm3
% L1/L2 = Om1/Om2 % to maintain the ratios (probability)
% L2/L3 = Om2/Om3
% L1 = 1/(1/Om1+1/Lm1) % to maintain the same total rate (!?!)
% L2 = 1/(1/Om2+1/Lm2)
% L3 = 1/(1/Om3+1/Lm3)
% Lm1 = Om1*Delta % to control the ratio, delta is a fixed parameter
% syms L1 L2 L3 Lm1 Lm2 Lm3 Om1 Om2 Om3 Delta positive
% ss = solve(L1/L2 - Om1/Om2, L2/L3 - Om2/Om3, L1 -1/(1/Om1+1/Lm1), L2 - 1/(1/Om2+1/Lm2), L3 - 1/(1/Om3+1/Lm3), Lm1 - Om1*Delta, Om1, Om2, Om3, Lm1, Lm2, Lm3);

% parameters
outDir = '../../';

% test case parameters
gamma = [1 1]; lambda = [1 2]; r = [1/4 1/6]; delta = [1/500 1/500]; Nsim = length(lambda);
X0 = 1; Nrun = 5e4;
xsimTg = {'Case 1', 'Case 2'};
xleg = cell(1, Nsim);
xcol = lines(Nsim);

% Run and plot results
for ii = 1:Nsim
    tg = dateTag(clock);
    
    % SIMULATION
    % ----------
    lambdaii = lambda(ii); rii = r(ii); gammaii = gamma(ii); deltaii = delta(ii);
        
    % population asymmetry - metastate
    [inputRawM, simOptionsM] = inputRead('PopulationAsymmetryM', {gammaii lambdaii rii X0 Nrun deltaii}, outDir);
    simModelM = inputSetup(inputRawM);
    % simulation
    tic
    simOutM = simCellStateLoop(simModelM, simOptionsM);
    toc
    
end
