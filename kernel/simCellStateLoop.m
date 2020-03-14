function simOut = simCellStateLoop(simModel, simOptions)
%  simOut = simCellStateLoop(simModel, simOptions)
% This function is a wrapper for the simCellState function. It runs two
% levels of loop over Nrs and Nrun (repeated simulation and number of
% simulations).
% The definition of the inital seed in each single run is computed based on
% the randSeedInitialization function.
% Results are optionally saved.
%  simModel     IN: simulation model structure
%  simOptions   IN: simulation option structure
%  simOut      OUT: output structure. Cell of size (1,Nrs)
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

% simulation options
Nrun = simOptions.Nrun;
Nrs = simOptions.Nrs;
xtime = simOptions.time;
Ntime = length(xtime);
Nrand = simOptions.Nrand;
T0 = simOptions.iniCondition.T0;
if isfield(simOptions.iniCondition, 'X0') % deterministic, scalar
    X0 = simOptions.iniCondition.X0;
    pX = [];
    X0v = repmat(X0, 1, Nrun);
elseif isfield(simOptions.iniCondition, 'X0v') % deterministic, vectorial
    X0v = simOptions.iniCondition.X0v;
    pX = [];
else % initial condition based on a distribution
    pX = simOptions.iniCondition.pX;
    X0T = simOptions.iniCondition.X0T;
end
if ~isfield(simOptions, 'mex')
    simOptions.mex = 0;
end
if simOptions.mex == 1
    fun = @wrapSimCellStateMex;
else
    fun = @simCellState;
end

% model parameters
Nstate = simModel.Nstate;
indxState = simModel.indxState;
sRate = simModel.sRate;
deltaState = simModel.deltaState;
indxActiveState = simModel.indxActiveState;

% random seed initialization
NrunMax = 150e3*200; % hardcoded parameter, for higher values memory issues/slow
[SDs, SDt, SDx] = randSeedInitialization(simOptions.seed0, Nrun, Nrs, ~isempty(pX), simOptions.simPar, NrunMax/Nrand);

% initialization
simOut = cell(1, Nrs);
for irs = 1:Nrs % loop on Nrs
    
    % setup initial condition
    if ~isempty(pX)
        rng(SDx(irs));
        X0v = zeros(Nstate, Nrun);
        for iX0T = 1:X0T
            xrs = rand(X0T, Nrun);
            for irun = 1:Nrun
                ixr = find(xrs(irun) - cumsum(pX) < 0, 1, 'first');
                X0v(ixr,irun) = X0v(ixr,irun) +1;
            end
        end
    end
    
    % random initialization
    if Nrun < NrunMax/Nrand % vectorized option
        rng(SDs(irs),'twister');
        xus = rand(Nrun, Nrand); % state change
        rng(SDt(irs),'twister');
        xut = rand(Nrun, Nrand); % time
        ivec = 1;
    else
        ivec = 0;
    end
    
    % initalization
    xstate = NaN(Nstate, Nrun, Ntime);
    timeExt = NaN(1,Nrun);
    istatus = NaN(1,Nrun);
    if simOptions.simPar > 0
        if ivec % vectorized
            parfor irun = 1:Nrun % loop on Nrun
                [xstate(:,irun,:), timeExt(irun), istatus(irun)] = fun(Nstate, X0v(:,irun), T0, xtime, ...
                    xut(irun, :), xus(irun, :), ...
                    sRate, indxState, deltaState, indxActiveState);
            end
        else
            parfor irun = 1:Nrun % loop on Nrun
                % initialize random
                rng(SDs(irun, irs),'twister');
                us = rand(1, Nrand); % state change
                rng(SDt(irun, irs),'twister');
                ut = rand(1, Nrand); % time
                [xstate(:,irun,:), timeExt(irun), istatus(irun)] = fun(Nstate, X0v(:,irun), T0, xtime, ...
                    ut, us, ...
                    sRate, indxState, deltaState, indxActiveState)
            end
        end
    else
        if ivec % vectorized random
            for irun = 1:Nrun % loop on Nrun
                [xstate(:,irun,:), timeExt(irun), istatus(irun)] = fun(Nstate, X0v(:,irun), T0, xtime, ...
                    xut(irun, :), xus(irun, :), ...
                    sRate, indxState, deltaState, indxActiveState);
            end
        else
            for irun = 1:Nrun % loop on Nrun
                % initialize random
                rng(SDs(irun, irs),'twister');
                us = rand(1, Nrand); % state change
                rng(SDt(irun, irs),'twister');
                ut = rand(1, Nrand); % time
                [xstate(:,irun,:), timeExt(irun), istatus(irun)] = fun(Nstate, X0v(:,irun), T0, xtime, ...
                    ut, us, ...
                    sRate, indxState, deltaState, indxActiveState);
            end
        end
    end
    simOut{irs} = struct('xstate', xstate, 'timeExt', timeExt, 'istatus', istatus, 'X0v', X0v);
    
    % save results
    if simOptions.isave
        outFolder = simOptions.outFolder;
        if ~exist(outFolder, 'dir')
            mkdir(outFolder)
        end
        save(fullfile(outFolder, simOptions.fileName), 'simOut', 'simModel', 'simOptions', '-v7.3')
    end
    
end