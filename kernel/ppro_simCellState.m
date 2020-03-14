function [simOut, meanODE, Qgf] = ppro_simCellState(simModel, simOptions, simOut)
%   [simOut, meanODE, Qgf] = ppro_simCellState(simModel, simOptions, simOut)
% This function postprocesses the simulation outputs. It computes: the
% mean, variance, distribution, survival cases based on the statAnalysis
% function.
% simModel   IN: simulation model structure
% simOptions IN: simulation option structure
% simOut     IN: simulation output cell of structure
% simOut    OUT: ppro field is added to each element, containing the
%               statAnalysis function results.
% meanODE   OUT: results of the integration of the differential equation
%               describing the mean (of all)
% Qgf       OUT: analytical estimation of the extinction probability as
%               function of time
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

% parameters
Ntime = length(simOptions.time); Nrs = simOptions.Nrs;

% STATISTICS
outStat_all = cell(1,Nrs); outStat_surv = cell(1,Nrs);
for irs = 1:Nrs
    outStat_all{irs} = statAnalysis(simOut{irs}.xstate, simModel.indxActiveState, []);
    outStat_surv{irs} = statAnalysis(simOut{irs}.xstate, simModel.indxActiveState, outStat_all{irs}.indxSurv);
    if any(simOut{irs}.istatus==-1)
        fprintf('run %d: %d %% NaN at the end of simulation\n', irs, (1-length(outStat_all{irs}.indxNotNan{end})/simOptions.Nrun)*100)
    end
end

% get initial conditions for mean and extinction estimation
if isfield(simOptions.iniCondition, 'X0')
    X0 = simOptions.iniCondition.X0;
else
    X0 = [];
    for irs = 1:Nrs
        X0 = unique([X0 simOut{irs}.X0v]', 'row')';
    end
end

% MEAN (ODE)
eventFnc = @(t, x) odeEvent(t, x, simModel.indxActiveState, -1, 1e6);
options = odeset('Event', eventFnc);
if ~isfield(simOptions, 'iode') || simOptions.iode
    if isnumeric(simModel.odeM) % old input file
        [time_ode, mu_ode] = ode45(@(t, x) simModel.odeM*x, [simOptions.iniCondition.T0 simOptions.time(end)], mean(simOut{irs}.X0v, 2), options);
    else
        [time_ode, mu_ode] = ode45(simModel.odeM, [simOptions.iniCondition.T0 simOptions.time(end)], mean(simOut{irs}.X0v, 2), options);
    end
    meanODE = struct('time', time_ode', 'mean', mu_ode');
else
    meanODE = [];
end

% EXTINCTION
Qs = zeros(Nrs, Ntime);
Qgf = struct('time', {}, 'Q', {}, 'Qa', {});    
if (~isfield(simOptions, 'iext') || simOptions.iext) && ~isempty(simModel.extEst)
    % simulation
    for irs = 1:simOptions.Nrs
        Qs(irs,:) = squeeze(sum(sum(simOut{irs}.xstate(simModel.indxActiveState,:,:), 1)==0, 2))';
    end
    Qs = Qs/simOptions.Nrun;
    % from generating function
    for iX0 = 1:size(X0, 2)
        [time, Qn, Qa] = extProbTime(@(s) repGenFnc(s, simModel.extEst.p0, simModel.extEst.p1, simModel.extEst.p2), ...
            simModel.extEst.rtot, X0(:,iX0), simOptions.time(end));
        Q = 1;
        for is = 1:simModel.Nstate
            Q = Q.*Qn(is,:);
        end
        Qgf(iX0) = struct('time', time, 'Q', Q, 'Qa', Qa);
    end
end

% add to output
for irs = 1:Nrs
    ppro = struct('outStat_all', outStat_all{irs}, 'outStat_surv', outStat_surv{irs}, ...
        'Qsim', Qs(irs,:));
    simOut{irs}.ppro = ppro;
end