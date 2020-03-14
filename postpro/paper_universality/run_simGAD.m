function [simOut] = run_simGAD(p, q, G, xTime, Nrun)

% This function run the simulation and postprocess the results for the GAD
% test case model, given the parameters p and q.
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

% model parameter
L1 = p*G; L2 = q*G;
r = [L1 L2 G]; r = r(r>0); r_min = min(r);

% setup simulation
Nstate = 3; X0 = zeros(Nstate,1); X0(1) = 1; T0 = 0;
inputRaw.inputT = [G 2 3]; inputRaw.inputS = [L1 1 1 2; L2 2 2 2];
inputRaw.Nstate = Nstate;
simOptions = struct('Nrun', Nrun, 'Nrs', 1, 'Nrand', 10e4, ...
    'seed0', 2*round(cputime*1e3)+1, ...
    'iniCondition', struct('X0', X0, 'T0', T0), 'time', xTime/r_min, ...
    'isave', 0, ...
    'simPar', 1, 'mex', 1, 'plot', []);
simModel = inputSetup(inputRaw);
simOut = simCellStateLoop(simModel, simOptions);

% postprocessing
simOut = ppro_simPaper(simModel, simOptions, simOut, 0);

end