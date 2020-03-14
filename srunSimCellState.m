clear all
close all

% Main driver to run the stochastic simulations.
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

% INPUT PARAMETERS
% ----------------
icase = 'out_ncons_20190409150503'; Nrun = 5e4;
scase = {Nrun};
simCase = 'genericNetwork';
simInFolder = 'GENERIC';
outDir = pwd;

% SD0 = 123189; rng(SD0);
SD0 = 2*round(cputime*1e3)+1; rng(SD0) % for repeatibility
jr = rand(1,1000);
ir = 225;
simCellStateClone(icase, scase, simInFolder, simCase, outDir, ir, jr(ir:end));

save(fullfile(outDir, 'io', 'OUT', simInFolder, icase, 'out'), 'SD0')