clear all
close all

% This script creates all the possible combination of transition matrix Ct
% and division matrix Cd for a given number of states (N/Nstate). A maximum
% number of state transitions for each state can be set (Nstate).
% Output is saved in a file.
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
N = 4; % number of states of the SCC
Nstate = N+1; % total number of states (N + 1 outgoing)
sMax = Inf; % max number of transition for each state

% DIVISION
% --------
% building all possible columns
Cdc = zeros(Nstate,1);
for ii = 1:Nstate
    for jj = 1:Nstate
        Cdc0 = zeros(Nstate,1);
        Cdc0(ii) = 1;
        Cdc0(jj) = Cdc0(jj) + 1;
        Cdc = [Cdc Cdc0];
    end
end
Cdc = unique(Cdc', 'rows')';
% combining columns
Cd = Cdc;
if N>1
    for ii = 2:N
        Cd = combvec(Cd,Cdc);
    end
end
Cd = [Cd; zeros(Nstate, size(Cd,2))];
% permutation could be removed, but it is feasible only if Nstate < 4

% TRANSITIONS
% -----------
% building all possible columns
Ctc = de2bi(0:1:2^Nstate-1)';
Ctc(:,sum(Ctc)>sMax) = []; % max number of transitions
if N>=5 % number of combination must be reduced (to half)
    Ctc(:,rand(1,size(Ctc,2))>0.5) = [];
end
% combining columns
Ct = Ctc;
if N>1
    for ii = 2:N
        Ct = combvec(Ct,Ctc);
    end
end
Ct = [Ct; zeros(Nstate, size(Ct,2))];
% removing transition to itself
for ii = 1:N
    indx = sub2ind([Nstate, Nstate], ii, ii);
    Ct(:,Ct(indx,:) == 1) = [];
end
% removing if there are more than sMax transitions for each state
indxrm = [];
for ii = 1:size(Ct,2)
    Ctii = reshape(Ct(:,ii), Nstate, Nstate);
    if any(sum(Ctii) > sMax)
        indxrm = [indxrm ii];
    end
end
Ct(:,indxrm) = [];

% sort for the number of transitions and divisions
[~, icd] = sort(sum(Cd)); Cd = Cd(:,icd);
[~, ict] = sort(sum(Ct)); Ct = Ct(:,ict);

% save results
fileName = sprintf('rawDT_N%d', N);
dout = fullfile(pwd, 'SCC');
if ~exist(dout, 'dir')
    mkdir(dout)
end
save(fullfile(dout, fileName), 'Cd', 'Ct')
