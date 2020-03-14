function [Cd, Ct] = attachNetworkCritical(xOut_crit, Cdsub, Ctsub)
%  [Cd, Ct] = attachNetworkCritical(xOut_crit, Cdsub, Ctsub)
% This function attach the critical SCC to the subcritical part of the
% generic network.
% The diagonal block is filled first. Then for each outgoing component
% the incoming one is selected randomly among the nodes of the subcritical
% SCC. The subctitical matrix is added at the end.
% 
% xOut_crit IN: cell containing the structure with the critical SCC outputs
% Cdsub     IN: Division matrix with the subcritical SCC
% Ctsub     IN: Transition matrix with the subcritical SCC
% Cd    OUT: Division matrix
% Ct    OUT: Transition matrix
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

% state of subcricitcal SCC
Nstate_crit = size(xOut_crit.Cd,1)-1;
Nstate_sub = size(Cdsub,1);
Nstate_tot = Nstate_crit+Nstate_sub;

% initialization
Cd = zeros(Nstate_tot, Nstate_tot);
Ct = zeros(Nstate_tot, Nstate_tot);

% insert the critcal SCC
Cd(1:Nstate_crit, 1:Nstate_crit) = xOut_crit.Cd(1:Nstate_crit, 1:Nstate_crit);
Ct(1:Nstate_crit, 1:Nstate_crit) = xOut_crit.Ct(1:Nstate_crit, 1:Nstate_crit);

% now fill (randomly) the outgoing components
kk = 1:Nstate_crit;
jj = Nstate_crit+1:Nstate_tot-1; % avoind linking to death
% assign randomly the outgoing components for the cell division
Cd = fillOutgoingCd(Cd, jj, kk, xOut_crit.Cd);
% assign randomly the outgoing components for the cell transition
Ct = fillOutgoingCt(Ct, jj, kk, xOut_crit.Ct);

% fill remaining part of Ct and Cd
Cd(Nstate_crit+1:Nstate_tot, Nstate_crit+1:Nstate_tot) = Cdsub;
Ct(Nstate_crit+1:Nstate_tot, Nstate_crit+1:Nstate_tot) = Ctsub;

return