function [n, Pn] = integrMEq(nAMax, nBMax, mEqFnc, time, indxT)
%  [n, Pn] = integrMEq(nAMax, nBMax, time, lambdaii, gammaii, indxT)
% This function integrates the master equation.
% 
% nAMax, nBMax      IN: maximum number of cells in A and B considered, 
%                       nAMax*nBMax gives the length of the state vector
% mEqFnc            IN: master equation function handle
% time              IN: time step for integration
% indxT             IN: index at which the output is returned
% n                 OUT: total number of cell
% Pn                OUT: probability at time(indxT)
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

% state
na = 0:1:nAMax; nb = 0:1:nBMax; % max numbers of cell in state A and B
[Na, Nb] = meshgrid(na, nb);
% initial condition
Pv0 = zeros(size(Na)); Pv0(nb==0, na==1) = 1; Pv0 = Pv0(:);

% integration
[~, Pv] = ode45(@(t, p) mEqFnc(t, p, Na, Nb), time, Pv0);

% output
P = reshape(Pv(indxT,:), size(Na));
n = (na(1)+nb(1)):1:(na(end)+nb(end));
Pn = zeros(size(n));
N = Na + Nb;
for jj = 1:length(n)
    Pn(jj) = sum(P(N==n(jj)));
end
% surviving only (rescale everything)
n = n(2:end);
Pn = Pn(2:end)/sum(Pn(2:end));

end