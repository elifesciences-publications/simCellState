function [xNout, tExt, istatus] = simCellState(Nstate, X0, T0, xtime, ut, us, sRate, indxState, deltaState, indxActiveState)
%  [xNout, tExt, istatus] = simCellState(Nstate, X0, T0, xtime, Ntime, ut, us, sRate, indxState, deltaState, indxActiveState)
% This function simulate a stochastic process. Process is stopped at the
% final time, in case of extinction or in case the maximum number of steps
% is reached (Nmax).
% Nstate    IN: number of states
% X0        IN: initial conditions (N,1) 
% T0        IN: initial time
% xtime     IN: simulation time steps (1,M)
% ut, us    IN: random numbers for time and for state (1,Nmax)
% sRate     IN: specific rates of each possible event (1, L)
% indxState IN: index of the state corresponding to the origin of each event (1,L)
% deltaState IN: state variation for each event (N,L)
% indxActiveState IN: index of active states
% xNout       OUT: state (N,M)
% tExt     OUT: extinction time. NaN if it extinctioni is not reached
% istatus  OUT: simulation status. If -1, NNmax has been reached
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

% initialization
Nrand = length(ut);
% Ntime = length(xtime);
% xNout = NaN(Nstate, Ntime);
xT = NaN(1, Nrand);
xN = NaN(Nstate, Nrand);
tExt = NaN;
istatus = 1; % check final condition. if -1, Nrand are not sufficient to cover the range of time

% loop on time
x0 = X0; t0 = T0;
xT(1) = T0; xN(:,1) = X0;
for ii = 1:Nrand
    % total rates
    if isnumeric(sRate)
        rate = sRate.*x0(indxState); % constant rates
    else
        rate = sRate(x0).*x0(indxState); % feedback
    end
    rtot = sum(rate);
    % compute time next event
    deltaT = log(1./ut(ii))/rtot;
    if t0 + deltaT > xtime(end)
        ii = ii-1;
        break
    end
    t1 = t0 + deltaT;
    % next event
    x1 = x0 + deltaState(:,find(us(ii) - cumsum(rate)/rtot < 0, 1, 'first')); % checking for 0 or 1...
    xT(ii+1) = t1; xN(:,ii+1) = x1;
    if all(x1(indxActiveState)==0)
        tExt = t1;
        t1=xtime(end);
        break
    end
    % reinitialize
    t0 = t1; x0 = x1;
end

% set output
if ii > 0
    xNout = interp1(xT(1:ii+1), xN(:,1:ii+1)', xtime, 'previous', 'extrap')';
else
    xNout = repmat(X0, 1, length(xtime));
end
if ii == Nrand
    xNout(:, xtime>t1) = NaN;
end

% check status
if ii == Nrand && t1<xtime(end)
    istatus = -1;
end

end
