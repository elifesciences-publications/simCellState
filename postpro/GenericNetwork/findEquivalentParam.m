function [lambda1Eq, lambda2Eq, gammaEq, p, q, mu2, JEq, xsol] = findEquivalentParam(simModel, simOptions)

% starting from a complex model, equivalent parameters of a simple AD models are found
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

% initial system
indxD = setdiff(1:1:simModel.Nstate, simModel.indxActiveState); 
inputS = simModel.inputRaw.inputS;
inputT = simModel.inputRaw.inputT;
Nstate = simModel.inputRaw.Nstate;
sysProp = sysProperties('TL', inputT, inputS, Nstate);
odeM = sysProp.J;
% X0 = zeros(1,Nstate); X0(indxR(1)) = 1;
X0 = simOptions.iniCondition.X0; 
tspan = [simOptions.time(1) simOptions.time(end)*3];

% first change the output of division into death into transition if necessary
ddD = find(all(inputS(:,3:4)==indxD, 2)); % 2 death states
if ~isempty(ddD) % change to transitions
    for jj = 1:length(ddD)
        [~, irow] = intersect(inputT(:, 2:3), inputS(ddD(jj), 2:3), 'row');
        if isempty(irow) % add a new transition
            inputT = [inputT; inputS(ddD(jj), 1:3)];
        else % transition already existing
            inputT(irow, 1) = inputS(ddD(jj), 1)+inputT(irow, 1);
        end
    end
    inputS(ddD,:) = [];
end
ddD = find(any(inputS(:,3:4)==indxD, 2)); % 1 in 1 death
if ~isempty(ddD) % change to transitions
    for jj = 1:length(ddD)
        if inputS(ddD(jj), 2) ~= setdiff(inputS(ddD(jj), 3:4), indxD) % transition not to itself
            [~, irow] = intersect(inputT(:, 2:3), [inputS(ddD(jj), 2) setdiff(inputS(ddD(jj), 3:4), indxD)], 'row');
            if isempty(irow) % add a new transition
                inputT = [inputT; [inputS(ddD(jj), 1:2) setdiff(inputS(ddD(jj), 3:4), indxD)]];
            else
                % transition already existing
                inputT(irow, 1) = inputS(ddD(jj), 1)+inputT(irow, 1);
            end
        end
    end
    inputS(ddD,:) = [];
end

% compute system property
sysProp = sysProperties('TL', inputT, inputS, Nstate);
odeMmod = sysProp.J;
% index renewing, committed
indxR = sysProp.sccIndx{cellfun(@(x) any(x == 1), sysProp.sccIndx) == 1};
% indxR = sysProp.sccIndx{sysProp.sccType == 0};
indxC = setdiff(1:1:Nstate, [indxR indxD]);
% find outgoing from renewing
[ii, ~] = find(odeMmod(indxC,indxR)); 
indxO = unique(indxC(ii));
% find incoming for death
indxI = find(odeMmod(indxD,:));

% estimate steady state mean: ode*x = 0; sum(n1) = 1;
xmeanmod = linsolve([odeMmod([indxR indxC],[indxR indxC]); ones(1, length(indxR)) zeros(1, length(indxC))],...
    [zeros(length([indxR indxC]),1); 1]);

% estimation of equivalent parameters
% lambda1
dndtR = odeMmod([indxR indxO], indxR)*xmeanmod(indxR);
lambda1Eq = sum(dndtR)/sum(xmeanmod(indxR));
% gamma
dndtD = odeMmod(indxD, indxI)*xmeanmod(indxI);
gammaEq = dndtD/sum(xmeanmod(indxC));
% lambda2
lambda2Eq = gammaEq - lambda1Eq./sum(xmeanmod(indxC));

%  mu, equivalent p and equivalent jacobian
mu2 = sum(xmeanmod(indxC));
p = lambda1Eq/gammaEq;
q = lambda2Eq/gammaEq;
JEq = [0 0 0; lambda1Eq lambda2Eq-gammaEq 0; 0 gammaEq 0];

if nargout == 8
    [time, xmean] = ode45(@(t,x) odeM*x, tspan, X0);
    n1 = sum(xmean(:,indxR),2); n2 = sum(xmean(:,indxC),2);
    [timemod, xmeanmod] = ode45(@(t,x) odeMmod*x, tspan, X0); 
    n1mod = sum(xmeanmod(:,indxR),2); n2mod = sum(xmeanmod(:,indxC),2);
    X0 = [1; 0; 0];
    [timeEq, xmeanEq] = ode45(@(t,x) JEq*x, tspan, X0); 
    n1Eq = xmeanEq(:,1); n2Eq = xmeanEq(:,2);
    xsol = struct('time', time, 'n1', n1, 'n2', n2, ...
        'timeMod', timemod, 'n1Mod', n1mod, 'n2Mod', n2mod, ...
        'timeEq', timeEq, 'n1Eq', n1Eq, 'n2Eq', n2Eq);
end

return