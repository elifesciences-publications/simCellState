function simModel = inputSetup(inputRaw)
%  simModel = inputSetup(inputRaw)
% This function reorganize the inputs parameters into the simulation model.
% Both state transition (S1 -> S2) and cell splitting (S1 -> S2+S3) are
% considered.
% 
%  inputRaw IN: Input parameter structure. See inputRead function. 
%  simModel OUT: Structure containing the following fields:
%       sRate: specific rate of each possible event (N+M,1)
%       deltaState: state variation corresponding to each event (NS,N+M) 
%           Each j-column is consitent with the j-row in sRate. 
%       indxState: State before transition or splitting for each possible
%            event in sRate (N+M,1)
%            This is used to compute the rates of each event from the
%            specific rates. Ex. Si -> Sj with sr_ij, r_ij = sr_ij*nSi
%       indxActiveState: index of active states - not death (1,L)
%       Nopt: total number of possible options for each event.
%       Nstate: number of states
%       odeM:  matrix containing the parameters to integrate the mean
%       extEst: structure containing the parameter for the anlaytical
%            estimation of the extinction probability
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

% inputs
Nstate = inputRaw.Nstate;
inputT = inputRaw.inputT;
inputS = inputRaw.inputS;
[A, J, symFlag, Tf, Lf, inputT, inputS] = TL2AJ(inputT, inputS, Nstate);

% transition matrix, transition specific rates and delta state
NT = size(inputT, 1); % number of possible transitions
indxStateT = inputT(:,2);
sRateT = Tf; % specific transition rate
deltaStateT = zeros(Nstate, NT);
for ii = 1:2
    indx = sub2ind([Nstate NT], inputT(:,1+ii)', 1:1:NT);
    if ii == 1 
        deltaStateT(indx) = -1; % before transition
    else
        deltaStateT(indx) = 1; % after transition
    end
end

% splitting matrix
NS = size(inputS, 1); % number of possible splitting
if NS > 0
    indxStateS = inputS(:,2);
    sRateS = Lf;
    deltaStateS = zeros(Nstate, NS);
    for ii = 1:3
        indx = sub2ind([Nstate NS], inputS(:,1+ii)', 1:1:NS);
        if ii == 1
            deltaStateS(indx) = -1; % before splitting
        else
            deltaStateS(indx) = deltaStateS(indx) + 1; % after splitting
        end
    end
else
    sRateS = [];
    deltaStateS = [];
    indxStateS = [];
end

% Combination of splitting and transition
if symFlag
    indxActiveState = find(sum(A(NaN),1)~=0); % not death
    sRate = @(x) [sRateS(sum(x(indxActiveState))); sRateT(sum(x(indxActiveState)))];
else
    sRate = [sRateS; sRateT];
    indxActiveState = find(sum(A,1)~=0); % not death
end
deltaState = [deltaStateS deltaStateT];
indxState = [indxStateS; indxStateT];

% ADDITIONAL PARAMETERS
if symFlag
    odeM = @(t,x) J(sum(x(indxActiveState)))*x;
    extEst = [];
else
    % input matrix for ODE
    odeM = J;
    % reproduction generating function
    [p0, p1, p2, rtot] = repGenFcnCoeff(inputT, inputS, Nstate);
    % extinction parameters
    extEst.rtot = rtot';
    extEst.p0 = p0;
    extEst.p1 = p1;
    extEst.p2 = p2;
end

% OUTPUT STRUCTURE
simModel = struct('sRate', sRate, 'deltaState', deltaState, 'indxState', indxState, ...
    'Nstate', Nstate, 'indxActiveState', indxActiveState, ...
    'odeM', odeM, 'extEst', extEst, 'inputRaw', inputRaw);

end