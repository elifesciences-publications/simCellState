function [inputRaw, simOptions] = inputRead(icase0, icase1, outDir)
%  [inputRaw, simOptions] = inputRead(icase0, iset1)
% This function select the set of inputs to be used.
%  icase0    IN: Slector of settings to be used
%  icase1   IN: internal slector of settings to be used
%  inputRaw OUT: Input parameter structure containing the following fields:
%              - Nstate: number of states
%              - inputT: parameters for state transition (NX3)
%               Each row contains the i and j indexes for the transitions,
%               corresponding respectively to the state before and after
%               the transition and the transition rate. Default [].
%              - inputS: parameters for splitting (MX3)
%                Default [].
%  simOptions  OUT: simulation option structure.
%              - Nrun, Nrs; number of simulations and number of repeated
%                loops of Nrun simulations
%              - Nrand: max number of time step (number of random)
%              - seed0: initial seed
%              - T0, X0, time: initial conditions and time
%              - outFolder, fileName, isave are the saving options
%              - simPar: if >0 parallel computation is used
%              - plot: structure containing plot options
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
% default parameters
Nstate = 0;
inputT = []; inputS = [];
c = clock;
icase = icase0;
if ~isempty(icase1) && ~iscell(icase1)
    icase = [icase0 '_' icase1];
end
fileName = [icase '_' dateTag(c)];
outFolder0 = fullfile(outDir, 'io', 'OUT');
isave = 1;
simPar = 0; imex = 1; % no parallel/use mex as default
% default simulation parameters
Nrun = 1e4; % number of runs
Nrs = 1; % number of repeated runs
NN = 5e4; % number of random numbers
SD0 = 2*round(cputime*1e3)+1; % odd number
T0 = 0;
% default plot options
plotOptions = struct('idist', 1, 'idistT', 1, 'indxtime', 'end', 'idistScaled', 0, 'idistRef', '', ... % flag and options for plotting the distributions
    'imean', 1, 'imeanT', 1, ... % flag for plotting the mean
    'iext', 1, ... % flag for plotting the extinction(1), survival (-1)
    'isave', 1, ... % flag for saving figures
    'title', sprintf('%s', icase), ... % title
    'irs', 'all'); % which repeated simulation plot
addStr = []; % additional info to be saved

% specific parameters definition
switch icase0
    case {'PopulationAsymmetryM'} % A->2*A/2*B/AB, B->death
        Nstate = 6; % number of states
        gamma = icase1{1}; lambda = icase1{2}; r = icase1{3}; Delta = icase1{6};
        L1 = lambda*r; L2 = lambda*(1-2*r); L3 = lambda*r;
        Lm1 = L1*(1 + Delta);
        Om1 = (L1*(Delta + 1))/Delta;
        Lm2 = L2*(1 + Delta);
        Om2 = (L2*(Delta + 1))/Delta;
        Lm3 = L3*(1 + Delta);
        Om3 = (L3*(Delta + 1))/Delta;        
        inputT = [gamma    2    3  ; % B -> 0
                  Lm1      1    4  
                  Lm2   1    5  
                  Lm3         1    6  ]; % *delta/(delta+lambda)
        inputS = [Om1 4  1  1; % metastate
                  Om2 5  1  2;
                  Om3 6  2  2]; 
        % initial condition and final time
        if icase1{4} == 1
            X0 = icase1{4}*[1 0 0 0 0 0]';
        else
            X0 = round(icase1{4}*[gamma/(lambda+gamma) lambda/(lambda+gamma) 0 0 0 0])';
        end
        NN = 10e4; % number of random numbers
        iniCondition = struct('X0', X0, 'T0', T0);
        ts = min([lambda gamma]);
        % xtime = (0:0.5:15)/ts;
        xtime = 0:1:10; % (0:0.5:15)/ts;
        Nrun = icase1{5}; simPar = 1;
        outFolder0 = fullfile(outFolder0, icase0);
        addStr.tscale = ts;
    case {'PopulationAsymmetryR'} % A->2*A/2*B/AB, B->death
        Nstate = 3; % number of states
        gamma = icase1{1}; lambda = icase1{2}; r = icase1{3};
        inputT = [gamma    2    3  ]; % B -> 0
        inputS = [lambda*r      1  1  1; % A -> A+A
                 lambda*(1-2*r) 1  1  2;
                lambda*r        1  2  2]; 
        % initial condition and final time
        if icase1{4} == 1
            X0 = icase1{4}*[1 0 0]';
        else
            X0 = round(icase1{4}*[gamma/(lambda+gamma) lambda/(lambda+gamma) 0])';
        end
        NN = 10e4; % number of random numbers
        iniCondition = struct('X0', X0, 'T0', T0);
        ts = min([lambda gamma]);
        % xtime = (0:0.5:15)/ts;
        % outFolder0 = fullfile(outFolder0, icase0);
        addStr.tscale = ts;
    case {'PopulationAsymmetry'} % A->2*A, A->B, B->death
        Nstate = 3; % number of states
        gamma = icase1{1}; lambda = icase1{2}; omega = icase1{3};
        inputT = [gamma    2    3  ; % B -> 0
            omega    1    2  ]; % A -> B
        inputS = [lambda    1  1  1]; % A -> A+A
        % initial condition and final time
        if icase1{4} == 1
            X0 = icase1{4}*[1 0 0]';
        else
            X0 = round(icase1{4}*[gamma/(lambda+gamma) lambda/(lambda+gamma) 0])';
        end
        NN = 10e4; % number of random numbers
        iniCondition = struct('X0', X0, 'T0', T0);
        ts = min([lambda omega gamma]);
        xtime = (0:0.5:15)/ts;
        Nrun = icase1{5}; simPar = 1;
        outFolder0 = fullfile(outFolder0, icase0);
        addStr.tscale = ts;
    case {'AsymmetricDivision'} % A->A+B, B->death
        Nstate = 3; % number of states
        gamma = icase1{1}; lambda = icase1{2};
        inputT = [gamma    2    3  ]; % B -> 0
        inputS = [lambda    1  1  2]; % A -> A+B
        % initial condition and final time
        if icase1{3} == 1
            X0 = icase1{3}*[1 0 0]';
        else
            X0 = round(icase1{3}*[gamma/(lambda+gamma) lambda/(lambda+gamma) 0])';
        end
        iniCondition = struct('X0', X0, 'T0', T0);
        ts = min([lambda gamma]);
        xtime = 0:1:10; % (0:0.5:15)/ts;
        Nrun = icase1{4}; simPar = 1;
        outFolder0 = fullfile(outFolder0, icase0);
        addStr.tscale = ts;
    case {'genericNetwork'}
        CdCtInput = load(fullfile(outDir, 'io', 'IN', 'GENERIC', icase1{1}, icase1{1}));
        Cd = CdCtInput.xCdCt{1, icase1{2}};
        Ct = CdCtInput.xCdCt{2, icase1{2}};
        if isempty(Cd)
            inputRaw = []; simOptions = [];
            return
        end
        [inputT, inputS, Nstate] = CdCt2TL(Cd, Ct);
        X0 = zeros(Nstate, 1);
        ind0 = ceil(icase1{3}*CdCtInput.NstateR(icase1{2}));
        X0(ind0) = 1; % random initialization in the critical SCC
        iniCondition = struct('X0', X0, 'T0', T0);
        Nrun = icase1{4}; NN = 10^5;
        fileName = ['R' num2str(icase1{2}) '_I' num2str(ind0) '_' dateTag(c)];
        outFolder0 = fullfile(outFolder0, 'GENERIC', icase1{1});
        xtime = [0 1 2:0.25:50]/min([inputT(:,1); inputS(:,1)]);% [0:1:4 5:0.1:15]/min([inputT(:,1); inputS(:,1)]);
        Nrs = 1; simPar = 1;
    case {'stochasticModelFeedback', 'stochasticModelFeedbackParam'}
        TLInput = load(fullfile(outDir, 'io', 'IN', 'StoModel', icase1{1}, ['feedback_' num2str(icase1{2}) '_', icase1{3}]));
        TLInput = TLInput.(icase1{4});
        inputT = symb2input(TLInput{1});
        inputS = symb2input(TLInput{2});
        Nstate = TLInput{3};
        fbTag = [icase1{4} '_'];
        if length(icase1) >= 5 % reinitialization of the simulation
            if strcmp(icase0, 'stochasticModelFeedback')
                X0 = icase1{5};
                T0 = icase1{6}; TF = T0+250;
                iniCondition = struct('X0v', X0, 'T0', T0);
            else % perturbation in the initial condition
                X0 = zeros(Nstate, 1); X0(:) = 100*(1+icase1{5});
                T0 = 0; TF = T0+150;
                iniCondition = struct('X0', X0, 'T0', T0);
            end
        else
            X0 = zeros(Nstate, 1); X0(:) = 100;
            T0 = 0; TF = T0+150;
            iniCondition = struct('X0', X0, 'T0', T0);
        end
        xtime = T0:1:TF;
        plotOptions.indxtime = [2 length(xtime)];
        plotOptions.idistRef = '';
        plotOptions.isave = 1;
        % Nrun = 1e2; Nrs = 1;
        Nrun = 1e4; Nrs = 1;
        outFolder0 = fullfile(outFolder0, 'StoModelFeedback', icase1{1}, ['R' num2str(icase1{2})]);
        fileName = ['R' num2str(icase1{2}) '_' fbTag dateTag(c)];
        simPar = 1;
        plotOptions.iext = 0;
    case {'test_Alcolea', 'test_AlcoleaDoupe'} % A, B or death
        Nstate = 4; % number of states
        switch icase1{1}
            case {'control', 'test_IRIDIS'}
                lambda = 1.9/7; r = 0.1; gamma = 3.5/7; sigma = 1/7; delta = 0;
            case 'fit' % fitting
                lambda = icase1{2}(1)/7; r = icase1{2}(2); gamma = icase1{2}(3)/7; sigma = icase1{2}(4)/7; delta = icase1{2}(5);
                isave = 0;
            case 'fitOpt' % optimum
                lambda = icase1{2}(1)/7; r = icase1{2}(2); gamma = icase1{2}(3)/7; sigma = icase1{2}(4)/7; delta = icase1{2}(5);
                isave = 1;
        end
        %         rate      i -> j
        inputT = [gamma    2    3   ; % D -> S
            sigma    3    4   ; % S -> death
            ];
        %         rate      i -> j1  j2 % splitting
        inputS = [lambda*r*(1+delta)    1  1  1; % P -> P+P
            lambda*(1-2*r)        1  1  2; % P -> P+D
            lambda*r*(1-delta)    1  2  2;]; % P -> D+D
        % initial condition and final time
        X0 = zeros(Nstate, 1); X0(1) = 1; % start from 1 P (if starting from D --> clone extinction)
        iniCondition = struct('X0', X0, 'T0', T0);
        outFolder0 = fullfile(outFolder0, icase0);
        switch icase0
            case 'test_Alcolea'
                xtime = 0:0.1:10;
            case 'test_AlcoleaDoupe'
                xtime = [0:0.1:10 15 30 45 60 75 90];
        end
        plotOptions.indxtime = [find(xtime == 3) find(xtime == 7) find(xtime == 10)];
        plotOptions.idistScaled = 1; plotOptions.idistRef = '';
        
        switch icase1{1}
            case 'test_IRIDIS'
                simPar = 1;
                SD0 = 206563;
                Nrun = 1e5;
            otherwise
                simPar = 1;
                Nrun = 10e3;
        end
    case {'AA_feedback'} % A, B or death
        Nstate = 2; % number of states
        switch icase1
            case 'all'
                lambda0 = 0; gamma0 = 0;
                Klambda = 10; Kgamma = 0.1;
            case 'test'
                lambda0 = 1; gamma0 = 0;
                Klambda = 0; Kgamma = 0.1;
            case 'test1'
                lambda0 = 0; gamma0 = 1;
                Klambda = 300; Kgamma = 0;
        end
        inputT = @(n) [gamma0 + Kgamma*n     1    2   ; % A -> 0
            ];
        inputS = @(n) [lambda0 + Klambda/n    1  1  1]; % A -> A+A
        % initial condition and final time
        X0 = zeros(Nstate, 1); X0(1) = 250;
        iniCondition = struct('X0', X0, 'T0', T0);
        outFolder0 = fullfile(outFolder0, icase0);
        xtime = 0:0.1:10;
        plotOptions.indxtime = length(xtime);
        plotOptions.idistScaled = 1; plotOptions.idistRef = 'exp';
        Nrun = 10e4;
    case {'AB_feedback'} % A, B or death
        Nstate = 3; % number of states
        switch icase1
            case 'all'
                omegaA0 = 0; lambda0 = 0; omegaB0 = 0; gamma0 = 0;
                KomegaA = 0.1; Klambda = 10; KomegaB = 0.1; Kgamma = 0.1;
            case 'testLambda'
                omegaA0 = 1; lambda0 = 0; omegaB0 = 1; gamma0 = 1;
                KomegaA = 0; Klambda = 10; KomegaB = 0; Kgamma = 0;
            case 'testLambdaLinear'
                omegaA0 = 1; lambda0 = 1.5; omegaB0 = 1; gamma0 = 1;
                KomegaA = 0; Klambda = 0.01; KomegaB = 0; Kgamma = 0;
        end
        inputT = @(n) [omegaB0 + KomegaB*n    1    2   ; % A -> B
            omegaA0 + KomegaA*n    2    1   ; % B -> A
            gamma0 + Kgamma*n     2    3   ; % B -> 0
            ];
        % inputS = @(n) [lambda0 + Klambda/n    1  1  2]; % A -> A+B
        inputS = @(n) [lambda0 - Klambda*n    1  1  2]; % A -> A+B
        % initial condition and final time
        X0 = zeros(Nstate, 1); X0(1) = 60;
        iniCondition = struct('X0', X0, 'T0', T0);
        outFolder0 = fullfile(outFolder0, icase0);
        xtime = 0:0.1:20;
        plotOptions.indxtime = length(xtime);
        plotOptions.idistScaled = 1; plotOptions.idistRef = '';
        Nrun = 10e4;
    case {'stochasticModel', 'stochasticModelLargeN', 'stochasticModelHighRes', 'stochasticModelHighRes150'}
        Rfld = []; stabTag = [];
        if length(icase1) >= 5 % TL
            CdCtInput = load(fullfile(outDir, 'io', 'IN', 'StoModel', icase1{1}, icase1{4}));
            inputT = CdCtInput.(icase1{5}){1}; inputT = symb2input(inputT);
            inputS = CdCtInput.(icase1{5}){2}; inputS = symb2input(inputS);
            Nstate = CdCtInput.(icase1{5}){3};
            tg = [icase1{5} '_']; Rfld = ['R' num2str(icase1{2}) '_' icase1{3}];
        elseif length(icase1) == 3 % original models
            CdCtInput = load(fullfile(outDir, 'io', 'IN', 'StoModel', icase1{1}, icase1{1}));
            Cd = CdCtInput.xCdCt{1, icase1{2}};
            Ct = CdCtInput.xCdCt{2, icase1{2}};
            if isempty(Cd)
                inputRaw = []; simOptions = [];
                return
            end
            [inputT, inputS, Nstate] = CdCt2TL(Cd, Ct);
        elseif length(icase1) == 4 % stabilized models
            CdCtInput = load(fullfile(outDir, 'io', 'IN', 'StoModel', icase1{1}, icase1{4}));
            stabTag = 's_';
            Cd = CdCtInput.xCdCt{1, icase1{2}};
            Ct = CdCtInput.xCdCt{2, icase1{2}};
            if isempty(Cd)
                inputRaw = []; simOptions = [];
                return
            end
            [inputT, inputS, Nstate] = CdCt2TL(Cd, Ct);
        end
        X0 = zeros(Nstate, 1);
        switch icase0
            case {'stochasticModelHighRes150'}
                if icase1{3} <0 % random initial condition
                    sysProp = sysProperties('TL', inputT, inputS, Nstate);
                    Ndev = length(sysProp.indxD);
                    rng(SD0-2)
                    if Ndev > 0
                        icase1{3} = sysProp.indxD(randi(Ndev,1));
                    else
                        Nren = length(sysProp.indxR);
                        icase1{3} = sysProp.indxR(randi(Nren,1));
                    end
                end
                X0(icase1{3}) = 1;
                fileName = ['R' num2str(icase1{2}) '_I' num2str(icase1{3}) '_' stabTag dateTag(c)];
                % dt = 1; TendStab = 150; TendUnstab = 10;
                dt = 0.5; TendStab = 250; TendUnstab = 10;
                Nrun = 10e3; NN = 15e4;
                outFolder0 = fullfile(outFolder0, 'StochasticModelHighRes150', icase1{1});
            case {'stochasticModel', 'stochasticModelHighRes'}
                X0(icase1{3}) = 1;
                fileName = ['R' num2str(icase1{2}) '_I' num2str(icase1{3}) '_' stabTag dateTag(c)];
                dt = 1; TendStab = 50; TendUnstab = 10;
                if strcmp(icase0, 'stochasticModelHighRes')
                    Nrun = 10e4;
                    outFolder0 = fullfile(outFolder0, 'StochasticModelHighRes', icase1{1});
                else
                    Nrun = 1e4;
                    outFolder0 = fullfile(outFolder0, 'StoModel', icase1{1});
                end
            case 'stochasticModelLargeN'
                X0(:) = 100;
                fileName = ['R' num2str(icase1{2}) '_' stabTag dateTag(c)];
                dt = 0.5; TendStab = 50; TendUnstab = 20;
                Nrun = 1e4;
                outFolder0 = fullfile(outFolder0, 'StoModelLargeN', icase1{1}, Rfld);
        end
        if length(icase1) > 5 % reinitialization of the simulation
            X0 = icase1{6};
            T0 = icase1{7}; TendStab = 100; TendUnstab = 100;
            iniCondition = struct('X0v', X0, 'T0', T0);
        else
            iniCondition = struct('X0', X0, 'T0', T0);
        end
        if isnumeric(inputT) && isnumeric(inputT)
            sysProp = sysProperties('TL', inputT, inputS, Nstate);
        else
            sysProp.stability.flag = 1; % feedback
        end
        if sysProp.stability.flag % stable
            xtime = T0:dt:TendStab;
            plotOptions.indxtime = [1 10 25 50]+1;
        else % unstable
            xtime = T0:dt:TendUnstab;
            plotOptions.indxtime = [1 10]+1;
        end
        plotOptions.indxtime = [2 length(xtime)];
        plotOptions.idistRef = 'geo';
        plotOptions.isave = 0;
        Nrs = 1;
        simPar = 1;
    case {'deterministicModel', 'deterministicModelLargeN', ...
            'deterministicModelHighMean', ...
            'deterministicModelHighMean01', 'deterministicModelHighMean02', 'deterministicModelHighMean03', 'deterministicModelHighMean04', ...
            'deterministicModelHighMean05' 'deterministicModelHighMean06'}
        Rfld = [];
        if length(icase1) >= 5
            CdCtInput = load(fullfile(outDir, 'io', 'IN', 'DetModel', icase1{1}, icase1{4}));
            inputT = CdCtInput.(icase1{5}){1};
            inputS = CdCtInput.(icase1{5}){2};
            Nstate = CdCtInput.(icase1{5}){3};
            tg = [icase1{5} '_']; Rfld = ['R' num2str(icase1{2}) '_' icase1{3}];
        elseif length(icase1) == 4 && isstr(icase1{4})
            CdCtInput = load(fullfile(outDir, 'io', 'IN', 'DetModel', icase1{1}, icase1{4}));
            tg = 'r_';
            Cd = CdCtInput.xCdCt{1, icase1{2}};
            Ct = CdCtInput.xCdCt{2, icase1{2}};
            [inputT, inputS, Nstate] = CdCt2TL(Cd, Ct);
        else
            CdCtInput = load(fullfile(outDir, 'io', 'IN', 'DetModel', icase1{1}, icase1{1}));
            tg = ''; outFolder = 'DetModel';
            Cd = CdCtInput.xCdCt{1, icase1{2}};
            Ct = CdCtInput.xCdCt{2, icase1{2}};
            [inputT, inputS, Nstate] = CdCt2TL(Cd, Ct);
        end
        if strcmp(icase0(1:26), 'deterministicModelHighMean') % modify rates and overwrite settings
            sysProp0 = sysProperties('CdCt', Cd, Ct);
            indx0 = sysProp0.indx0;
            if sysProp0.stability.flag == 0 || isempty(sysProp0.indxC)
                inputRaw = []; simOptions = [];
                return
            end
            if strcmp(icase0(27:end), '05') || strcmp(icase0(27:end), '06')
                SD1 = icase1{4}; rng(SD1);
                xrr = rand(Nrun, Nstate^2 + Nstate + 1); % all possible death rates + division rates + initial condition
                xrr = xrr(icase1{2},:); % select only the applicable values for the simulation
                % reduce death rate
                randw = 1./((reshape(xrr(1:Nstate^2),Nstate, Nstate))*(10-2)+2); % between 1/2 and 1/10
                Ct(indx0,:) = Ct(indx0,:).*randw(indx0,:);
                % increase division rate
                randl = repmat(xrr(Nstate^2+1:Nstate^2+Nstate), Nstate,1)*(10-2)+2; % between 2 and 10
                Cd = Cd.*randl;
            else
                for ii = 1:length(indx0)
                    Ct(indx0(ii), Ct(indx0(ii),:) > 0) = 1/4; % reduce death rate
                end
                Cd = 15*Cd; % increase division rate
            end
            [inputT, inputS, Nstate] = CdCt2TL(Cd, Ct);
            tg = 's_'; % to be consistent with stochastic
            outFolder = 'DeterministicModelHighMean';
            if length(icase0) > 26 % modifications
                outFolder = [outFolder icase0(27:end)];
            end
        end
        X0 = zeros(Nstate, 1);
        switch icase0
            case {'deterministicModel', 'deterministicModelHighMean', 'deterministicModelHighMean01', ...
                    'deterministicModelHighMean02', 'deterministicModelHighMean03', 'deterministicModelHighMean04'} % clonal dyn
                X0(icase1{3}) = 1;
                fileName = ['R' num2str(icase1{2}) '_I' num2str(icase1{3}) '_' tg dateTag(c)];
                dt = 1; TendStab = 50; TendUnstab = 10;
                Nrun = 1e3;
                outFolder0 = fullfile(outFolder0, outFolder, icase1{1});
            case 'deterministicModelHighMean05'
                indxIC = [sysProp0.indxD, sysProp0.indxR, sysProp0.indxO]; % pick up random between develop, renewing and other
                ii = ceil(xrr(end)*length(indxIC)); ii(ii==0) = 1;
                X0(indxIC(ii)) = 1;
                fileName = ['R' num2str(icase1{2}) '_I' num2str(indxIC(ii)) '_' tg dateTag(c)];
                dt = []; xtime = [0:5:29 30:0.25:50];
                Nrun = 1e3;
                outFolder0 = fullfile(outFolder0, outFolder, icase1{1});
            case 'deterministicModelHighMean06'
                dt = []; xtime = [0:1:4.5 5:0.05:10]/min([inputT(:,1); inputS(:,1)]);
                indxActiveState = find(sum(sysProp0.A,1)~=0);
                nmax = zeros(1,Nstate);
                [~, J] = TL2AJ(inputT, inputS, Nstate);
                for ii = 1:Nstate
                    X0i = X0; X0i(ii) = 1;
                    if any(ii-indxActiveState== 0)
                        [~,y] = ode45(@(t,x) J*x, xtime, X0i);
                        nmax(ii) = sum(y(end,indxActiveState));
                    end
                end
                [~,indxIC] = max(nmax);
                X0(indxIC) = 1;
                fileName = ['R' num2str(icase1{2}) '_I' num2str(indxIC) '_' tg dateTag(c)];
                Nrun = 1e3;
                outFolder0 = fullfile(outFolder0, outFolder, icase1{1});
            case 'deterministicModelLargeN' % tissue dyn
                X0(:) = 100;
                fileName = ['R' num2str(icase1{2}) '_' tg dateTag(c)];
                dt = 0.5; TendStab = 50; TendUnstab = 20;
                Nrun = 1e4;
                outFolder0 = fullfile(outFolder0, 'DetModelLargeN', icase1{1}, Rfld);
        end
        if length(icase1) > 5 % reinitialization of the simulation
            X0 = icase1{6};
            T0 = icase1{7}; TendStab = 100; TendUnstab = 100;
            iniCondition = struct('X0v', X0, 'T0', T0);
        else
            iniCondition = struct('X0', X0, 'T0', T0);
        end
        sysProp = sysProperties('TL', inputT, inputS, Nstate);
        if sysProp.stability.flag % stable
            if ~isempty(dt)
                xtime = T0:dt:(TendStab+T0);
            end
            plotOptions.indxtime = [1 10 25 50]+1;
        else % unstable
            xtime = T0:dt:(TendUnstab+T0);
            plotOptions.indxtime = [1 10]+1;
        end
        plotOptions.idistRef = 'poi';
        plotOptions.isave = 1;
        Nrs = 1;
        simPar = 1;
    case {'MATH6140_AB'} % A, B or death
        Nstate = 3; % number of states
        switch icase1
            case 'hom_1'
                omegaA = 1; lambda = 1; omegaB = 1; gamma = 1;
            case 'hom_13'
                omegaA = 1/3; lambda = 1; omegaB = 1; gamma = 1/3;
            case 'hom_1neq'
                omegaA = 1; lambda = 1; omegaB = 1; gamma = 1;
        end
        %         rate      i -> j
        inputT = [omegaB    1    2   ; % A -> B
            omegaA    2    1   ; % B -> A
            gamma     2    3   ; % B -> 0
            ];
        %         rate      i -> j1  j2 % splitting
        inputS = [lambda    1  1  2]; % A -> A+B
        % initial condition and final time
        rho = omegaA/(omegaA+omegaB);
        Omega = omegaB/rho*(1+omegaB/lambda)/(1+omegaB/(lambda*rho))^2;
        X0T = 1;
        switch icase1
            case 'hom_1neq'
                pb = (1-rho)*0.1; pa = 1-pb;
            otherwise
                pb = (1-rho); pa = 1-pb;
        end
        p = [pa pb 0];
        iniCondition = struct('X0T', X0T, 'pX', p, 'T0', T0);
        outFolder0 = fullfile(outFolder0, icase0);
        xtime = unique([0 0.1 1:1:max(25, 10/Omega) 1/Omega 10/Omega]);
        plotOptions.indxtime = [2 find(xtime == 1/Omega) find(xtime == 10/Omega)];
        plotOptions.idistScaled = 1; plotOptions.idistRef = 'geo';
        Nrun = 10e4;
    case {'MATH6140_A'} % A or death
        Nstate = 2; % number of states
        switch icase1
            case 'crit'
                omega = 1; lambda = 1;
            case 'sub'
                omega = 1; lambda = 0.6;
            case 'sup'
                omega = 0.8; lambda = 1;
        end
        %         rate      i -> j
        inputT = [omega    1    2   ; % A -> 0
            ];
        %         rate      i -> j1  j2 % splitting
        inputS = [lambda    1  1  1]; % A -> A+A
        % initial condition and final time
        X0 = [1 0]'; iniCondition = struct('X0', X0, 'T0', T0);
        xtime = [0:1:25 50];
        plotOptions.indxtime = [length(xtime)-1 length(xtime)];
        plotOptions.idistT = 0; % no distribution plot
        plotOptions.imeanT = 0; % mean of total
        outFolder0 = fullfile(outFolder0, icase0);
        Nrun = 10e4;
    case 'test_expo'
        Nstate = 7; % number of states
        %         rate      i -> j
        inputT = [  1      2    1;
            1      1    2;
            1      2    3;
            1      2    4;
            1      3    4;
            1      5    2;
            1      6    2;
            1      1    7];
        %         rate      i -> j1  j2 % splitting
        inputS = [  1      4    1   2;
            1      3    1   2;
            1      1    5   6;];
        % initial condition and final time
        switch icase1
            case 'dA1B0'
                X0 = [1 0 0 0 0 0 0]';
            otherwise
                X0 = [15 5 0 1 0 0 0]';
        end
        xtime = 0:1:15;
        plotOptions.indxtime = [2 length(xtime)];
        iniCondition = struct('X0', X0, 'T0', T0);
        Nrs = 1; % Nrun = 1e2;
        % plotOptions.idistScaled = 1;
        plotOptions.idistRef = 'geo';
    case {'test_3s_mts0' 'test_3s_mts3'}
        Nstate = 3; % number of states
        gamma = 1; omegaa = 1; lambda = 1; omegab = 1; % default parameters
        %         rate      i -> j
        inputT = [gamma    2    3;
            omegab   1    2;
            omegaa   2    1;
            ];
        %         rate      i -> j1  j2 % splitting
        inputS = [lambda    1    1   1;
            lambda    1    2   2;
            lambda    1    1   2];
        % initial condition and final time
        switch icase1
            case 'dA1B0'
                X0 = [1 0 0]';
                iniCondition = struct('X0', X0, 'T0', T0);
        end
        % add metastate
        switch icase0
            case 'test_3s_mts3'
                Nstate = Nstate + 3;
                X0 = [X0; 0; 0; 0];
                iniCondition = struct('X0', X0, 'T0', T0);
                inputT = [inputT;
                    inputS(1,1) 1   4;
                    inputS(2,1) 1   5;
                    inputS(3,1) 1   6;
                    ];
                inputS(:,1) = inputS(:,1)*100;
                inputS(:,2) = [4 5 6]';
        end
        xtime = 0:1:15;
        plotOptions.indxtime = [2 length(xtime)-1];
    case {'test_3s_hom_1' 'test_3s_hom_08' 'test_3s_hom_12' 'test_3s_hom_03'}
        Nstate = 3; % number of states
        gamma = 1; omegaa = 1; lambda = 1; omegab = 1; % default parameters
        switch icase0
            case 'test_3s_hom_03'
                gamma = omegab/3; omegaa = omegab/3;
            case 'test_3s_hom_08'
                gamma = 0.8; omegaa = 0.8;
            case 'test_3s_hom_12'
                gamma = 1.2; omegaa = 1.2;
        end
        %         rate      i -> j
        inputT = [omegaa    2    1   ;
            omegab    1    2   ;
            gamma     2    3   ;];
        %         rate      i -> j1  j2 % splitting
        inputS = [lambda    1    1   2];
        % time
        ro = omegaa/(omegaa+omegab);
        Omega = omegab/ro*(1+omegab/lambda)/(1+omegab/(lambda*ro))^2;
        xtime = unique([0:1:10 20 1/Omega 10/Omega]);
        plotOptions.indxtime = [find(xtime == 1/Omega) find(xtime == 10/Omega)];
        plotOptions.idistScaled = 1; plotOptions.idistRef = 'exp';
        % initial condition distribution
        switch icase1
            case 'dA1B0'
                X0 = [1 0 0]';
                iniCondition = struct('X0', X0, 'T0', T0);
            case 'pA1B0'
                X0T = 1;
                p = [1 0 0];
                iniCondition = struct('X0T', X0T, 'pX', p, 'T0', T0);
            case 'dA0B1'
                X0 = [0 1 0]';
                iniCondition = struct('X0', X0, 'T0', T0);
            case 'pA0B1'
                X0T = 1;
                p = [0 1 0];
                iniCondition = struct('X0T', X0T, 'pX', p, 'T0', T0);
            case 'dA1B1'
                X0 = [1 1 0]';
                iniCondition = struct('X0', X0, 'T0', T0);
            case 'pA1B1'
                X0T = 2;
                p = [0.5 0.5 0];
                iniCondition = struct('X0T', X0T, 'pX', p, 'T0', T0);
            case 'pABeq'
                X0T = 1;
                pb = gamma/lambda/(1+gamma/lambda); pa = 1-pb;
                p = [pa pb 0];
                iniCondition = struct('X0T', X0T, 'pX', p, 'T0', T0);
        end
    case {'SD0_3098563_20171127121322_s01', 'SD0_3098563_20171127121322_s02', 'SD0_3098563_20171127121322_s03', 'SD0_3098563_20171127121322_s04', ...
            'SD0_3098563_20171127121322_s11', 'SD0_3098563_20171127121322_s12', 'SD0_3098563_20171127121322_s13', 'SD0_3098563_20171127121322_s14', ...
            } % check statistics
        if str2double(icase0(end-1)) > 0
            Nstate = 5;
        else
            Nstate = 3;
        end
        AB = load(fullfile('..', 'rfKMC', 'sdynHet1', 'testCase_AB', icase0(1:end-4)));
        Nrun = AB.Nrun; Nrs = AB.Nrs;
        NN = AB.NN;
        inputT = [AB.omegaa    2    1   ; % B->A
            AB.omegab    1    2   ; % A->B
            AB.gamma     2    3   ]; % B->0
        inputS = [AB.lambda    1    1   2]; % A->A+B
        SD0 = AB.SD0;
        X0 = zeros(Nstate,1);
        X0(1:2) = AB.X0; T0 = AB.T0;
        xtime = AB.xtime;
        isave = 0;
        % variations
        if str2double(icase0(end)) > 1
            SD0 = 2*round(cputime)+1; % different inital seed
        end
        if str2double(icase0(end)) > 2 % flip inputs
            inputT = flip(inputT);
        end
        if str2double(icase0(end)) > 3 % reshuffle states
            if str2double(icase0(end-1)) > 0
                indx = [4 2 5 3 1];
            else
                indx = [2 3 1];
            end
            inputS1 = inputS; inputT1 = inputT; X01 = X0;
            for istate = 1:Nstate
                [is,js] = find(inputS1(:,2:end) == istate);
                inputS(sub2ind([size(inputS1,1) 4], is, js+1)) = indx(istate);
                [is,js] = find(inputT1(:,2:end) == istate);
                inputT(sub2ind([size(inputT1,1) 3], is, js+1)) = indx(istate);
                X0(indx(istate)) = X01(istate);
            end
        end
        iniCondition = struct('X0', X0, 'T0', T0);
    case {'SD0_3098563_20171127121322'}
        Nstate = 3;
        AB = load(fullfile('..', 'rfKMC', 'sdynHet1', 'testCase_AB', icase0));
        Nrun = AB.Nrun; Nrs = AB.Nrs;
        NN = AB.NN;
        inputT = [AB.omegaa    2    1   ; % B->A
            AB.omegab    1    2   ; % A->B
            AB.gamma     2    3   ]; % B->0
        inputS = [AB.lambda    1    1   2]; % A->A+B
        SD0 = AB.SD0;
        X0 = [AB.X0; 0]; T0 = AB.T0;
        xtime = AB.xtime;
        isave = 0;
        iniCondition = struct('X0', X0, 'T0', T0);
    case {'testCase_AA_20171127113837', 'testCase_AA_20171127115412', 'testCase_AA_2017112712712' 'testCase_AA_201712271521'}
        Nstate = 2;
        Nrs = 1;
        AA = load(fullfile('..', 'rfKMC', 'rfkMC_1', 'testCase_AA', icase0));
        Nrun = AA.Nrun;
        NN = AA.NN;
        inputT = [AA.r2 1 2];
        inputS = [AA.r1 1 1 1];
        SD0 = AA.SD0;
        X0 = [AA.k0; 0];
        T0 = AA.t0;
        xtime = AA.tf;
        isave = 1;
        iniCondition = struct('X0', X0, 'T0', T0);
    otherwise
        % default empty options
        xtime = []; iniCondition = [];
end

% problem param
inputRaw = struct('Nstate', Nstate, 'inputT', inputT, 'inputS', inputS, 'addStr', addStr);
% plot option
if strcmp(plotOptions.indxtime, 'end')
    plotOptions.indxtime = length(xtime);
end
if strcmp(plotOptions.irs, 'all')
    plotOptions.irs = 1:Nrs;
end

% simulation options
simOptions = struct('Nrun', Nrun, 'Nrs', Nrs, 'Nrand', NN, ...
    'seed0', SD0, ...
    'iniCondition', iniCondition, 'time', xtime, ...
    'outFolder', fullfile(outFolder0, fileName), 'fileName', fileName, 'isave', isave, ...
    'simPar', simPar, 'mex', imex, 'plot', plotOptions);

end