clear all
close all

% This function build a set of SCC. The size of SCC (NSCC) is used to load
% the raw data file containing all the possible combination of transition
% and division matrices. 
% The output is a matrix of cells in which the SCCs are ordered as function
% of the number of division (rows) and number of transitions (cols).
% If, for a given number of divisions and transitions there are more than
% the NPMax parameter, then NNMax combinations are randomly selected.
% Parameter icons controls the conserved/not conserved networks:
% pre-filtering the pool of Cd,Ct if conserved and/or filtering after the
% generation of the network.
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
NSCC = 4;
icons = 0;
NPMax = 50e4; % max number of networks that can be processed
NNMax = 1e4; % number of cases randomly picked
dOut = 'SCC';
SD0 = 301; rng(SD0)

% load data
fileName = sprintf('rawDT_N%d.mat', NSCC);
load(fullfile(dOut,fileName))

% filter data if only conserved are wanted
if icons == 1
    Ct(:,any(Ct(NSCC+1:NSCC+1:end,:))) = []; % remove solutions with transitions outside the SCC
    Cd = Cd(:,any(Cd(NSCC+1:NSCC+1:end,:)==1)); % keep solutions having at least one asymmetric division
    Cd(:,any(Cd(NSCC+1:NSCC+1:end,:)==2)) = []; % remove solutions with double division going outside the SCC (equiv. to transition)
    % remove unstable ones
    unstFlag = zeros(1,size(Cd, 2));
    for ii = 1:size(Cd, 2)
        Cdii = reshape(Cd(:,ii), NSCC+1, NSCC+1);
        if any(sum(Cdii(1:NSCC,1:NSCC)) > 1)
            unstFlag(ii) = 1;
        end
    end
    Cd(:, unstFlag==1) = [];
    consTg = '_C';
else
    consTg = '_nC';
end

% number of combinations
sCd = sum(Cd)/2; sCt = sum(Ct);
ND = unique(sCd); NT = unique(sCt);

xOutF = cell(length(ND), length(NT)); xOut = cell(length(ND), length(NT));
xProp = struct('NSCC', [], 'ND', [], 'NT', [], 'Nall', [], 'Nallc', [], 'NN', [], 'indxType1', []);
dateTg = '';
tic
for it = 1:length(NT)
    for id = 1:length(ND)
        Cdid = Cd(:, sCd == ND(id));
        Ctit = Ct(:, sCt == NT(it));
        if size(Cdid, 2) > 0 && size(Ctit,2) > 0
            if size(Cdid, 2)*size(Ctit,2) <= NPMax % process all networks
                [xOutF{id,it}, aux]= processCdCtLoop1(Cdid, Ctit, icons);
                NN = size(xOutF{id,it}, 1);
                NNc = size(aux, 1); 
                if NN > 0
                    xtype1 = cellfun(@(x) x==1, xOutF{id,it}(:,4));
                    indxType1 = find(xtype1);
                else
                    indxType1 = [];
                end
            else
                indxCd = randi(size(Cdid, 2), [1, NNMax]);
                indxCt = randi(size(Ctit, 2), [1, NNMax]);
                [xOutr, aux] = processCdCtLoop2(Cdid(:, indxCd), Ctit(:, indxCt), icons);
                xOutF{id,it} = xOutr; 
                NN = size(xOutF{id,it}, 1); 
                NNc = size(aux, 1); 
                if NN > 0
                    xtype1 = cellfun(@(x) x==1, xOutF{id,it}(:,4));
                    indxType1 = find(xtype1); % add while loop!?! I need at least 10 in case of large SCC
                else
                    indxType1 = [];
                end
                dateTg = ['_' dateTag(clock)];
            end
        end
        xProp(id,it).NSCC = NSCC; 
        xProp(id,it).ND = ND(id); 
        xProp(id,it).NT = NT(it);
        xProp(id,it).Nall = size(Cdid, 2)*size(Ctit,2);
        xProp(id,it).Nallc = NNc;
        xProp(id,it).NN = NN;
        xProp(id,it).indxType1 = indxType1;
    end
end
toc

% save in file
outfileName = sprintf('SCC_N%d%s%s', NSCC, consTg, dateTg);
save(fullfile(dOut, outfileName), 'xOutF', 'xProp')

% for checking
% ------------
xOutF = xOutF(:);
xProp = xProp(:);

% Summary
totNall = sum(arrayfun(@(x) x.Nall,xProp));
totNallc = sum(arrayfun(@(x) x.Nallc,xProp));
totNN = sum(arrayfun(@(x) x.NN,xProp));
totNN1 = sum(arrayfun(@(x) length(x.indxType1),xProp));

disp([length(sCd), length(sCt), totNall, totNallc, totNN1, totNN])