clear all
% close all

% This script process all the simulations. The following cell array is saved:
% - xout = {isim, simModel, simOption, simOutPostpro}
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
d0 = fullfile('..', '..', 'IO', 'OUT', 'GENERIC');
tg = '20190409150503_p30_V01';
maxExtRate = []; xTime = []; % initialization
outFolder = ['out_cons_' tg]; icons = 1; tgo = ''; xTime = [5 10 15 20 25];
% outFolder = ['out_ncons_' tg]; icons = 0; maxExtRate = [50 60 70 80 90 95 98 99]/100; tgo =  ''; % ['_' num2str(round(maxExtRate*100))];
NN = 1e3; % number of simulations

% initialization
xout = cell(NN,4);
% loop
for ii = 1:NN
    % load simulation
    isim = ii;
    simDir = dir(fullfile(d0, outFolder, ['R' num2str(isim) '_I*']));
    simCase = simDir.name;
    xOut = load(fullfile(d0, outFolder, simCase, simCase));
    out = ppro_simPaper(xOut.simModel, xOut.simOptions, xOut.simOut, maxExtRate, xTime);
    if exist(fullfile(d0, outFolder, simCase, [simCase '_1.mat']), 'file')
        xOut1 = load(fullfile(d0, outFolder, simCase, [simCase '_1']));
        out1 = ppro_simPaper(xOut1.simModel, xOut1.simOptions, xOut1.simOut, maxExtRate, xTime);
        out1.indxT = out1.indxT + out.inan-1;
        clear xOut1
        ff = setdiff(fieldnames(out), {'indxR', 'indxC', 'indxD', 'extRate', 'inan', 'nn_sim', 'indxT'});
        for jj = 1:length(ff) % remove nan from prev out
            out.(ff{jj})(out.inan+1:end) = [];
        end
%         indx1 = find(xOut.simOptions.time < xOut1.simOptions.time(1));
%         xOut.simOut{1}.xstate = cat(3, xOut.simOut{1}.xstate(:,:,indx1), xOut1.simOut{1}.xstate);
%         xOut.simOptions.time = [xOut.simOptions.time(indx1), xOut1.simOptions.time];
    else
        out1 = [];
    end
    if exist(fullfile(d0, outFolder, simCase, [simCase '_1_2.mat']), 'file')
        xOut2 = load(fullfile(d0, outFolder, simCase, [simCase '_1_2']));
        out2 = ppro_simPaper(xOut2.simModel, xOut2.simOptions, xOut2.simOut, maxExtRate, xTime);
        out2.indxT = out2.indxT + out.inan-1 + out1.inan-1;
        clear xOut2
        ff = setdiff(fieldnames(out), {'indxR', 'indxC', 'indxD', 'extRate', 'inan', 'nn_sim', 'indxT'});
        for jj = 1:length(ff) % remove nan from prev out
            out1.(ff{jj})(out1.inan+1:end) = [];
        end
%         indx1 = find(xOut.simOptions.time < xOut1.simOptions.time(1));
%         xOut.simOut{1}.xstate = cat(3, xOut.simOut{1}.xstate(:,:,indx1), xOut1.simOut{1}.xstate);
%         xOut.simOptions.time = [xOut.simOptions.time(indx1), xOut1.simOptions.time];
    else
        out2 = [];
    end
    if ~isempty(out1)
        ff = setdiff(fieldnames(out), {'indxR', 'indxC', 'indxD'});
        for jj = 1:length(ff)
            if any(ismember({'indxT', 'extRate', 'nn_sim'}, ff{jj}))
                out.(ff{jj})(end) = [];
            end
            out.(ff{jj}) = [out.(ff{jj}) out1.(ff{jj})(2:end)];
            if ~isempty(out2)
                if any(ismember({'indxT', 'extRate', 'nn_sim'}, ff{jj}))
                    out.(ff{jj})(end) = [];
                end
                out.(ff{jj}) = [out.(ff{jj}) out2.(ff{jj})(2:end)];
            end
        end        
    end
    
    % set output 
    xout{ii, 3} = xOut.simOptions;
    xout{ii, 2} = xOut.simModel;
    xout{ii, 1} = isim;
    if icons % add equivalent parameters
        [lambda1Eq, lambda2Eq, gammaEq, p, q, mu, ~, xsol] = findEquivalentParam(xOut.simModel, xOut.simOptions);
        out.lambda1Eq = lambda1Eq; out.lambda2Eq = lambda2Eq; out.gammaEq = gammaEq;
        out.p = p; out.q = q; out.mu = mu;
        out.xsol = xsol;
    end
    xout{ii, 4} = out;
end

% save
pproFolder = fullfile(d0, ['ppro_' outFolder]);
if ~exist(pproFolder, 'dir')
    mkdir(pproFolder)
end
save(fullfile(pproFolder, ['out_ppro' tgo]), 'xout')