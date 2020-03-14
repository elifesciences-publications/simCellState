function inputRaw1 = findModifiedRates(p_tg, inputRaw, indxR1, Rmax, simOptions)

% This function search for the modified parameters (stochastic search)
% A Particle Swarm Optimization algorithm (PSO) is used.
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

Ytg = 1;
if Rmax < 2e3
    [~, indxS] = intersect(inputRaw.inputS(:,2), indxR1);
    lia = ismember(inputRaw.inputT(:,2), indxR1); indxT = find(lia);
    fun = @(x) updateModel(x, inputRaw, indxT, indxS, simOptions);
    opt = optimoptions('ga');
    [x,Ygbest] = ga(@(x, iter) abs(fun(x) - p_tg),length(indxS)+length(indxT),[],[],[],[],...
        0.1*ones(1, length(indxS)+length(indxT)),Rmax*ones(1, length(indxS)+length(indxT)),[],[],opt);
    if Ygbest > Ytg
        [p1, q1] = updateModel(x, inputRaw, indxT, indxS, simOptions);
        disp([p1 q1 p1/(1-q1)])
        inputRaw1 = findModifiedRates(p_tg, inputRaw, indxR1, Rmax*2, simOptions);
    else
        [p1, q1, inputRaw1] = updateModel(x, inputRaw, indxT, indxS, simOptions);
        disp([p1 q1 p1/(1-q1)])
    end  
else
    inputRaw1 = [];
end