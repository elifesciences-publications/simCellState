function [simOut, simOutFull] = run_simGAD1(p, q, a, b, c, GEq, xTime, Nrun)
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

% MODEL:
% X1, X2 critical SCC, X3 subcritical, X4 is death
% X1 -> X1 + X3, L1; X2 -> X2 + X3, L2
% X1 -> X2, O12; X2 -> X1, O21
% X3 -> X3 + X3, L3
% X3 -> X4, G
% PARAMETERS
% GEq = G; 
% a = L1/L2, b = O12/O21, 
% hatOmega = c = O12/GEq
% L1eq = p*GEq
% L2Eq = q*GEq
L2Eq = q*GEq;
L1Eq = p*GEq;
hatOmega = c;
% internal parameter
G = GEq;
L3 = L2Eq;
O12 = hatOmega*GEq; O21 = O12/b; 
% J = [-O12 O21; O12 -O21; L1 L2];
% -O12*n1+O21*n2 = 0
% n1+n2 = 1 --> n1 = 1-n2
% (O12*n1-O21*n2 = 0)
% -O12*(1-n2)+O21*n2 = 0, -O12+(O12+O21)*n2 = 0, 
% n2 = O12/(O12+O21);
% n1 = 1-n2 = O21/(O12+O21);
% L1eq = (L1*O21/(O12+O21)+L2*O12/(O12+O21))/1
% L2 = L1/a
% L1Eq = L1*O21/(O12+O21)+L1/a*O12/(O12+O21)
% L1Eq = L1*(O21+O12/a)/(O12+O21)
L1 = L1Eq*(O12+O21)/(O21+O12/a);
L2 = L1/a;
% % eqns = [p*G == L1*O21/(O12+O21)+L2*O12/(O12+O21), a == L1/L2, b == O12/O21, c == O12/L1];
% % vars = [L1 L2 O12 O21];
% % [solL1, solL2, solO12, solO21] = solve(eqns, vars)
% % solL1 = (a*(G*p + G*b*p))/(a + b); solL2 = (G*p + G*b*p)/(a + b); 
% % solO12 = (G*a*c*p + G*a*b*c*p)/(a + b); solO21 = (G*a*c*p + G*a*b*c*p)/(b^2 + a*b)
% % L1 = G*a*p*(b + 1)/(a + b); L2 = G*p*(b + 1)/(a + b);
% % O12 = G*a*c*p*(b + 1)/(a + b); O21 = G*a*c*p*(b + 1)/(b*(a + b));
% 
% % c == O12/G
% O12 = c*G; O21 = O12/b; 
% L1 = (a*(G*O12*p + G*O21*p))/(O12 + O21*a);
% L2 = (G*O12*p + G*O21*p)/(O12 + O21*a);

r = [L1 L2 L3 O12 O21 G]; r = r(r>0); r_min = min(r);

% setup simulation
Nstate = 4; X0 = zeros(Nstate,1); X0(1) = 1; T0 = 0;
inputRaw.inputT = [O12 1 2; O21 2 1; G 3 4]; inputRaw.inputS = [L1 1 1 3; L2 2 2 3; L3 3 3 3];
inputRaw.Nstate = Nstate;
simOptions = struct('Nrun', Nrun, 'Nrs', 1, 'Nrand', 10e4, ...
    'seed0', 2*round(cputime*1e3)+1, ...
    'iniCondition', struct('X0', X0, 'T0', T0), 'time', xTime/r_min, ...
    'isave', 0, 'simPar', 1, 'mex', 1, 'plot', []);
simModel = inputSetup(inputRaw);
simOutFull = simCellStateLoop(simModel, simOptions);

% postprocessing
simOut = ppro_simPaper(simModel, simOptions, simOutFull, 0);
[lambda1Eq, lambda2Eq, gammaEq, p, q] = findEquivalentParam(simModel, simOptions);

% add additional parameters
outParam = struct('L1', L1, 'L2', L2, 'L3', L3, 'O12', O12, 'O21', O21, 'G', G, ...
    'p1', L1/G, 'p2', L2/G, ...
    'mu_tg', p/(1-q), 'mu', simOut.mean_n2_sim(simOut.indxT), 'mu1', L1/(G-L3), 'mu2', L2/(G-L3), ...
    'a', L1/L2, 'b', O12/O21, 'c', O12/G, 'lambda1Eq', lambda1Eq, 'lambda2Eq', lambda2Eq, 'gammaEq', gammaEq, ...
    'p', p, 'q', q);
outParam.mix_tg = (outParam.mu_tg-outParam.mu2)/(outParam.mu1-outParam.mu2); % target
outParam.mix = (outParam.mu-outParam.mu2)/(outParam.mu1-outParam.mu2); % final
outParam.nn = min(simOut.nn_sim{1}):1:max(simOut.nn_sim{1});
[outParam.g1nn, ss1] = gnnpdf(outParam.nn, outParam.mu1, outParam.p1, q); % use nn2 to shift by one
[outParam.g2nn, ss2] = gnnpdf(outParam.nn, outParam.mu2, outParam.p2, q);
outParam.gnn = outParam.mix*outParam.g1nn + (1-outParam.mix)*outParam.g2nn;
% rescaling of the bimodal distribution
bm_mu = outParam.mix*outParam.mu1 + (1-outParam.mix)*outParam.mu2;
outParam.xx = outParam.nn/bm_mu;
outParam.gxx = outParam.gnn*bm_mu;
bm_sigma = sqrt(outParam.mix*(ss1 + (outParam.mu1-bm_mu)^2) + (1-outParam.mix)*(ss2 + (outParam.mu2-bm_mu)^2)); % standard dev
outParam.xxtilde = (outParam.xx-1)*bm_mu/bm_sigma;
outParam.gxxtilde = outParam.gxx*bm_sigma/bm_mu;

simOut.outParam = outParam;

end


function [gnn, ss] = gnnpdf(nn, mu, p, q)

    if p > 15 % normal distr
        gnn = normpdf(nn, mu, mu*sqrt(1/p));
        ss = mu*sqrt(1/p);
    elseif q < 0.1
        gnn = poisspdf(nn, mu);
        ss = mu;
    else % full analytical distribution (cont but converted in nn?)
        gnn = NaN;
        ss = NaN;
    end
    
end