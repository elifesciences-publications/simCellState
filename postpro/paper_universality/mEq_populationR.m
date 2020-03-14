function dPvdt = mEq_populationR(~, Pv, Na, Nb, lambda, r, gamma)

% na -> na+na, lambda*r
% na -> na+nb, lambda*(1-2*r)
% na -> nb+nb, lambda*r
% nb -> death, gamma

% NA, NB matrix

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

% dimensions
sX = size(Na);

% reahape inputs
P = zeros(sX);
for indx = 1:length(Pv)
    [indb, inda] = ind2sub(sX, indx);
    P(indb, inda) = Pv(indx);
end

% rates
Wap = fncWap(Na, Nb, lambda*r);
Wbp = fncWbp(Na, Nb, lambda*r);
Wbm = fncWbm(Na, Nb, gamma);
Wabp = fncWambp(Na, Nb, lambda*(1-2*r));
za = zeros(sX(1),1);
zb = zeros(1, sX(2));

% Master equation
dPdt = [za Wap(:,1:end-1).*P(:,1:end-1)] + ...
    [zb; zb; [Wbp(1:end-2,2:end).*P(1:end-2,2:end) za(1:end-2)]] + ...
    [Wbm(2:end,:).*P(2:end,:); zb] + ...
    [zb; Wabp(1:end-1,:).*P(1:end-1,:)] + ...
    -(Wap + Wbp + Wbm + Wabp).*P;

% reshape output
dPvdt = dPdt(:);

end


function Wap = fncWap(na, ~, r) % transition from nanb to na+1
    Wap = r*na;
end

function Wbp = fncWbp(na, ~, r) % transition from nanb to na-1nb+2
    Wbp = r*na;
end

function Wbm = fncWbm(~, nb, r) % transition from nanb to nb-1
    Wbm = r*nb;
end

function Wabp = fncWambp(na, ~, r) % transition from nanb to nanb+1
    Wabp = r*na;
end
