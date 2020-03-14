function [xOut, ierr] = modifyRates(xOut, xProp, Rmean, Rmin, eigTg, eigTol)

% This function applies some random rates to the non-zero element of the
% division and transition matrices (xOut{jj, 1/2}).
% Rates are exponentially distributed with mean Rmean and greater than
% Rmin.
% In case eigTol is set, then the rates are further tuned to meet the
% requirement on the maximum eigenvalue (<=eigTol).
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

itg = 0; imax = 0;
if nargin > 4
    if eigTg <0
        imax = 1;
    else
        itg = 1;
    end
end

% initialization
ierr = zeros(1, length(xOut));
% loop on size of SCC (ii) and networks (jj)
for ii = 1:length(xOut)
    Ns = xProp(ii).NSCC+1;
    for jj = 1:size(xOut{ii},1)
        % modify cell division rates
        Cd = xOut{ii}{jj,1}; indxCd = find(sum(Cd)>0);
        xr = log(1-rand(1,length(indxCd)))/(-Rmean-Rmin)+Rmin;
        rr = ones(1, Ns); rr(indxCd) = xr; rr = repmat(rr, Ns, 1);
        Cd1 = Cd.*rr; % division
        % modify cell transition rates
        Ct = xOut{ii}{jj,2}; indxCt = find(Ct>0);
        xr = log(1-rand(1,length(indxCt)))/(-Rmean-Rmin)+Rmin;
        rr = ones(Ns, Ns); rr(indxCt) = xr;
        Ct1 = Ct.*rr; % transition
        % updating max eigenvalue
        [~, J] = CdCt2AJ(Cd1, Ct1);
        eg = eig(J(1:Ns-1,1:Ns-1)); [~, im] = max(real(eg));
        egMax1 = eg(im);
        if itg || (imax && egMax1 > eigTg)
            [Cd2, Ct2, egMax2] = modifyRateEig(Cd1, Ct1, eigTg, Rmin);
            % ierr
            if abs(egMax2-eigTg) > eigTol
                ierr(ii) = -1;
            end
        else
            Cd2 = Cd1; Ct2 = Ct1; egMax2 = egMax1;
        end
        % update output
        xOut{ii}{jj,1} = Cd2;
        xOut{ii}{jj,2} = Ct2;
        xOut{ii}{jj,3} = egMax2;
        if any(Cd2(Cd2>0) < Rmin) || any(Cd2(Cd2>0) < Rmin)
            ierr(ii) = -2;
        end            
    end
end

end