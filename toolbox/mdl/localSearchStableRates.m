function [fstab, dCd, dCt, exitflag] = localSearchStableRates(Cd, Ct, indxCd, indxCt, inode, x0, varargin)
%  [fstab, dCd, dCt, exitflag] = localSearchStableRates(Cd, Ct, indxCd, indxCt, inode, x0)
% This function find the perturbations to be added to Cd and Ct matrices in
% the indxCd and indxCt components.
%  Cd     IN: deterministic matrix
%  Ct     IN: stochastic matrix
%  indxCd IN: index of the column of the deterministic matrix to be perturbed
%  indxCt IN: index of stochastic matrix perturbations
%  inode IN: index of the strongly connected component nodes
%  x0    IN: initial value of the perturbation
%  varargin IN: optional argument for calling stabilitySCCper function
%  fstab OUT: stability measure, prod(abs(eig)), 0 means it is stable. If
%       perturbed rates are equal or less than zero or exitflag is not 1
%       then fstab is set to Inf
%  dCd    OUT: deterministic matrix perturbation
%  dCt    OUT: stochastic matrix perturbation
%  exitflag OUT: 1 fmincon converged, stable, positive rates
%                0 fmincon converged, stable, rate larger than 100 s-1
%                -1 fmincon converged, stable, negative or zero rate
%                -2 fmincon converged, unstable solution
%                -3 fmincon not converged
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

% update here
% fminfunc
tolFun = 1e-15; 
options = optimset('TolFun', tolFun, 'Display', 'off');
dstabFnc = @(x) stabilitySCCper(Cd, Ct, inode, indxCd, indxCt, x, varargin{:});
[x,fstab,exitflag] = fminsearch(dstabFnc, x0, options);
% [x,fstab,exitflag] = fmincon(dstabFnc, x0, [], [], [], [], x0, [], [], options);

% matrix perturbations
dCd = zeros(size(Cd));
for im = 1:length(indxCd)
    dCd(Cd(:,indxCd(im))~=0, indxCd(im)) = x(im);
end
dCt = zeros(size(Ct)); dCt(indxCt) = x(1+length(indxCd):end);
Cd = Cd + dCd;
Ct = Ct + dCt;

% exit flag
tolStab = 1e-10; tolRate = 1e-3; maxRate = 100;
if exitflag == 1
    if fstab < tolStab
        if any(Ct(indxCt) >= maxRate) || any(sum(Cd(:,indxCd)) > maxRate)
            exitflag = 0;
        end
        if any(Ct(indxCt) <= tolRate) || any(sum(Cd(:,indxCd)) <= tolRate)
            exitflag = -1;
        end
    else
        exitflag = -2;
    end
else
    exitflag = -3;
end

end

