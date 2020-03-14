function [dstab, eigJ] = stabilitySCCper(Cd, Ct, inode, indxCd, indxCt, dr, varargin)
%  dstab = stabilitySCCper(Cd, Ct, dCd, dCt, inode, dr)
% This function computes the product of the module of the
% eigenvalues of the part of the Jacobian related to a strongly connected
% component. 
% Equilibrium of the SCC results in case dstab is zero (meaning that 1
% eigenvalue is zero).
%  Cd    IN: Division matrix
%  Ct    IN: Transition matrix
%  inode IN: index of the strongly connected component nodes
%  indxCd IN: index of the column of the deterministic matrix to be perturbed 
%  indxCt IN: index of stochastic matrix perturbations
%  dr    IN: perturbation
%  varargin IN: optional argument for calling stabilitySCC function
%  dstab OUT: distance from stable condition
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
% perturbed inputs
dCd = zeros(size(Cd)); 
for im = 1:length(indxCd)
    dCd(Cd(:,indxCd(im))~=0, indxCd(im)) = dr(im);
end
dCt = zeros(size(Ct)); dCt(indxCt) = dr(1+length(indxCd):end);
Cd = Cd + dCd;
Ct = Ct + dCt;

% stability
[~, J] = CdCt2AJ(Cd, Ct);
eigJ = eig(J(inode, inode));

% stabilty
dstab = stabilitySCC(eigJ, varargin{:});

end