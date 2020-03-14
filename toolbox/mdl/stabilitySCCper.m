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