function outStruct = findStableRates(Cd, Ct, iScc)
%   outStruct = findStableRates(Cd, Ct, iScc)
% This function finds a set of parameters that stabilizes a given system. 
% Assumptions:
% - all the parameters are perturbed
% - only one strongly connected component is analyzed: this has to be generalized
% 
%  Cd    IN: Division matrix (N,N)
%  Ct    IN: Transition matrix (N,N)
%  iScc  IN: strongly connected component
%  outStruct  OUT: structure containing the solutions. The fields are:
%        - exitflag: exit flag condition (see stabilityMO)
%        - fstab: stability measure (0 is stable)
%        - xMO: cell containing the updated M and O matrices

% system properties
sysProp = sysProperties('CdCt', Cd, Ct);

% non zero entries of the M and O matrices
indxT = find(Ct);
nT = length(indxT); % 0-N*(N-1)
indxD = find(sum(Cd)); % 0-N
nD = length(indxD);

if sysProp.NSSC == 1 % assuming just 1 SCC and unstable
    inode = sysProp.SCC == iScc; % just one SCC
    Nvar = nT + nD;
    
    % global optimization
    if ~sysProp.stability.flag
        options = optimoptions(@ga, 'display', 'off');
        [x,fval] = ga(@(x) stabilitySCCper(Cd, Ct, inode, indxD, indxT, x), Nvar, [], [], [], [], ...
            -0.9*ones(1,Nvar),1.5*ones(1, Nvar), [], options);
        % local search
        [fstab, dCd, dCt, exitflag] = localSearchStableRates(Cd, Ct, indxD, indxT, inode, x);
        Cts = Ct+dCt;
        Cds = Cd+dCd;
        xCdCt = {Cds; Cts};
    else
        fstab = 0; exitflag = 1;
        xCdCt = {Cd; Ct}; % already stable
    end
else
    xCdCt = [];
    exitflag = NaN;
    fstab = NaN;
end

% output structure
outStruct = struct('exitflag', exitflag, 'fstab', fstab);
outStruct.xCdCt = xCdCt;

return