function outStruct = findRelevantParameters(Cd, Ct, iScc, ipar, varargin)
%   outStruct = findRelevantParameters(Cd, Ct, iScc, ipar, varargin)
% This function finds the minimum set of parameters that stabilizes a given
% system. All the combinations (ipar number of parameters) are tested.
% Assumptions:
% - all the parameters are perturbed
% - only one strongly connected component is analyzed: this has to be generalized
% 
%  Cd    IN: Division matrix (N,N)
%  Ct    IN: Transition matrix (N,N)
%  iScc  IN: strongly connected component
%  ipar  IN: number of parameters to be perturbed. If empty, all the
%           combinations with any number of parameter are tested.
%  varargin  IN: additional optional input to stabilitySCCper and
%           localSearchStableRates functions
%  outStruct  OUT: structure containing the solutions. The fields are:
%        - exitflag: exit flag condition (see stabilityMO)
%        - fstab: stability measure (0 is stable)
%        - xMO: cell containing the updated M and O matrices

% optional inputs
ishift = 0;
if nargin > 4
    ishift = 1;
end

% system properties
sysProp = sysProperties('CdCt', Cd, Ct);

% number of parameters
indxT = find(Ct);
nT = length(indxT); % 0-N*(N-1)
indxD = find(sum(Cd)); % 0-N
nD = length(indxD);
Nvar = nT + nD;
if isempty(ipar)
    ipar = 1:Nvar; % all possiblility are tested
end

% initialization
xCdCt = []; xexitflag = []; xfstab = [];

if sysProp.NSSC == 1 % assuming just 1 SCC and unstable
    inode = sysProp.SCC == iScc; % just one SCC
    if ~sysProp.stability.flag || ishift
        for ipar = ipar
            if ~isempty(xexitflag) % just the minimum set of parameters
                break
            end
            C = nchoosek([-indxD'; indxT], ipar);            
            for isol = 1:size(C,1)
                Csol = C(isol,:);
                % global search
                LB = -0.9*ones(1,ipar);
                UB = 1.5*ones(1,ipar);
                options = optimoptions(@ga, 'display', 'off');
                [x,fval] = ga(@(x) stabilitySCCper(Cd, Ct, inode, -Csol(Csol<0), Csol(Csol>0), x, varargin{:}), ipar, [], [], [], [], ...
                    LB, UB, [], options);
                % local search
                [fstab, dCd, dCt, exitflag] = localSearchStableRates(Cd, Ct, -Csol(Csol<0), Csol(Csol>0), inode, x, varargin{:});
                if exitflag == 1
                    Cts = Ct+dCt;
                    Cds = Cd+dCd;
                    xCdCt = [xCdCt {Cds; Cts}];
                    xexitflag = [xexitflag exitflag];
                    xfstab = [xfstab fstab];
                end
            end
        end
    else
        xfstab = 0; xexitflag = 1;
        xCdCt = {Cd; Ct}; % already stable
    end
else
    xCdCt = [];
    xexitflag = NaN;
    xfstab = NaN;
end

% output structure
outStruct = struct('exitflag', xexitflag, 'fstab', xfstab);
outStruct.xCdCt = xCdCt;

return