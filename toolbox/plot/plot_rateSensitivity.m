function hfg = plot_rateSensitivity(Cd, Ct)
%  hfg = plot_rateSensitivity(Cd, Ct)
% This function plots the sensitivity of the stability of the system to
% each rate. Only the first strongly connected component is analyzed.
% 
%  Cd    IN: division matrix
%  Ct    IN: transition matrix
%  hfg   OUT: figure handle

% non zero entries of the Cd and Ct matrices
indxCt = find(Ct);
nCt = length(indxCt); % 0-N*(N-1)
indxCd = find(sum(Cd)); % 0-N
nCd = length(indxCd);

% system properties
[A, J] = CdCt2AJ(Cd, Ct);
sysProp = sysProperties('CdCt', Cd, Ct);
inode = sysProp.SCC == 1; % just one SCC

% rate perturbation
dr = -10:0.01:100;
dstabFnc = zeros(size(dr));
xleg = cell(1, nCt+nCd);
h = zeros(1, nCt+nCd);
hfg = figure; hold on; grid on
for io = 1:nCt
    for ii = 1:length(dr)
        dstabFnc(ii) = stabilitySCCper(Cd, Ct, inode, [], indxCt(io), dr(ii));
    end
    R0 = Ct(indxCt(io));
    h(io) = plot(dr+R0, dstabFnc);
    [I,J] = ind2sub(size(Cd), indxCt(io));
    xleg{io} = sprintf('%s_{%d%d}', '\omega', J, I);
end
for im = 1:nCd
    for ii = 1:length(dr)
        dstabFnc(ii) = stabilitySCCper(Cd, Ct, inode, indxCd(im), [], dr(ii));
    end
    R0 = sum(Cd(:,indxCd(im)))/2; Rs = 3-sum(Cd(:,3)~=0); % ii->2*jj or ii->jj+kk
    h(im+nCt) = plot(dr/Rs+R0, dstabFnc);
    xleg{im+nCt} = sprintf('%s_{%d}', '\lambda', indxCd(im));
end
xlabel('rate [s^{-1}]'); ylabel('min(abs(eig))')
legend(h, xleg);
