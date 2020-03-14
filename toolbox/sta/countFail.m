function [xmus, nfail, ntot] = countFail(musurv, xp, dmu, pThrs)
%  [xmus, nfail, ntot] = countFail(musurv, xp, dmu, pThrs)
% This function count the number of failed cases as function of the mean
% value of the survived cases. Samples are binned toghether based on the
% mean.
%  musurv  IN: mean of surviving clones matrix (Nsim X Ntime)
%  xp      IN: p-value matrix (Nsim X Ntime)
%  dmu     IN: output mean step
%  pThrs   IN: p-value threshold
%  xmus    OUT: mean (N)
%  nfail   OUT: number of failing cases (N)
%  ntot    OUT: number of total cases (N)

% evaluate number of failed cases
xmus = max(min(min(floor(musurv))),2):dmu:max(max(ceil(musurv))); % remove mu = 1 (initial condition)
musurvr = discretize(musurv, [0 xmus(1:end-1)+diff(xmus)/2 xmus(end)], xmus); % binned values
% initialization
ntot = NaN(size(xmus));
nfail = NaN(size(xmus));
% loop
for imu = 2:length(xmus) % first point removed
    ii = find(musurvr == xmus(imu));
    ntot(imu) = length(ii);
    nfail(imu) = length(find(xp(ii) < pThrs));
end

end