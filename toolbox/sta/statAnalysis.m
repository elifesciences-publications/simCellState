function outStat = statAnalysis(kf, indxActiveState, indxSurvIn)
% STATISTICS ANALYSIS
% outStat = statAnalysis(kf, indxActiveState, indxSurvIn)
%
% This function computes the statistics of the given set of data kf. NaN
% (representing data not reaching the final time) are always excluded.
% If data is 3D, last dimension is considered as time. 
%
%  kf         IN: data, each row is a variable, each column is a sample.
%                 Third dimension is time (M,N,L)
%  indxActiveState IN: active state index used to computed the survived
%                 cases (1, K)
%  indxSurvIn IN: cell containing the index of the survived sample at each time (1,L)
%  outStat    OUT: structure containing
%                 - mean and std (M,L)
%                 - meanT and stdT corresponding to the total of the active states (M,L)
%                 - distr: structure containing xx, pdf, cdf profiles (M,L)
%                 - distrT: structure containing xx, pdf, cdf profiles of
%                 the totality of the active states (1,L )
%                 - index of not-NaN cases (N,L)
%                 - index of survived cases (N,L)

% set dimensions
sk = size(kf);
Nvar = sk(1);
if length(sk) == 3
    Ntime = sk(3);
else
    Ntime = 1;
end

% initialization
mm = NaN(Nvar, Ntime); ss = NaN(Nvar, Ntime);
mmT = NaN(1, Ntime); ssT = NaN(1, Ntime);
indxNotNan = cell(1,Ntime); indxSurv = cell(1,Ntime);
distr = struct('xx', {}, 'pdf', {}, 'cdf', {});
distrT = struct('xx', {}, 'pdf', {}, 'cdf', {});

% indentify notNaN and survived cases
for itime = 1:Ntime
    % indexes
    indx = find(~isnan(sum(kf(:,:,itime), 1)));
    indxNotNan{itime} = indx;
    if isempty(indxSurvIn)
        indxSurv{itime} = find(~(sum(kf(indxActiveState,:,itime), 1) == 0)); % survived & not NaN
    else
        indxSurv{itime} = indxSurvIn{itime}; % from input
        indx = intersect(indxNotNan{itime}, indxSurv{itime});
    end
    
    % compute mean and std: all (excluding NaN)
    if ~isempty(indx)
        % mean and std
        mm(:,itime) = mean(kf(:, indx, itime),2);
        ss(:,itime) = std(kf(:, indx, itime),0,2);
        % distribution
        for im = 1:size(kf,1)
           [distr(im, itime).xx, distr(im, itime).pdf, distr(im, itime).cdf] = computeDitribution(kf(im, indx, itime));
        end
        % mean, std and distribution of the sum
        xsum = sum(kf(indxActiveState, indx, itime), 1);
        mmT(itime) = mean(xsum,2);
        ssT(itime) = std(xsum,0,2);
        [distrT(itime).xx, distrT(itime).pdf, distrT(itime).cdf] = computeDitribution(xsum);
    end
    
end

% set output
outStat.indxNotNan = indxNotNan;
outStat.indxSurv = indxSurv;
outStat.mean = mm;
outStat.std = ss;
outStat.distr = distr;
outStat.distrT = distrT;
outStat.meanT = mmT;
outStat.stdT = ssT;

end


% buildin function
function [xx, pdf, cdf] = computeDitribution(kf)

% definition of the bins
% edges = floor(min(kf))-0.5:1:max(kf)+0.5;
dbin = max([ceil(kf/50) 1]);
edges = floor(min(kf))-0.5:dbin:max(kf)+dbin; % max 50 bins to improve the resolution, at least 1 bin

% pdf
pdf = histcounts(kf,edges,'Normalization','pdf');
cdf = histcounts(kf,edges,'Normalization','cdf');
% xx
xx = edges(1:end-1)+diff(edges)/2;

end
