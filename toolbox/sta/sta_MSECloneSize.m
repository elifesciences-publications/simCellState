function [mse, xmse] = sta_MSECloneSize(Data, Model, varargin)
%   [logL, xlogL] = sta_likelihoodCloneSize(Data, Model, varargin)
% This function estimate the Mean Square Error (MSE) of the observed VS
% model data, which represent clonal size distributions at different
% time points. At the i-time point the MSE is computed; output values are
% the sum of the MSE at different time points.
% Confidence interval is also estimated if, for each time step, the size of
% Model is greater than 1. In this case the output corresponds to the
% percentiles [50-CI/2 50 50+CI/2]
% 
%  Data     IN: Observed clonal size data at each time point, cell(1,N)
%  Model    IN: Model (simulated clonal size data) at each time point, cell(Nsim,N)
%  varargin{1} IN: Confidence interval (%), optional, default 90%
%  mse     OUT: Mean Square Error. If Nsim > 1, the output corresponds to
%               [50-CI/2 50 50+CI/2]%ile
%  xmse    OUT: Mean Square Error. (all values in case  Nsim >1)

% initialization
Nt = length(Data);
Nsim = size(Model,1);
xxMSE = NaN(Nsim, Nt);
if nargin == 2
    CI = 90;
else
    CI = varargin{1}; % Confidence interval (%)
end

% loop on time point
for jj = 1:Nsim
    for ii = 1:Nt
        
        % clonal size data at the ii-time step
        xm = Model{jj,ii};
        xd = Data{ii};
        
        % define common bin edges
        edges = floor(min([xm; xd]))-0.5:1:ceil(max([xm; xd]))+0.5;
        
        % probability model/frequency data
        pm = histcounts(xm, edges, 'Normalization', 'probability');
        pd = histcounts(xd, edges, 'Normalization', 'probability');
        
        % mean square error
        xxMSE(jj, ii) = immse(pm, pd);
        
    end
end

% sum MSE in time
mse = sum(xxMSE, 2);

if Nsim > 1 % mean and confidence interval
    pval = 50 + [-CI/2 0 CI/2];
    xmse = mse;
    mse = prctile(xmse, pval);
else
    xmse = [];
end

end