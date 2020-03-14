function logL = sta_likelihoodCloneSize(Data, Model)
%   logL = sta_likelihoodCloneSize(Data, Model)
% This function estimate the Likelihood function, which is the probability
% that a model reproduces the observed data. 
% The model and the data represent clonal size distributions at different
% time points. At the i-time point the likelihood follows a multinomial
% distribution given by:
%   L_i(Model_i|Data_i) = P(Data_i|Model_i) = Sum(f_n)!/Prod(f_n!)*Prod(p_n^f_n)
% where Data_i = {f_n}, corresponding to the frequency of cloned of size n in the observed data,
% and Model_i = {p_n}, corresponding to the probability of cloned of size n in the simulated data.
% Likelihood is then computed as
%   L(Model|Data) = Prod(L_i(Model_i|Data_i))
% Logarithmic scale is used for simplify the computation and for returning
% the output.
% 
%  Data     IN: Observed clonal size data at each time point, cell(1,N)
%  Model    IN: Model (simulated clonal size data) at each time point, cell(Nsim,N)
%  varargin{1} IN: stat parameters structure
%  logL    OUT: Likelihood of the model in log10 scale (NsimX1)

% initialization
Nt = length(Data);
Nsim = size(Model,1);
logLi = zeros(Nsim, Nt);

% loop on time point
for jj = 1:Nsim
    for ii = 1:Nt
        
        % clonal size data at the ii-time step
        xm = Model{jj,ii};
        xd = Data{ii};
        
        % define common bin edges
        % edges = floor(min([xm; xd]))-0.5:1:ceil(max([xm; xd]))+0.5;
        edges = floor(min(xd))-0.5:1:ceil(max(xd))+0.5; % only in the range of the observations
        
        % probability model/frequency data
        pm = histcounts(xm, edges, 'Normalization', 'probability');
        fd = histcounts(xd, edges, 'Normalization', 'count');
        
        % Likelihood
        % ----------
        % a = factorial(sum(fd)) = 1*2*...*Sfd
        % log(a) = log(1) + log(2) + ... + log(Sfd)
        Sfd = sum(fd);
        loga = sum(log10(1:1:Sfd));
        % b = prod(factorial(fd)) = (1*2*...*fd(1))*(1*2*...*fd(2))*...*(1*2*...*fd(end))
        % logb = log(1*2*...*fd(1)) + log(1*2*...*fd(2)) + ... + log(1*2*...*fd(end))
        %      = log(1)+log(2)+...+log(fd(1))+...+log(1)+log(2)+...+log(fd(end))
        logb = 0;
        for ifd = 1:length(fd)
            logb = logb + sum(log10(1:1:fd(ifd)));
        end
        % c = prod(pm.^fd);
        % logc = sum(log(pm.^fd)) = sum(fd*log(pm))
        inz = fd>0; % we need to exclude the points in which we don't have observations otherwise likelihood is 0
        if isempty(find(inz>0,1)) || any(pm == 0) % if there are no points in common or in one point the model probability is zero, returns NaN
            logc = NaN;
        else
            logc = sum(fd(inz).*log10(pm(inz)));
        end
        % log of likelihood
        logLi(jj,ii) = loga - logb + logc;
    end
end

% merge all the contributions at different time points
logL = sum(logLi, 2);

end