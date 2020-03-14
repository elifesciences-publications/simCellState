function [pd, h, p, st] = sta_fitTestDistrib(n, type, iplot, min_thrs, ibin, nbin, eminpool)

% This function fit and run the chi2test to data (n) based on a given
% distribution (type = poi, norm, geo, exp). Optionally also plot the
% fitting results. The output is the distribution of the observed data
% (pd), and the chi2gof test outputs (h, p, st).
% min_thrs is a cutoff in case of expo
% if ibin == 1, bin size is 1, otherwise nbin has to be specified, where
% ibin = max(1, floor(mean/nbin))
% eminpool cutoff if observation < eminpool (use at least 5 for limits in chi2test)

% mean and std of data
muii = mean(n);
stdii = std(n);

% create fitted distribution
switch type
    case 'poi'
        pd = makedist('Poisson', 'lambda', muii);
    case 'nor'
        pd = makedist('Normal', 'mu', muii, 'sigma', stdii);
    case 'geo'
        pd = makedist('NegativeBinomial','R',1,'p',1/muii);
    case 'exp'
        pd = makedist('Exponential', 'mu', muii);
    otherwise
        error('check type input')
end

% observation and expected counting (pooling data when expected count < eminpool)
[binsp, edgesp] = sta_poolBins(n, pd, min_thrs, ibin, nbin, eminpool);
if length(binsp) >= 2
    [obsCounts, expCounts] = sta_obsExpCounts(n, pd, binsp, edgesp);
    
    if max(expCounts) > eminpool
        % adjust settings depending on reference distribution
        emin = 0; % not using internal pooling function --> controlled bins
        Nparams = 1; % default
        switch type
            case 'norm'
                Nparams = 2;
        end
        % chi2 test
        [h, p, st] = chi2gof(binsp, 'Edges', edgesp, 'Frequency', obsCounts, 'expected', expCounts, 'NParams', Nparams, 'emin', emin);
        
        % plot
        if iplot == 1
            xx = st.edges(1:end-1)+ diff(st.edges)/2;
            switch type
                case {'geo', 'poi'}
                    xx = round(xx);
            end
            figure; hold on; grid on
            plot(xx, expCounts, 'b', binsp, obsCounts, 'b*')
            plot(xx, st.E, 'r', xx, st.O, 'r.')
            xlabel('x'); ylabel('N count'); legend('Expected (in)', 'Observed (in)', 'Expected (out)', 'Observed (out)')
        elseif iplot > 1
            xcol = lines(iplot-1); xcol = xcol(end,:); 
            xx = st.edges(1:end-1)+ diff(st.edges)/2;
            h1 = plot(xx, st.E/length(n)*100, 'k-');
            h2 = plot(xx, st.O/length(n)*100, '.', 'MarkerSize', 10); set(h2, 'color', xcol);
            xlabel('n'); ylabel('Frequency [%]');
            if p<0.01
                tg = sprintf('%s = %.0f, p-val = %.2e', '\mu', muii, p);
            else
                tg = sprintf('%s = %.0f, p-val = %.2f', '\mu', muii, p);
            end
            legendHandle = legend([h2, h1], 'Observed', 'Expected');
            title(legendHandle, tg)
        end
    else
        h = NaN;
        p = NaN;
        st = [];
    end
else
    h = NaN;
    p = NaN;
    st = [];
end
end

function [binsp, edgesp] = sta_poolBins(n, pd, min_thrs, ibin, nbin, emin)

% bins covering the whole range
if ~isempty(ibin)
    bins1 = min(n):1:max(n);
else % variable binning mean/ibin
    db = max(1, floor(mean(n)/nbin));
    bins1 = min(n):db:max(n);
end

% real count
exp = histcounts(n, bins1); % length(n)*pd.pdf(bins1);
% find index for pooling tails
i1 = find(exp>emin, 1, 'first');
i2 = find(exp>emin, 1, 'last');
% edges
if ~isempty(i1)
    edgesp = bins1(i1:i2-1)+0.5;
else
    edgesp = bins1;
end
% pooled bin (centres)
binsp = edgesp(1:end-1) + diff(edgesp)/2;

% remove additional parts
switch pd.DistributionName % pdf integration
    case {'Negative Binomial', 'Exponential'}
        edgesp(binsp<min_thrs*mean(n)) = [];
        binsp = edgesp(1:end-1) + diff(edgesp)/2;
end

end

function [obsCounts, expCounts] = sta_obsExpCounts(n, pd, binsp, edgesp)

% observed counts
obsCounts = histcounts(n, edgesp); % 'BinWidth',1, 'BinLimits',[min(n)-0.5,max(n)+0.5]);

% expected counts
expCounts = zeros(1,length(binsp));
switch pd.DistributionName % pdf integration
    case 'Poisson'
        for ibin = 1:length(binsp)
            xx = unique(ceil(edgesp(ibin):1:floor(edgesp(ibin+1))));
            expCounts(ibin) = sum(pd.pdf(xx));
        end
    case 'Negative Binomial'
        for ibin = 1:length(binsp)
            xx = unique(ceil(edgesp(ibin):1:floor(edgesp(ibin+1))));
            expCounts(ibin) = sum(pd.pdf(xx-1)); % shift by 1
        end
    case 'Normal'
        for ibin = 1:length(binsp)
            xx = edgesp(ibin):0.01:edgesp(ibin+1);
            expCounts(ibin) = trapz(xx,pd.pdf(xx));
        end
    case 'Exponential'
        for ibin = 1:length(binsp)
            xx = edgesp(ibin):0.01:edgesp(ibin+1);
            expCounts(ibin) = trapz(xx,pd.pdf(xx)); % trapz(xx,pd.pdf(xx-0.5)); % shift by 0.5
        end
end
expCounts = expCounts*length(n);

end
