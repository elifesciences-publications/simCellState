function [nn1, nnpmf] = getPmf(nn_sim, dn)

% This function estimate the pmf of the nn_sim data.

% histcounts
[nnc, nncedges] = histcounts(nn_sim, 'BinLimits', [0-dn/2 ceil(max(nn_sim)/dn)*dn+dn/2], 'BinWidth', dn, 'Normalization', 'count');

% normalize
nn1 = (nncedges(1:end-1)+diff(nncedges)/2);
nnpmf = nnc/length(nn_sim);

end
