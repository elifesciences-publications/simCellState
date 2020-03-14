function [xx1, xxpdf] = getPdf(nn_sim, dx)

% This function estimate the continuous normalized pdf based on the
% cumulative count of the input nn_sim.
% Output is given at steps dx.

% cumulative count
[nncc, nnccedges] = histcounts(nn_sim, 'BinLimits', [0 max(nn_sim)], 'BinWidth', 1,'Normalization', 'cumcount');

% normalize
xx0 = (nnccedges(1:end-1)+diff(nnccedges)/2)/mean(nn_sim);
xx1 = xx0:dx:max(xx0); 
nncdf1 = interp1(xx0, nncc/max(nncc), xx1);

% pdf (based on finite difference 2nd order)
xxpdf = varDer2(xx1, nncdf1);
xxpdf(xxpdf<0) = 0; % numerical errors

end
