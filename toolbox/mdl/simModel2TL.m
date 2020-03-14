function [T, L, Nstate] = simModel2TL(simModel)
%  [T, L, Nstate] = simModel2TL(simModel)
% This function is a wrapper to generate the sparse matrix (deterministic
% model) compatible with the inputRead function outputs.
% It is assumed that each colum can sum at max 2 (corresponding to cell
% splitting option). All rates are qual to 1 s^-1.
%  simModel     IN: model info
%  inputT      OUT: Transitions
%  inputS      OUT: Splitting
%  Nstate      OUT: Number of state

% parameters
Nstate = simModel.Nstate;
indxState = simModel.indxState;
deltaState = simModel.deltaState;
sRate = simModel.sRate;

% initialization
T = []; L = [];
% addDS = zeros(Nstate); addDS(sub2ind([Nstate Nstate], indxState, indxState)) = 1;
% deltaState = deltaState + addDS;
for irate = 1:length(sRate)
    iIn = indxState(irate);
    dsti = deltaState(:,irate);
    dsti(iIn) = dsti(iIn)+1;
    iOut = find(dsti);
    if length(iOut) == 2
        L = cat(1, L, [sRate(irate) iIn iOut']);
    else
        T = cat(1, T, [sRate(irate) iIn iOut]);
    end
end

end