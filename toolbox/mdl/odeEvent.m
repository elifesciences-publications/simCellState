function [value,isterminal,direction] = odeEvent(t, x, indexActiveState, nMin, nMax)
% For the I-th event function: 
%   VALUE(I) is the value of the function, ISTERMINAL(I)=1 if the integration 
%   is to terminate at a zero of this event function and 0 otherwise. 
%   DIRECTION(I)=0 if all zeros are to be computed
% This specific function stop the integration if the mean total number of
% cells is below nMin or above Nmax

direction = 0; isterminal = 1;
if sum(x(indexActiveState)) > nMax || sum(x(indexActiveState)) <= nMin
    value = 0;
else
    value = 1;
end