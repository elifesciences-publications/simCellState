function [value,isterminal,direction] = feedback_odeEvent(t, x, indexActiveState)

% For the I-th event function: 
%   VALUE(I) is the value of the function, ISTERMINAL(I)=1 if the integration 
%   is to terminate at a zero of this event function and 0 otherwise. 
%   DIRECTION(I)=0 if all zeros are to be computed

direction = 0; isterminal = 1;
if sum(x(indexActiveState)) > 1e4 || sum(x(indexActiveState)) <= 0.1
    value = 0;
else
    value = 1;
end