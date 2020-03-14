function [value,isterminal,direction] = odeEvent(t, x, indexActiveState, nMin, nMax)
% For the I-th event function: 
%   VALUE(I) is the value of the function, ISTERMINAL(I)=1 if the integration 
%   is to terminate at a zero of this event function and 0 otherwise. 
%   DIRECTION(I)=0 if all zeros are to be computed
% This specific function stop the integration if the mean total number of
% cells is below nMin or above Nmax
% 
% Cristina Parigini, 14/03/2020
% 
% Copyright 2020 Cristina Parigini
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

direction = 0; isterminal = 1;
if sum(x(indexActiveState)) > nMax || sum(x(indexActiveState)) <= nMin
    value = 0;
else
    value = 1;
end