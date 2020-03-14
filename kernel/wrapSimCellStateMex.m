function [xNout, tExt, istatus] = wrapSimCellStateMex(~, X0, T0, xtime, ut, us, sRate, indxState, deltaState, indxActiveState)

% wrapper function for calling the mex version of simCellState function.
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

% call mex
[xNout, tF, istatus] = simCellState_mex(T0, xtime, X0', ut, us, sRate', indxState', deltaState, indxActiveState);
if istatus == -2 % extinction
    tExt = tF;
    istatus = 0;
else
    tExt = NaN;
end
if istatus == -1 % all rand used
   xNout(:, xtime>tF) = NaN; 
end

end