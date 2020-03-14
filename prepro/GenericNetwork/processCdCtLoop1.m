function [xOutF, xOut] = processCdCtLoop1(Cd, Ct, icons)

% This function calls processCdCt and store the outputs in matrix form
% where rows corresponds to Cd and columns to Ct. The full output (xOut) is then
% filtered (xOutF).
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

N = sqrt(size(Cd, 1))-1;
Ncd = size(Cd, 2); Nct = size(Ct,2);
xOut = cell(Ncd*Nct, 4);
kk = 0;
for ii = 1:Ncd
    Cdii = reshape(Cd(:,ii), N+1, N+1);
    for jj = 1:Nct
        Ctjj = reshape(Ct(:,jj), N+1, N+1);
        kk = kk+1;
        xOut(kk,:) = processCdCt(Cdii, Ctjj);
    end
end
% filter solutions
xOutF = filterOut(xOut, icons);

return
