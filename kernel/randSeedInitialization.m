function [SDs, SDt, SDx] = randSeedInitialization(SD0, Nrun, Nrs, iX0r, ipar, NrunMax)
%   [SDs, SDt] = randSeedInitialization(SD0, Nrun, Nrs, iX0r, ipar, NrunMax)
% This function generate the inital seed for initializing random values in
% each simulation. Random seeds corresponds to a series of odd values that
% starts from SD0.
%
%  SD0   IN: inital seed
%  Nrun   IN: number of repeated Nrun
%  Nrs  IN: number of simulations
%  iX0r  IN: flag for random intial condition, if empty SDx is zero
%  ipar  IN: flag for parallel computation
%  NrunMax IN: max size allowed for vectorized version in case of parallel (Nrun overwritten to 1)
%  SDs  OUT: matrix for initializing random number for state (Nrun, Nrs)
%  SDt  OUT: matrix for initializing random number for time (Nrun, Nrs)
%  SDx  OUT: matrix for initializing random number for initial condition (1, Nrs)
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
% initialization
if Nrun < NrunMax
    Nrun = 1;
end
SDs = zeros(Nrun, Nrs);
SDt = zeros(Nrun, Nrs);
SDx = zeros(Nrun, Nrs);

SD = SD0;
for irs = 1:Nrs % loop on Nrs
    for irun = 1:Nrun % loop on Nrun
        SD = SD+2;
        SDs(irun, irs) = SD;
        SD = SD+2;
        SDt(irun, irs) = SD;
    end
    if iX0r
        SD = SD+2;
        SDx(irs) = SD;
    end
end