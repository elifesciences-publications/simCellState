function dstab = stabilitySCC(eigJ, varargin)
%  dstab = stabilitySCC(eigJ, varargin)
% This function measure the stability of a SCC based on the eigenvalues of
% the jacobian. The output is the absolute value of the eigenvalue with
% maximum module with sign. This means that:
% - if all eigenvalues have negative real part, it returns the module of the largest one.
% - if all are negative and one is positive, it returns the module of the positive one.
% - if there is one (or more) zero eigenvalue and all negative, it returns zero.
% - if there is one (or more) zero eigenvalue, one (or more) positive, it returns the largest positive one.
% Optionally a shift can be passed as input (to force the SCC being
% subcritical).
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

if nargin == 2
    shift = varargin{1};
    if isempty(varargin{1})
        shift = 0;
    end
else
    shift = 0;
end
% output
dstab = abs(max(real(eigJ))-shift);