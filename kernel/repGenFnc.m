function F = repGenFnc(s, p0, p1, p2)
%  F = repGenFnc(s, p0, p1, p2)
% This function evaluate the extinction generating function defined by the
% parameters p0, p1 and p2.
%  p0  IN: constant component (N,1)
%  p1  IN: linear component (N,N)
%  p2  IN: quadratic component rearranged to have dimensions N, and N*N (N,N*N)
%  F  OUT: reproduction generating function (N,1)
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

% zero and first order term
s2 = s*s';
F = p0 + p1*s + p2*s2(:);