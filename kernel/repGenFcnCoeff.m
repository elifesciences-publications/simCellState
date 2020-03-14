function [p0, p1, p2, rtot] = repGenFcnCoeff(T, L, N)
%   [p0, p1, p2] = repGenFcnCoeff(T, L, N)
% This function computes the coefficients of the reproduction generating
% function. It is assumed that F(s) = p0 + p1*s + p2*s*s' (see
% repGenFnc.m).
% 
%  T    IN: Transition options. The first column contains the rates, the
%           second/third the incoming/outcoming states (M,3)
%  L    IN: Division options. The first column contains the rates, the
%           second the incoming states and the third and fourth the outcoming states (L,4)
%  N    IN: Number of states
%  p0  OUT: constant component (N,1)
%  p1  OUT: linear component (N,N)
%  p2  OUT: quadratic component rearranged to have dimensions N, and N^2 (N,N^2)
%  rtot  OUT: total rate for each state (1,N)
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

% total rates
rtot = zeros(1,N);
for is = 1:N
   rtot(is) = sum(T(T(:,2)==is, 1)) + sum(L(L(:,2)==is, 1));
end

% initialization
p1 = zeros(N,N);
p2 = zeros(N,N^2);
% add transition (s^1)
for it = 1:size(T,1)
    indx = sub2ind([N N], T(it,2), T(it,3));
    p1(indx) = p1(indx) + T(it,1)/rtot(T(it,2));
end
% add splitting (s^2)
for is = 1:size(L,1)
	indx = sub2ind([N N], L(is,3), L(is,4));
    p2(L(is,2), indx) = p2(L(is,2), indx) + L(is,1)/rtot(L(is,2));
end
% p0
p0 = 1-p1*ones(N,1)-p2*ones(N*N,1);

return