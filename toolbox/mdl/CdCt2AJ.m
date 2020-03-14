function [A, J] = CdCt2AJ(Cd, Ct)
%  [A, J] = CdCt2AJ(Cd, Ct)
% This function buld the adiacency and jacoian matrices starting from the
% division and transition matrices.
% Assumptions:
% - in Ct all the state transitions
% - in Cd only cell division, ii --> jj+jj is also possible
% - only one cell division option per state
% - each element is the rate
% 
% Cd    IN: Division matrix (N,N)
% Ct    IN: Transition matrix (N,N)
% A     OUT: Adjacency (N,N)
% J     OUT: Jacobian matrix (N,N)
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

% Adjacency matrix
A = Cd + Ct;

% Jacobian
D = diag(sum(Ct) + sum(Cd/2));
J = A-D;

end