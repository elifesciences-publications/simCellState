function [indexS, indexL] = findSourceLeak(Cd, Ct, inode)
%  [indexS, indexL] = findSourceLeak(Cd, Ct, inode)
% This function identify the source and leak nodes inside a strongly
% connected component of a given system.
% Assumptions:
% - Source nodes are those that potentially increase the number of cells.
%   (i.e. cell division where both children remain within the SCC).
% - Leak nodes are those that potentially decrease the number of cells.
%   (i.e. nodes forcing a cell to go outside the SCC and decreasing the
%   global number of cell inside the SCC. This means that the case of
%   division where one cell remains inside the SCC and one go outside is
%   not a leak; if both children are going outside it is a leak).
%  Cd    IN: Division matrix (N,N)
%  Ct    IN: Transition matrix (N,N)
%  inode IN: index of the strongly connected component nodes
%  indexS OUT: source nodes
%  indexL OUT: leak node
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

% deterministic division within the SCC
ss = sum(Cd~=0);
indexS = inode(sum(Cd(inode, inode)>0) == ss(inode) & ss(inode)>0); % two division inside the scc

% number of outgoing nodes
Cd1 = Cd;
Cd1(:,inode(sum(Cd(inode, inode)>0) == 1)) = 0; % remove division of the type: one inside and one outside
A1 = Cd1 + Ct; % just transition or in case all the divisions are going outside
indexL = inode(sum(A1(setdiff(1:1:length(A1),inode),inode), 1) > 0);

end