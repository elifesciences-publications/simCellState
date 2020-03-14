function [eigA, eigJ] = eigSCC(A, J, S, C)
%  [eigA, eigJ] = eigSCC(A, J, S, C)
% This function returns the eigenvalues of the adiacency and jacobian
% matrices for each strongly connected component.
%
% A     IN: Adjacency (N,N)
% J     IN: Jacobian (N,N)
% S     IN: number of strongly connected components
% C     IN: id of the strongly connected components for each state (1,N)
% eigA  OUT: eigenvalues of A for each strongly connected component, cell(1,S)
% eigJ  OUT: eigenvalues of J for each strongly connected component, cell(1,S)
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

% eigenvalues
eigA = cell(1, S);
eigJ = cell(1, S);
for isi = 1:S
    idinode = find(C == isi);
    eigA{isi} = eig(A(idinode, idinode));
    eigJ{isi} = eig(J(idinode, idinode));
end

end