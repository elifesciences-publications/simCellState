function save2eps(figPath)

% This function savesa copy of all the fig-files in eps (color)
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

dd = dir(fullfile(figPath, '*.fig'));

for ii = 1:length(dd)
    hh = openfig(fullfile(figPath, dd(ii).name));
    [~,name,~] = fileparts(dd(ii).name);
    saveas(hh, fullfile(figPath, name), 'epsc')
    close(hh)
end