function dydt = varDer2(t, y)

% This function computes the second order finite difference of y(t).
% Extreme points are modified (first derivative is applied to the middle point)
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

n=length(y);
dydt = zeros(1, n);
% dydt(1)=(y(2)-y(1))/(t(2)-t(1));
for i=2:n-1
    dydt(i)=(y(i+1)-y(i-1))/(t(i+1)-t(i-1));
end
% dydt(n)=(y(n)-y(n-1))/(t(n)-t(n-1));

dydt(1) = varDerFirstLast(t(1:3), y(1:3));
dydt(n) = varDerFirstLast(t(end-2:end), y(end-2:end));

end

function dydt = varDerFirstLast(t, y)

dydt = diff(y)./diff(t); 
dydt = interp1(t(1:2)+diff(t)/2, dydt, t(1), 'linear', 'extrap');

end