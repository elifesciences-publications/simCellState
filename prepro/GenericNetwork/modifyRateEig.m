function [Cd1, Ct1, egMax1] = modifyRateEig(Cd1, Ct1, eigTg, Rmin)
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

% parameters/initialization
Ns = size(Cd1, 1);
inode = 1:Ns-1;

% identification of indexes and bounds
indxCd = find(sum(Cd1));
indxCt = find(Ct1);
Nvar = length(indxCd) + length(indxCt);
dr0 = zeros(Nvar,1);
sCd = sum(Cd1>0); gam = ones(size(sCd)); gam(sCd~=0) = 2;
LB = Rmin-[sum(Cd1(:,indxCd))./gam(indxCd) Ct1(indxCt)']; % Rmin = R+LB
UB = 10*ones(1,Nvar);

% fmincon
options = optimoptions(@fmincon, 'display', 'off');
[dr,fval,exFlag] = fmincon(@(x) stabilitySCCper(Cd1, Ct1, inode, indxCd, indxCt, x, eigTg), dr0, [], [], [], [], ...
                    LB, UB, [], options);

% output
if ~isempty(indxCd)
    dCd = zeros(size(Cd1)); 
    for im = 1:length(indxCd)
        dCd(Cd1(:,indxCd(im))~=0, indxCd(im)) = dr(im);
    end
    Cd1 = Cd1 + dCd; dr(1:length(indxCd)) = [];
end
if ~isempty(indxCt)
    Ct1(indxCt) = Ct1(indxCt)+dr(1:length(indxCt));
end
[~, J] = CdCt2AJ(Cd1, Ct1);
eg = eig(J(1:Ns-1,1:Ns-1)); [~, im] = max(real(eg)); 
egMax1 = eg(im);

% % check here
% if min(min(Ct1(Ct1>0))) < 0.3
%     disp('Check min Ct');
% end
% if any(sum(Cd1(:,sCd>0)) < 0.6)
%     disp('Check min Cd');
% end

return

% % try first just upregulating outgoing or division rates
% deltaeg = egMax1-eigTg;
% if deltaeg > 0 % need to increase outgouing rates
%     if any(Ct1(end,:)>0)
%         indxCt = sub2ind([Ns Ns], Ns, find(Ct1(end,:)>0));
%     end
%     if any(Cd0(end,:)==2) % before changing rate
%         indxCd = find(Cd0(end,:)==2); % index of column
%     end
%     ipar = length(indxCd) + length(indxCt);
%     LB = 0*ones(1,ipar);
%     UB = 1.5*ones(1,ipar);
% else % need to increase division rates
%     Cd0s = Cd0; Cd0s(end,:) = 0;
%     if any(any(Cd0s==2)) % before changing rate
%         [~, indxCd] = ind2sub([Ns Ns], find(Cd0s));
%     end
%     ipar = length(indxCd) + length(indxCt);
%     LB = 0.3-[sum(Cd1(:,indxCd))/2 Ct1(indxCt)];
%     UB = 1.5*ones(1,ipar);
% end

% options = optimoptions(@ga, 'display', 'off');
% [x,fval] = ga(@(x) stabilitySCCper(Cd1, Ct1, inode, indxCd, indxCt, x, eigTg), ipar, [], [], [], [], ...
%                     LB, UB, [], options)