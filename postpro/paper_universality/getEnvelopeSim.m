function [xint, yint, N, Nfilt, prcOut] = getEnvelopeSim(xout, xint, minOut, maxOut, xfilt, iextr)

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

% remove those with just one bin
sxout = cellfun(@(x) size(x, 2), xout); xout(sxout == 1) = []; 
% interpolation
if iextr == 1
    yy = cell2mat(cellfun(@(y) interp1([-10 y(1,:) 1e3], [y(2,1) y(2,:) y(2,end)], xint), xout, 'UniformOutput', false));
else
    yy = cell2mat(cellfun(@(y) interp1([y(1,:) 1e3], [y(2,:) y(2,end)], xint), xout, 'UniformOutput', false)); % only at higher values
end
% saturation
if ~isempty(minOut)
    yy(yy<minOut) = minOut;
else
    minOut = min(min(yy));
end
if ~isempty(maxOut)
    yy(yy>maxOut) = maxOut;
else
    maxOut = max(max(yy));
end

% find percentile
Nout = length(xout);
Ny = max([round(Nout/2) 500]); dy = (maxOut-minOut)/Ny;
yEdges = minOut-dy/2:dy:maxOut+dy/2;
yint = yEdges(1:end-1)+diff(yEdges)/2;
N = NaN(length(yEdges)-1, length(xint));
Nfilt = NaN(length(yEdges)-1, length(xint));
Xout = zeros(Nout, length(xint));
for ii = 1:length(xint)
    indxNotNan = find(~isnan(yy(:,ii)));
    nn = histcounts(yy(indxNotNan,ii), yEdges, 'normalization', 'cumcount');
    if length(indxNotNan) > 1
        % normalize count
        nn = nn'/length(indxNotNan)*100;
        N(:,ii) = nn;
        % cut below min and above 100
        imin = find(nn<xfilt(1), 1, 'last');
        if imin > 1
            nn(1:imin-1) = NaN;
            ymin = yint(imin);
        else
            ymin = -Inf;
        end
        imax = find(nn>xfilt(2), 1, 'first');
        if imax < length(nn)
            nn(imax+1:end) = NaN;
            ymax = yint(imax);
        else
            ymax = Inf;
        end
        % matrix
        Nfilt(:,ii) = nn;
        % indentify outliers
        indxOut = find(yy(:,ii) < ymin | yy(:,ii) > ymax);
        if length(indxOut) < Nout
            Xout(indxOut,ii) = 1;
        else
            Xout(indxOut,ii) = NaN;
        end
    end
end

% percentage of points outside the 90%ile range
prcOut = sum(Xout,2, 'omitnan')./sum(~isnan(Xout),2)*100;

end