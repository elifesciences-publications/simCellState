function hfig = plotEnvelopeSim(xint, yint, N, xfilt, xaxis, iplot, yprof)

% This function prepare the figures with the envelope, 90%ile range, and
% the selected and reference profiles.
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

switch iplot
    case 'dist_norm'
        hfig = figure; setFigureProp(hfig, 10); % SI
    otherwise
        hfig = figure; setFigureProp(hfig, 20); % main text
end
hold on; grid on
% adjust colormap
map0 = colormap(hfig,'gray'); smap = size(map0, 1);
map0(1:round(smap*0.1),:) = []; smap = size(map0, 1); % less darker shade
xmap = unique([linspace(-1,0,smap) linspace(0,1,smap)]);
ymap = [flip(map0,1); map0(2:end,:)];
cmap = interp1(xmap,ymap,linspace(-1, 1, 101)');
colormap(hfig,cmap);

% plot
[M,c] = contourf(xint, yint, N, 0:0.1:100);
c.LineStyle = 'none';
hcb = colorbar;
set(get(hcb,'Title'),'String','%ile')
[cmin, hmin] = contour(xint, yint, N, xfilt(1)*[1 1]);
hmin.LineColor = 'k'; hmin.LineWidth = 2;
[cMax, hMax] = contour(xint, yint, N, xfilt(2)*[1 1]);
hMax.LineColor = 'k'; hMax.LineWidth = 2;
caxis([0 100])
set(gca, 'Layer', 'top')
axis(xaxis)
xlegAllNtw = '5-95%ile';
href = []; xlegref = [];
switch iplot
    case 'dist_norm'
        href = plot(xint, normpdf(xint, 1, sqrt(1/30)), 'r', 'LineWidth', 2);
        xlegref = 'Normal(1, 1/30))';
        xlab = '$x$'; ylab = '$P(x)$';
        xpos = 'northeast';
    case 'dist_normSigma'
        href = plot(xint, normpdf(xint, 0, 1), 'r', 'LineWidth', 2);
        xlegref = 'Normal(0, 1)';
        xlab = '$\tilde x$'; ylab = '$P(\tilde x)$';
        if xaxis(3) > 0
            set(gca, 'YScale', 'log')
        end
        xpos = 'northeast';
    case 'dist_exp'
        href = plot(xint, exppdf(xint, 1), 'r', 'LineWidth', 2);
        xlegref = 'Exp(1)';
        xlab = '$x$'; ylab = '$P(x)$';
        if xaxis(3) > 0
            set(gca, 'YScale', 'log')
        end
        xpos = 'northeast';
    case 'mean'
        hleg = hMax; xleg = {xlegAllNtw};
        xlab = 't/$\tau$'; ylab = '$\bar{n}_s$';
        xpos = 'northwest';
end
for ii = 1:length(yprof)
    hprof = plot(yprof{ii}(1,:), yprof{ii}(2,:), 'color', 'b', 'LineWidth', 1.5); % handle overwritten
end
xlabel(xlab); ylabel(ylab)
legend([hMax hprof href], [xlegAllNtw {'Selected models'} xlegref], 'Location', xpos)

end