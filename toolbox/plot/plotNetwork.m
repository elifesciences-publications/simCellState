function h = plotNetwork(sysProp, xtitle, varargin)
% h = plotNetwork(sysProp, xtitle, varargin)
% This function plot the input network. Strongly connected components are
% highlighted with the same color (jet).
% 
%  sysProp   IN: input system property structure (see detSysProperties)
%  xtitle    IN: title strings
%  varargin  IN: optional inputs for plot digraph, ex. {'layout', 'layered'}
%  h        OUT: figure handle

% system parameters
G = sysProp.digraph;
N = sysProp.N;
S = sysProp.NSSC;
C = sysProp.SCC;

% plot graph
if ismember(varargin(1:2:end), 'EdgeLabel')
    inPlot = varargin;
else
    inPlot = [{'EdgeLabel', round(G.Edges.Weight*1e3)/1e3} varargin];
end
% h = figure; hold on; grid off
h = [];
hg = plot(G, inPlot{:});
highlight(hg, 1:1:N, 'NodeColor', 'k')

% remove tick and axes
set(gca, 'XTick', [], 'YTick', [], 'Box', 'on')

% adjust colors
highlight(hg, G, 'EdgeColor', 'k','LineWidth',1.5)
colors = lines(7); colors = colors([1 7 3],:);% [0 0 1; 0 0.5 0; 1 0 0]; 
colorLeg = [0; 1; -1];
for is = 1:length(sysProp.sccIndx)
    idnode = sysProp.sccIndx{is};
    icol = sysProp.sccType(is)==colorLeg;
    highlight(hg, idnode, 'NodeColor', colors(icol,:))
    for jsi = 1:length(idnode)
        highlight(hg, idnode(jsi), intersect(idnode, successors(G, idnode(jsi))), 'EdgeColor', colors(icol,:))
    end
end

% adjust style of the edges
for iedge = 1:height(G.Edges)
    itr = sysProp.Ct(G.Edges.EndNodes(iedge, 2), G.Edges.EndNodes(iedge, 1));
    idv = sysProp.Cd(G.Edges.EndNodes(iedge, 2), G.Edges.EndNodes(iedge, 1));
    if itr && idv % transition+ division
        highlight(hg, G.Edges.EndNodes(iedge, 1), G.Edges.EndNodes(iedge, 2), 'LineStyle', '-.')
    elseif itr && ~idv % only transition
        highlight(hg, G.Edges.EndNodes(iedge, 1), G.Edges.EndNodes(iedge, 2), 'LineStyle', ':')
    end
end

% change committed and developing nodes
highlight(hg, find(sysProp.vertexType == 3), 'Marker', 's') % committed
highlight(hg, find(sysProp.vertexType == 1), 'Marker', 'v') % develop
highlight(hg, find(sysProp.vertexType == 0), 'Marker', 'x') % death
highlight(hg, find(sysProp.vertexType == -1), 'Marker', 'd') % other

% % add box with info
% if S > 0
%     msg = printStabilityBox(sysProp, S);
%     annotation('textbox',[0.7 0.8 0.01 0.01], 'String',msg,'FitBoxToText','on')
% end

% add title
if ~isempty(xtitle)
    title(xtitle)
end

end


% buildin function
% ----------------
function msg = printStabilityBox(sysProp, S)

if sysProp.stability.flag
    stbText = 'stable';
else
    eigTxt = [];
    if sysProp.stabFlags(1) == 1
        eigTxt = 'stable';
        if all(sysProp.stabFlags(2:3) == 0)
            eigTxt = [eigTxt ' (connected)'];
        end
    end
    if sysProp.stabFlags(2) == 1
        if isempty(eigTxt)
            eigTxt = 'growing';
        else
            eigTxt = [eigTxt '/growing'];
        end
    end
    if sysProp.stabFlags(3) == 1
        if isempty(eigTxt)
            eigTxt = 'decaying';
        else
            eigTxt = [eigTxt '/decaying'];
        end
    end
    stbText = ['unstable - ' eigTxt];
end
msg = sprintf('SCC: %d \n %s', S, stbText);
end

