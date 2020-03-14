function save2eps(figPath)

% This function savesa copy of all the fig-files in eps (color)

dd = dir(fullfile(figPath, '*.fig'));

for ii = 1:length(dd)
    hh = openfig(fullfile(figPath, dd(ii).name));
    [~,name,~] = fileparts(dd(ii).name);
    saveas(hh, fullfile(figPath, name), 'epsc')
    close(hh)
end