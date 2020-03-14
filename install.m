function install()
% path to be added
simPath = {'kernel' 'toolbox\mdl'  'toolbox\sta'  'toolbox\plot'};
for ii = 1:length(simPath)
    addpath(fullfile(pwd, simPath{ii}));
end
