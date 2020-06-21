clear all

% uncomment one of the following
% outFolder = 'out_cons_20190409150503_p30 - Copy';
% outFolder = 'ppro_out_cons_20190409150503 - Copy';
outFolder = 'ppro_out_ncons_20190409150503 - Copy';

% Load files
P1 = load(fullfile(outFolder, 'partI'));
P2 = load(fullfile(outFolder, 'partII'));

% merge output
xout = [P1.xout; P2.xout];

% save file
save(fullfile(outFolder, 'out_ppro'), 'xout')