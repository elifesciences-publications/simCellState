clear all
close all

% INPUT PARAMETERS
% ----------------
icase = 'out_ncons_20190409150503'; Nrun = 5e4;
scase = {Nrun};
simCase = 'genericNetwork';
simInFolder = 'GENERIC';
outDir = pwd;

% SD0 = 123189; rng(SD0);
SD0 = 2*round(cputime*1e3)+1; rng(SD0) % for repeatibility
jr = rand(1,1000);
ir = 225;
simCellStateClone(icase, scase, simInFolder, simCase, outDir, ir, jr(ir:end));

save(fullfile(outDir, 'io', 'OUT', simInFolder, icase, 'out'), 'SD0')