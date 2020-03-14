
% TOOL INSTALLATION
Run install.m

% GENERIC NETWORK GENERATION
Run the following scripts under prepro folder:
 - srun_rawDT.m
 - srun_singleSCC.m
 - srun_poolSCC.m
 - sbuild_genericNetwork.m

% RUN AND POSTPROCESS STOCHASTIC SIMULATIONS
Run the following script under the main and postpro/GenericNetwork folders:
 - srunSimCellState.m
 - srun_simGADModifiedRates.m
 - srun_simGADpsens.m

% FIGURE
Figures for the main text
 - splot_randSimPaper.m
 - splot_pSensitivity.m

Figures for the SI
 - splot_testCases.m
 - stest_GIA0.m
 - stest_EquivalentModel.m
 - splot_randSimExpConvPaper.m
 - splot_randSimExpConvPaper.m