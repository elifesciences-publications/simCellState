### This code is associated with the paper from Parigini and Greulich, "Universality of clonal dynamics poses fundamental limits to identify stem cell self-renewal strategies". eLife, 2020. http://dx.doi.org/10.7554/eLife.56532



% TOOL INSTALLATION
Run install.m

% 1) GENERIC NETWORK GENERATION
Run the following scripts under prepro folder:
 - srun_rawDT.m
 - srun_singleSCC.m
 - srun_poolSCC.m
 - sbuild_genericNetwork.m
This generates the random networks which are the input for the stocastic simulations. 
The current input data files (i.e. networks) used for running the simulations are currently saved under io\IN\GENERIC folder:
 - GIA networks: out_cons_20190409150503\out_cons_20190409150503.mat
 - GPA networks: out_ncons_20190409150503\out_ncons_20190409150503.mat

% 2) RUN AND POSTPROCESS STOCHASTIC SIMULATIONS
To run and postprocess simulations, use the following scripts under the main and postpro\GenericNetwork folders:
 - srunSimCellState.m
 - sprocess_genericSimPaper.m
Additional scripts:
 - srun_simGADModifiedRates.m - to update the network parameters and re-run the simulations
 - srun_asymmetricDivision.m, srun_populationAsymmetryR.m and srun_populationAsymmetryM.m - to run simulations for SI
Simulation output data files - postprocessed - are currently saved under io\OUT\GENERIC folder (if needed use smergeOutFiles.m to build the complete output files).

% 3) FIGURE
The scripts for producing the figures of the main text are:
 - splot_randSimPaper.m (Fig.2/3 and Fig.4 (a))
 - splot_pSensitivity.m (Fig.4 (a))

Figures for the SI
 - splot_testCases.m (Fig.1-4)
 - stest_GIA0.m (Fig.5-8)
 - stest_EquivalentModel.m (Fig.9-13)
 - splot_randSimPaper.m (Fig.14)
 - stestSimulationGIAB.m (Fig.15)
 - splot_randSimExpConvPaper.m (Fig.16)
