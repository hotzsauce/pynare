"Uncertainty Shocks in a Model of Effective Demand"

By Susanto Basu and Brent Bundick

November 2016

README File

**All figures will be saved as .pdf files; they will not be displayed in Matlab. The final figures used in the paper can be found in the folders labeled "FinalFigures."**

To reproduce the baseline model from Sections 4-5 as well as the figures in Appendices B and D, run the plotfigures.m file in the BaselineModel folder.

To estimate the model with impulse response and moment matching, run the matchimpulseresponses.m file in the BaselineModel\Calibration\ImpulseResponseMomentMatching folder.  Figure 1 of the main text will appear in the folder after the estimation routine converges to a solution.

To replicate the calculations for the empirical and model-implied unconditional moments, run datamoments.m in the UnconditionalMoments\Data folder and simulatesmallsamplemoments.m in the UnconditionalMoments\Model folder.

To reproduce the figure from Section 6, run simulateestimatevarsmallsamples.m in the EstimateVARSimulatedDataMarkups folder.

The replication code for the empirical results in Section 2 is located in the VAREvidence folder.

The baseline empirical results from Section 2 can be reproduced by running the estimatebaselinevar.m file in the Quarterly\Baseline folder.  To replicate the structural uncertainty shocks from the VAR (shown in Appendix A), run the plotvxoshocks.m file in the same folder.

The empirical results in the Appendix A can be replicated by running one of the estimate*.m files in the Robustness folder, where the * will represent the particular code being run.

Replicating the zero lower bound results from Section 7 requires a Fortran compiler with the IMSL numerical libraries.
The decision rules are initialized and written to a file by Matlab.
Then, the Fortran code solves the model and re-writes the decision rules to a file.
The model simulations are then completed in Matlab.

To compile the model:

$F90 $F90FLAGS $LINK_FNL -o basubundickuncertaintyeffectivedemand basubundickuncertaintyeffectivedemand.f90 -O4

To run the model:

./basubundickuncertaintyeffectivedemand

After solving the model, run plotzlbfigures.m to replicate the zero lower bound results.
