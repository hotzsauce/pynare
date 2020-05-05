Peter N. Ireland
"A New Keynesian Perspective on the Great Recession"
Journal of Money, Credit, and Banking MS#10-167

These files contain notes, data, and MATLAB programs that will allow you to reproduce the econometric work from my paper.

notes.pdf - Contains notes on solving and estimating the model.

data.xlsx - Contains the data in spreadsheet form.

grp.dat - Contains the data in text form, which can be read by MATLAB.

solve.m - Solves the model.

imp.m - Produces impulse responses using the solution found by solve.m.

est.m - Estimates the model.

bsse.m - Estimates the model and computes bootstrapped standard errors.

llfn.m - Contains the log likelihood function that gets maximizes by est.m and bsse.m.

tstar.mat - Vector of parameter estimates produced by bsse.m.

tstarbs.mat - Matrix of bootstrapped parameter estimates produced by bsse.m.

sevec.mat - Vector of standard errors produced by bsse.m.

vardec.m - Uses the parameter estimates found by est.m to compute variance decompositions.

ksmooth.m - Uses the parameter estimates found by est.m to construct smoothed estimates of the model's unobserved shocks.

countersim.m - Runs counterfactual simulations.

More detailed descriptions of each MATLAB program are given in the comments lines at the top of each file.


