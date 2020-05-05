% THIS FUNCTION COMPUTES VARIANCE DECOMPOSITIONS FROM A SVAR WHERE
% BcompMC is the companion form matrix of the reduced form VAR without the mean and/or trend
% VC_epsMC is the covariance matrix of the reduced form (not in companion form - length(VC_epsMC) = nvars = Nbig
% var_decomp_opt is a row vector defining which forecasts to compute var decomps for (i.e. var_decomp_opt = [1; 5; 10; 15; 20])
function [var_decomp] = get_var_decomp(Bcomp, VC_eps, var_decomp_opt);
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This routine computes forecast error variance decompostions
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% Bcomp:   matrix with structure of estimated reduced-form coefficients in 
%           companion form excluding all the deterministic terms. 
%           This includes the estimated parameters in the partition 
%           n x (n x (n x vlag)). The remaining partition is made up of a 
%           diagonalic identity matrix that pins down the lead-lag relation 
%           between variables. The size is (n x vlag) x (n x vlag).  
%
% VC_eps:  covariance matrix of reduced-form residuals, size n x n
%
% Notes: 
% A. VC_eps denotes the variance-covariance decomposition of shocks of
% reduced-form model;
% B. Nbig denotes the number of variables;
% C. VC_eps_chol denotes the Choleski decomposition of the
% variance-covariance matrix of shocks of reduced-form model;
% D. stdepsvec denotes the standard deviations of shocks of reduced-form
% model.
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% var_decomp:     forecast error variance decompositions num_var_decomp x nvars
%
%--------------------------------------------------------------------------
%
% Author:  Lee Smith June 7, 2013
%
%--------------------------------------------------------------------------
%%
Nbig = length(VC_eps);
Nbigcomp = length(Bcomp);
%%
%% 1.) COMPUTE THE CHOLESKY DECOMPOSITION OF THE REDUCED FORM COVARIANCE MATRIX AND COMPANION FORM
%%
VC_epschol = chol(VC_eps)';
%%
VC_epschol_comp = [VC_epschol; zeros(Nbigcomp-Nbig,Nbig)];
%%
%% 2.) NOW COMPUTE THE FORECAST ERRORS FOR THE DESIRED LENGTH
%%
max_forecast 	= max(var_decomp_opt);
forecast_error 	= zeros(Nbig,Nbig,max_forecast);
%%
for ii = 1:1:Nbig;
	for jj = 1:1:Nbig;
		forecast_error(ii,jj,1) = VC_epschol_comp(ii,jj);
	end;
end;
%%
for ff = 2:1:max_forecast;
	Bcomp_ff = (Bcomp^(ff-1));
		for ii = 1:1:Nbig;
			for jj = 1:1:Nbig;
			
				forecast_error(ii,jj,ff) = Bcomp_ff(ii,:)*VC_epschol_comp(:,jj) + forecast_error(ii,jj,ff-1);
				
			end;
		end;
end;
%%
%% 3.) NOW COMPUTE THE VARIANCE DECOMPOSITIONS OF THE FORECAST ERROR
%%
lvd = length(var_decomp_opt);
var_decomp = zeros(lvd,Nbig,Nbig);
%%
for ff = 1:1:lvd;
	forecast = var_decomp_opt(ff,1);
	for ii = 1:1:Nbig;
		for jj = 1:1:Nbig;
			FE_for_variable 	= forecast_error(ii,:,forecast);
			FEV_for_variable	= FE_for_variable*FE_for_variable';
			var_decomp(forecast,jj,ii) = (FE_for_variable(1,jj)^2)/FEV_for_variable;
		end; 	%% for jj
	end;		%% for ii
end;			%% for ff