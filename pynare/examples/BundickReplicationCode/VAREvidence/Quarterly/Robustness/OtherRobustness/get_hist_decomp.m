% THIS FUNCTION COMPUTES THE HISTORICAL DECOMPOSITION FROM A SVAR WHERE
% BcompMC is the companion form matrix of the reduced form VAR without the mean and/or trend
% VC_epsMC is the covariance matrix of the reduced form (not in companion form - length(VC_epsMC) = nvars = Nbig
% Residmat is a (Tbig-vlag) x Nbig matrix of the reduced form residuls from the estimated VAR
% hist_decomp_opt is a 2 x Nbig matrix defining the following:
% hist_decomp_opt(1,1) is an integer which:
%	 if <= vlag computes the forecast beginning at the earliest possible period
%	 elseif Tbig > > vlag computes the forecast from period hist_decomp_opt(1,1)
% hist_decomp_opt(1,2) equals either 1 or 0 with:
%	 if == 1 then subtract out base forecast (forecast with no shocks)
%	 if == 0 then DON'T subtract out base forecast
% hist_decomp_opt(1,3) equals either 1 or 0 with:
%	 if == 1 then plot the decomps
%	 if == 0 then DON'T plot the decomps
% hist_decomp_opt(1,4) equals either 1 or 0 with:
%	 if == 1 then DON'T compute hist_decomps with confidence bands in confid_int
%	 if == 0 then compute hist_decomps with confidence bands in confid_int
% hist_decomp_opt(2,j) equals either 1 or 0 for j  1:Nbig with:
%	 if equals 1 then sum these decomps
% 	 if equals 0 then don't sum (so if all are 0 there will be no combined decomps)
function [hist_decomp,start_date,struct_shocks] = get_hist_decomp(Bcomp, VC_eps, Y, cvec, hist_decomp_opt)
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This routine computes historical of structural shocks decompostions
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
% C. VC_epschol denotes the Choleski decomposition of the
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
check = 0;
%%
Nbig = length(VC_eps);
Nbigcomp = length(Bcomp);
vlag = Nbigcomp/Nbig;
Tbig = length(Y(:,1));
resid_mm = Tbig - vlag;
resid_nn = Nbig;
Y_LHS_comp = zeros(Nbig*vlag,1);
%%
mm_cvec_comp = Nbig*(vlag-1);
cvec_comp = [cvec(1:Nbig,1); zeros(mm_cvec_comp,1)];
%%
for ll = 1:1:vlag;
	upper = ll*Nbig;
	lower = upper - Nbig + 1;
	Y_LHS_comp(lower:upper,1) = Y(vlag - ll + 1,:)';	
end;
%%
for tt = 1:1:resid_mm;
Y_RHS_comp = Y_LHS_comp;
	for ll = 1:1:vlag;
		upper = ll*Nbig;
		lower = upper - Nbig + 1;
		Y_LHS_comp(lower:upper,1) = Y(tt + 1 + vlag - ll,:)';
	end;
Resid_comp = Y_LHS_comp - Bcomp*Y_RHS_comp - cvec_comp;
Residmat(tt,:) = Resid_comp(1:Nbig,1)';
end;
%%

%%
Tbig = length(Residmat) + vlag; 
%% DEFINE NUMBER OF ROWS
if hist_decomp_opt(1,1) <= vlag;
	n_rows = Tbig - vlag;
	start_date = vlag+1;
	resid_start_date = start_date - vlag;
elseif hist_decomp_opt(1,1) > length(Residmat)
	error('Forecast Horizon Too Short');
else
	n_rows = Tbig - hist_decomp_opt(1,1);
	start_date = hist_decomp_opt(1,1)+1;
	resid_start_date = start_date + 1 - vlag;
end;
%% DEFINE NUMBER OF COLUMNS
combine_hist = hist_decomp_opt(2,:);
sum_combine_hist = sum(combine_hist);
if sum_combine_hist == 0;
	n_cols = Nbig;
else
	n_cols = Nbig + 1;
end;

hist_decomp = zeros(n_rows,n_cols+1,Nbig);

hist_decomp_comp = zeros(Nbig*vlag,n_rows);

forecast_horizon = size(hist_decomp,1);

%%
%%
%% 1.) COMPUTE THE CHOLESKY DECOMPOSITION OF THE REDUCED FORM COVARIANCE MATRIX AND GENERATE THE STRUCTURAL RESIDUALS
%%
VC_epschol = chol(VC_eps)';
inv_VC_epschol = inv(VC_epschol);
%%
VC_epschol_comp = [VC_epschol; zeros(Nbigcomp-Nbig,Nbig)];
%%
%% IF e_t (Nbig x 1) represents the reduced form VAR residual and eps_t denote the structural VAR residual  at time t,
%% then with the recursive identification scheme imposed by the Choleski decomposition, we have that at each time t
%% e_t = VC_epschol*eps_t. So we can recover eps_t from e_t using the equation eps_t = inv(VC_epschol)*e_t. 
struct_shocks = zeros(Tbig-vlag,Nbig);

for t = 1:1:(Tbig-vlag);
	struct_shocks(t,:) = (inv_VC_epschol*Residmat(t,:)')';
end;

%% EACH Nbig x Nbig matrix (for all time periods) decomposes the contribution of the structural shocks
%% to the reduced form VAR residuals. The first column is the first structural shock mapped backed into
%% reduced form residual, the second column ... etc. Notice, the sum of the columns at any time t will 
%% generate a Nbig x 1 vector that equals the reduced form residual at time t.
struct_shocks_contribution = zeros(Nbigcomp,Nbig,(Tbig-vlag)); 

for t = 1:1:(Tbig-vlag);
	for ii = 1:1:Nbig;
		struct_shocks_contribution(1:Nbig,ii,t) = VC_epschol(:,ii)*struct_shocks(t,ii);
	end;
end;

%%-----CHECK--------
if check == 1;
	sum_shocks = zeros(Nbig,(Tbig-vlag));
		for tt = 2:1:(Tbig-vlag)
			for jj = 1:1:Nbig;
				sum_shocks(:,tt) = sum_shocks(:,tt) +struct_shocks_contribution(1:Nbig,jj,tt);
			end;
		end;
sum_shocks == Residmat'
end;
%%---------------------
Y_init = Y;
[Tbig,Nbig] = size(Y_init);
%for ii = 1:1:Nbig;
%	mean_vec = ones(Tbig,1)*cvec(ii);
%	Y_init(:,ii) = Y_init(:,ii) - mean_vec;
%end;

for ll = 1:1:vlag;
	upper = ll*Nbig;
	lower = upper - Nbig + 1;
	hist_decomp_comp(lower:upper,1) = Y_init(start_date+1-ll,:)';	
end;	

Y_init(start_date,:);
hist_decomp_comp(:,1);

%% 2.) NOW COMPUTE THE FORECAST (SIMULATE THE MODEL) WITH THE STRUCTURAL SHOCKS
%%
for jj = 1:1:Nbig
	for tt = 2:1:forecast_horizon;
	resid_date = resid_start_date + tt - 2;
	hist_decomp_comp(:,tt) = cvec_comp + Bcomp*hist_decomp_comp(:,tt-1) + struct_shocks_contribution(:,jj,resid_date);
		for ii = 1:1:Nbig;
			hist_decomp(tt,jj,ii) = hist_decomp_comp(ii,tt);
		end;
	end;
end;
%%
if n_cols > Nbig;
	struct_shocks_contribution_sum = zeros(Nbigcomp,1,(Tbig-vlag)); 
	for ii = 1:1:Nbig;
		if hist_decomp_opt(2,ii) == 1;
		struct_shocks_contribution_sum = struct_shocks_contribution_sum + struct_shocks_contribution(:,ii,:);
		end;
	end;
	%%
	for tt = 2:1:forecast_horizon;
	resid_date = resid_start_date + tt - 2;
	hist_decomp_comp(:,tt) = cvec_comp + Bcomp*hist_decomp_comp(:,tt-1) + struct_shocks_contribution_sum(:,:,resid_date);
		for ii = 1:1:Nbig;
			hist_decomp(tt,Nbig+1,ii) = hist_decomp_comp(ii,tt);
		end;
	end;
end;
%%
for tt = 2:1:forecast_horizon;
resid_date = resid_start_date + tt - 2;
hist_decomp_comp(:,tt) = cvec_comp + Bcomp*hist_decomp_comp(:,tt-1);
	for ii = 1:1:Nbig;
		hist_decomp(tt,n_cols+1,ii) = hist_decomp_comp(ii,tt);
	end;
end;

%% 3.) SUBTRACT THE BASE FORECAST IF DESIRED
if hist_decomp_opt(1,2) == 1;
	
	hist_decomp_comp_base = zeros(Nbig*vlag,n_rows);
	hist_decomp_comp_base(:,1) = hist_decomp_comp(:,1);
	
	%%
	for ii = 1:1:Nbig;
			hist_decomp(1,n_cols+1,ii) = hist_decomp_comp_base(ii,1); 
		for jj = 1:1:n_cols;
			hist_decomp(1,jj,ii) = hist_decomp(1,jj,ii) - hist_decomp_comp_base(ii,1);
		end;
	end;
	%%
	
	for tt = 2:1:forecast_horizon;
	hist_decomp_comp_base(:,tt) = cvec_comp + Bcomp*hist_decomp_comp_base(:,tt-1);
		for ii = 1:1:Nbig;
			hist_decomp(tt,n_cols+1,ii) = hist_decomp_comp_base(ii,tt); 
			for jj = 1:1:n_cols;
				hist_decomp(tt,jj,ii) = hist_decomp(tt,jj,ii) - hist_decomp_comp_base(ii,tt);  
			end;
		end;
	end;	
	
end;
hist_decomp = hist_decomp(2:end,:,:);
%% 3.) PLOT IF DESIRED
if hist_decomp_opt(1,3) == 1;
%	for ii = 1:1:Nbig;
%		mean_vec = ones(horizon,1)*mean(ii);
%		Y_plot(:,ii) = Y_plot(:,ii) - mean_vec;
%	end;
%Y_plot(1,:)
Y_plot = Y(start_date+1:end,:);
[horizon,Nbig] = size(Y_plot);
n_shocks = max(Nbig,n_cols);
xaxis     = 1:horizon;
zerol     = zeros(horizon, 1);
plotcount = 1;   
%%
%% COMPUTE SUM OF SQUARED DEVIATIONS IF DESIRED
ssr = zeros(Nbig,1);
shock_check = 4;
resid_hist = [];
for ii = 1:1:Nbig;
	resid_hist = [resid_hist; (Y_plot(:,ii) - hist_decomp(:,shock_check,ii))]/100;
	%ssr(ii,1) = resid_hist'*resid_hist;
end;
ssr = resid_hist'*resid_hist;
ssr
%%
figure;

  
%% NOTE TO SELF:
%% ii loops through rows - top to bottom - (which are variables)
%% jj loops through columns - left to right - (which are shocks)

    for ii=1:Nbig;
        for jj=1:n_shocks;
		size(Y_plot);
		size(hist_decomp);
		plot_vector = [Y_plot(:,ii), hist_decomp(:,jj,ii)];
        yaxmin = 1.1*min(plot_vector);
        yaxmax = 1.1*max(plot_vector);
			%
			%
            subplot( Nbig, n_shocks, plotcount );  
            plot( xaxis, hist_decomp(:,jj,ii)', '-b', xaxis, Y_plot(:,ii)', '-r', xaxis, zerol, '-k', 'LineWidth', 1.5 );  
            %axis( [0 (horizon-1) yaxmin yaxmax] ); 
            %title( ['shock' int2str(jj) ' to var ' int2str(ii)] );
            title( ['var' int2str(ii) ' to shock ' int2str(jj)] );
			xlim([1 horizon]);
			plotcount = plotcount+1;
            
        end
    end


end;
%%-----------------------------