%--------------------------------------------------------------------------
%
% Original Author:  Paolo Z., September 2011
% EDITED BY LEE SMITH 5/20/2013 TO CHANGE LABELING OF GRAPHS
% EDITED BY LEE SMITH 6/4/2013  TO:
% ADD (OPTIONAL) COMPUTATION OF: AUTO-CORRLATIONS, MOMENTS, AND VARIANCE DECOMPOSITIONS
% BAYESIAN COMPUTATION OF POSTERIOR CONFIDENCE BANDS
%--------------------------------------------------------------------------


function [IRpoint, IRup_2, IRup_1, IRlo_1, IRlo_2, IRmed, IRstd, Bevmean, Bevstd, auto_corr, moments, var_decomp] = ...
	confid_int(Bcomp,VC_eps,cvec,dvec,IRhoriz,IRtype,sizesho,NMC,Tbig,CIperc,ispl,num_auto_corr,no_moments,var_decomp_opt,hist_decomp_opt,Y,plot_options,save_graphs,save_marker,var_labels,shock_labels);

%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This routine computes the impulse response functions, their confidence 
% intervals and associated statistics
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
% VC_eps:  covariance matrix of reduced-form residuals, size n x n.
%
% cvec:    matrix with estimated constants; this returns the scalar 0 
%           if no constant is included in the model; 
%           otherwise, this is a matrix of size n x 1.
%
% dvec:    matrix with estimated parameters on the deterministic trends;
%           this returns the scalar 0 if no constant is included 
%           otherwise, this is a matrix of size n x 1.
%
% IRhoriz: number of periods for which point impulse responses are computed
%
%
% IRtype:  this string variable can be assigned the values'c' or 'g'
%
% sizesho: this is a singleton that is used to normalize the Choleski 
%           transformation of the variance-covariance matrix of 
%           reduced-form residuals. If it has a value equal to 0, 
%           alternative normalizations of the Choleski are used,
%           according to the following table
%
%
% Table 1: structural shock vector as a function of the two key inputs
% ---------------------------------------------------------------------------------------------------------
% |            |                                                                                          |
% |            |              IRtype=c                                    IRtype=g                        |
% |            |------------------------------------------------------------------------------------------|
% |            |                                                                                          |
% | sizesho=0  |         VC_eps_chol*eye(Nbig)                    VC_eps*diag(stdepsvec.^(-1))            |
% |            |                                                                                          |
% |sizesho~=0  |   VC_eps_chol*diag(sizesho./stdepsvec)   VC_eps_chol*diag(stdepsvec.^(-2))*diag(sizesho) |
% |            |                                                                                          |
% |            |                                                                                          |
% ---------------------------------------------------------------------------------------------------------
% Notes: 
% A. VC_eps denotes the variance-covariance decomposition of shocks of
% reduced-form model;
% B. Nbig denotes the number of variables;
% C. VC_eps_chol denotes the Choleski decomposition of the
% variance-covariance matrix of shocks of reduced-form model;
% D. stdepsvec denotes the standard deviations of shocks of reduced-form
% model.
%
%
% NMC:     number of Monte Carlo simulations for the density of simulated
%           impulse responses
%
% Tbig:    number of simulation periods after burn-in
%
% CIperc:  this is a 2 x 1 vector with lower and upper percentiles for 
%           the confidence itervals that are computed and, eventually,
%           plotted
%
% ispl:    1 to plot the impulse responses; 0 for no plotting 
%
% burnin:  number of observations used for the burn-in/initialization of 
%           the estimation of the VAR models
%
% num_auto_corr: number of autocorrelations to compute
%
% no_moments: 1 to NOT compute any variances and correlations; 0 otherwise computes these
%
% var_decomp_opt: if rows of var_decomp_opt > 1 then compute variance decomp; otherwise DON'T compute 
%
% var_labels: nvars x 1 vector of strings with variable names that print on subplot (DEFAULT IS NO LABELS)
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% IRpoint: point estimate of impulse responses 
%
% IRup:    upper bound of impulse responses 
%
% IRlo:    lower bound of impulse responses
%
% IRmed:   median of impulse response distribution (at each horizon point)
%
% IRstd:  standard deviation of impulse response distribution (at each horizon point)
%
% Bevmean: mean of companion matrix coefficients over the Monte Carlo runs
%
% Bevstd:  standard deviation of companion matrix coefficients over 
%           the Monte Carlo runs
%
% auto_corr: 	IF num_auto_corr > 0: 
%				num_auto_corr x nvars x 3 MATRIX with: 
%			 	auto_corr(:,:,1) MATRIX is the LOWER confidence bound
%				auto_corr(:,:,2) MATRIX is the MEDIAN confidence bound
%				auto_corr(:,:,3) MATRIX is the UPPER confidence bound
%
% moments:		IF no_moments == 0: 
%				nvars x nvars x 3 UPPER DIAGONAL MATRIX with:
%				moments(:,:,1) MATRIX is the LOWER confidence bound
%				moments(:,:,2) MATRIX is the MEDIAN confidence bound
%				moments(:,:,3) MATRIX is the UPPER confidence bound
%
% var_decomp    IF length(var_decomp_opt) > 1: COMPUTES THE VARIANCE DECOMPOSITIONS AT THE MEAN COEFF ESTIMATES:
%				var_decomp = length(var_decomp_opt) x nvars x nvars with: 
%					var_decomp(:,:,1) MATRIX is the FEVD for the first variable
%					var_decomp(:,:,2) MATRIX is the FEVD for the second variable
%					.
%					.
%					.
%					var_decomp(:,:,nvars) MATRIX is the FEVD for the last variable
%
% hist_decomp	IF no_hist_decomp == 0:
%				forecast_length x Nshocks x Nbig x 3 with: (See the options of get_hist_decomp for more info on size of 
%				forecast_length and Nshocks
%					hist_decomp(:,:,:,1) MATRIX is the LOWER confidence bound
%					hist_decomp(:,:,:,2) MATRIX is the MEDIAN confidence bound
%					hist_decomp(:,:,:,3) MATRIX is the UPPER confidence bound
%-----------------------------------------------------------------------------------------------------------------------------------
%% CHECK IF THE USER DEFINED LABELS
if nargin == 19;
lvlab = 0;
lslab = 0;
h = 0;
elseif nargin == 20;
lvlab = length(var_labels);
lslab = 0;
elseif nargin == 21;
lvlab = length(var_labels);
lslab = length(shock_labels);
end;
%%----------------------------------------------------------------------------------------------------------------------------------
%% BACK OUT THE SIZE OF VARIOUS OBJECTS OF INTEREST
trendh = 1:Tbig;
Nbig   = length(VC_eps);
Nbigcomp = length(Bcomp);
vlag   = length(Bcomp)/Nbig;
BIG_X 		= zeros(Tbig-vlag,Nbig*vlag+1);
BIG_X(:,1) 	= ones(Tbig-vlag,1);
for lags = 1:1:vlag;
	lower_index_row		= vlag+1-lags;
	upper_index_row		= Tbig-lags;
	upper_index_col		= Nbig*lags+1;
	lower_index_col		= upper_index_col - Nbig + 1;
	BIG_X(:,lower_index_col:upper_index_col) = Y(lower_index_row:upper_index_row,:); 	
end;

BIG_Y = Y(vlag+1:Tbig,:);

istr  = sum(abs(dvec))>0;
iscon = sum(abs(cvec))>0;

Bevsum    = zeros(Nbig*vlag+iscon+istr, Nbig);
Bevsqsum  = zeros(Nbig*vlag+iscon+istr, Nbig);
VC_epssum = zeros(Nbig,Nbig);

allIRs 	 = zeros(IRhoriz+1, Nbig, Nbig, NMC);
IRup_1   = zeros(IRhoriz+1, Nbig, Nbig);
IRlo_1   = IRup_1;
IRup_2	 = IRup_1;
IRlo_2	 = IRup_1;
IRmed  	 = IRup_1;
IRmean   = IRup_1;

IRpoint=VAR_irf(Bcomp, VC_eps, IRhoriz, IRtype, sizesho, 0);


if istr==0;
    dvec = zeros(Nbig, 1);
end

if iscon==0;
    cvec = zeros(Nbig, 1);
end

if vlag>1;
    cvec = [cvec; zeros((vlag-1)*Nbig, 1)];
    dvec = [dvec; zeros((vlag-1)*Nbig, 1)];
end

%% MAKE SURE THE INPUTS MAKE SENSE
[po_rows,po_cols] = size(plot_options);
if po_rows ~= Nbig || po_cols ~= Nbig
	display('Plot Options Matrix Must be Nbig x Nbig -- Plotting Everything Instead');
	plot_options = ones(Nbig,Nbig);
end;
%%--------------------------------------------------------------
%% PRE-ALLOCATE FOR AUTO-CORRELATIONS (IF NEEDED)
l_CI = length(CIperc);
if num_auto_corr > 0;
	auto_corr = zeros(num_auto_corr, Nbig, l_CI+1);
	auto_corr_draws = zeros(num_auto_corr, Nbig, NMC);
	comp_auto_corr = 1;
else
	auto_corr = 0;
end
%% PRE-ALLOCATE FOR MOMENTS (IF NEEDED)
if no_moments == 0;
	moments = zeros(Nbig, Nbig, l_CI+1);
	moments_draws = zeros(Nbig, Nbig, NMC);
	comp_moments = 1;
else
	moments = 0;
end
%% PRE-ALLOCATE FOR VAR DECOMPS (IF NEEDED)
lvd = length(var_decomp_opt);
if lvd > 1;
	var_decomp = zeros(lvd, Nbig, Nbig);
else
	var_decomp = 0;
end
%% PRE-ALLOCATE FOR HIST DECOMPS (IF NEEDED)
no_hist_decomp = hist_decomp_opt(1,4);
if no_hist_decomp == 0;
	%% DEFINE NUMBER OF ROWS
	if hist_decomp_opt(1,1) <= vlag;
		n_rows = Tbig - vlag;
		start_date = vlag+1;
		resid_start_date = start_date - vlag;
	%elseif hist_decomp_opt(1,1) > length(Residmat)
	%	error('Forecast Horizon Too Short');
	else
		%n_rows = Tbig - vlag - hist_decomp_opt(1,1)
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
hist_decomp_draws = zeros(n_rows,n_cols,Nbig,NMC); 
hist_decomp = zeros(n_rows,n_cols,Nbig, l_CI + 1);
end;
%%--------------------------------------------------------------
Idmatcomp      	= eye(vlag*Nbig);
uncondmeancomp 	= inv(Idmatcomp-Bcomp)*cvec;
VC_epscompchol 	= zeros(Nbig*vlag, Nbig*vlag);
VC_epscompchol(1:Nbig, 1:Nbig) = chol(VC_eps)'; %diag(diag(chol(VC_eps)'));
BIG_A 			= (BIG_X'*BIG_X)\BIG_X'*BIG_Y;  % Nbig*vlag+1 x Nbig
mean			= BIG_A(1,:)';
vec_BIG_A 		= reshape(BIG_A,Nbig*(1+Nbig*vlag),1);
wishdof 		= Tbig;
scale			= Tbig-1;
COVAR 			= inv(BIG_X'*BIG_X);
%% DO MC RUNS
for mm=1:NMC;
	draw = mm;
	workbar(mm/NMC,'Performing Monte Carlo Integration...','Progress');
	%% HERE I EMPLOY ANTITHETIC SAMPLING WHICH SPEEDS UP THE CONVERGENCE
	if mod(mm,2) == 1
		% odd draw again
		draw_again = 0;
		while draw_again == 0;
			VC_eps_draw 		= iwishrnd(scale*VC_eps,wishdof);
			vec_BIG_A_COVAR		= kron(VC_eps_draw,COVAR);
			vec_A_draw_u		= mvnrnd(zeros(Nbig*(1+Nbig*vlag),1)',vec_BIG_A_COVAR);
			vec_BIG_A_draw		= vec_BIG_A + vec_A_draw_u'; 
			BIG_A_draw			= reshape(vec_BIG_A_draw,1+Nbig*vlag,Nbig);
			Bcomp_draw			=[BIG_A_draw(2:end,:)'; kron(eye(vlag-1), eye(Nbig)), ... 
									zeros((vlag-1)*Nbig, Nbig)];
			e_1					= eigs(Bcomp_draw);
			norm_1				= norm(e_1(1,1));
			vec_BIG_A_draw		= vec_BIG_A - vec_A_draw_u'; 
			BIG_A_draw			= reshape(vec_BIG_A_draw,1+Nbig*vlag,Nbig);
			Bcomp_draw			=[BIG_A_draw(2:end,:)'; kron(eye(vlag-1), eye(Nbig)), ... 
									zeros((vlag-1)*Nbig, Nbig)];
			e_2					= eigs(Bcomp_draw);
			norm_2				= norm(e_2(1,1));
			if norm_1 < 1 && norm_2 < 1;
				draw_again = 1;
			end;
			%%
		end;
		vec_BIG_A_draw		= vec_BIG_A + vec_A_draw_u'; 
		BIG_A_draw			= reshape(vec_BIG_A_draw,1+Nbig*vlag,Nbig);
	else
		% even draw - take reflection
		vec_BIG_A_draw		= vec_BIG_A - vec_A_draw_u'; 
		BIG_A_draw			= reshape(vec_BIG_A_draw,1+Nbig*vlag,Nbig);
	end;
	%%
	Bpl_ev_draw			= flipud(BIG_A_draw);
	Bevsum				= Bevsum + Bpl_ev_draw;
	Bevsqsum 			= Bevsqsum + Bpl_ev_draw.^2;
	VC_epssum 			= VC_epssum + VC_eps_draw;
	%%
	Bcomp_draw			=[BIG_A_draw(2:end,:)'; kron(eye(vlag-1), eye(Nbig)), ... 
									zeros((vlag-1)*Nbig, Nbig)];
	IRloc    			= VAR_irf(Bcomp_draw, VC_eps_draw, IRhoriz, IRtype, sizesho, 0);   
	allIRs(:, :, :, mm) = IRloc;
	[auro_corr_sample] = get_auto_corr(Bcomp_draw,VC_eps_draw,num_auto_corr);
	auto_corr_draws(:,:,mm) = auro_corr_sample;
	[moments_sample]	= get_moments(Bcomp_draw,VC_eps_draw);
	moments_draws(:,:,mm) = moments_sample;
end
%% COMPUTE VARIOUS STATISTICS
Bevmean = Bevsum/NMC;
Bevstd  = (Bevsqsum/NMC-Bevmean.^2).^0.5;
VC_epsmean = VC_epssum/NMC;
Bcompmean = [fliplr(Bevsum(1:Nbig*vlag,1:Nbig)'); kron(eye(vlag-1), eye(Nbig)), ... 
                                        zeros((vlag-1)*Nbig, Nbig)];
%% COMPUTE VARIANCE DECOMPOSTION IF DESIRED
if lvd > 1;
var_decomp = get_var_decomp(Bcomp,VC_eps,var_decomp_opt);
end;
%% GET LOWER, MEDIAN AND UPPER VARIABLES OF INTEREST FROM SIMULATED DISTRIBUTION
%% AUTOCORRELATIONS
if num_auto_corr > 0;
%%
for acs = 1:1:num_auto_corr;
	
	for jj = 1:1:Nbig;
		slice_auto_corr = squeeze(auto_corr_draws(acs,jj,:));
		if l_CI == 2;
		auto_corr_loc   = quantile(slice_auto_corr,[CIperc(1);0.5;CIperc(2)]);
		elseif l_CI == 4;
		auto_corr_loc   = quantile(slice_auto_corr,[CIperc(1);CIperc(2);0.5;CIperc(3);CIperc(4)]);
		end;
		auto_corr(acs,jj,:) = auto_corr_loc';
	end;

end;
%%
end;
%% MOMENTS
if no_moments == 0;
	%%
	for ii = 1:1:Nbig;
		for jj = 1:1:Nbig;
			slice_moments = squeeze(moments_draws(ii,jj,:));
			if l_CI == 2;
				moments_loc	  = quantile(slice_moments,[CIperc(1);0.5;CIperc(2)]);
			elseif l_CI == 4;
				moments_loc	  = quantile(slice_moments,[CIperc(1);CIperc(2);0.5;CIperc(3);CIperc(4)]);
			end;
			moments(ii,jj,:) = moments_loc';
		end;
	end;
	%%
end;
%% IRFS
%% CHECK SIGN
sign_check = plot_options == abs(plot_options);
%%
for hh=1:IRhoriz+1;
    
    for ii=1:Nbig;
        
        for jj=1:Nbig;
			
			flip = all(sign_check(:,jj));
			if flip == 0;
				multiplier = -1;
			elseif flip == 1;
				multiplier = 1;
			end;
			
            sliceloc   =  multiplier*squeeze(allIRs(hh,ii,jj,:));
			
			if l_CI == 2;
            IRpercsloc = quantile(sliceloc,[CIperc(1);0.5;CIperc(2)]);
            IRlo_2(hh,ii,jj)   = IRpercsloc(1);
			IRmed(hh,ii,jj)    = IRpercsloc(2);
            IRup_2(hh,ii,jj)   = IRpercsloc(3);
			elseif l_CI == 4;
			IRpercsloc = quantile(sliceloc,[CIperc(1);CIperc(2);0.5;CIperc(3);CIperc(4)]);
            IRlo_1(hh,ii,jj)   = IRpercsloc(1);
			IRlo_2(hh,ii,jj)   = IRpercsloc(2);
            IRmed(hh,ii,jj)    = IRpercsloc(3);
			IRup_2(hh,ii,jj)   = IRpercsloc(4);
			IRup_1(hh,ii,jj)   = IRpercsloc(5);
			end;
			IRstd(hh,ii,jj)    = std(sliceloc);
			%IRmeanloc  = mean(sliceloc);
			%IRmean(hh,ii,jj) = IRmeanloc;
        end
        
    end
    
end
%%--------------------
plot_options = abs(plot_options);
%% PLOT IF DESIRED - LABELING AND PLOTTING FUNCTIONS SUPPORTED FOR (UP TO) 6 VARIABLE VAR
if lvlab >= 1 && lslab == 0;
	variables = var_labels;
	shocks    = var_labels;
elseif lvlab >= 1 && lslab >= 1;
	variables = var_labels;
	shocks    = shock_labels;
else
	variables.a = '1';
	variables.b = '2';
	variables.c = '3';
	variables.d = '4';
	variables.e = '5';
	variables.f = '6';
	shocks = variables;
end;
%%
vn = fieldnames(variables);
sn = fieldnames(variables);
%%
if ispl == 0;
	h = 0;
elseif ispl==1;
%%
plot_all = plot_options == ones(Nbig,Nbig);
	if plot_all == ones(Nbig,Nbig);
		plot_all_shocks(l_CI,IRhoriz,IRmed,IRlo_2,IRup_2,variables,shocks,vn,sn,IRlo_1,IRup_1,save_graphs,save_marker);
	else
		for jj = 1:1:Nbig;
			if sum(plot_options(:,jj)) > 0;
				plot_one_shock(l_CI,IRhoriz,plot_options,jj,IRmed,IRlo_2,IRup_2,variables,shocks,vn,sn,IRlo_1,IRup_1,save_graphs,save_marker);
			end;
		end;
	end;
end; % if/elseif ispl