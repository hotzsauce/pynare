% THIS FUNCTION COMPUTES AUTO-CORREATIONS FROM A VAR WHERE
% Bcomp is the companion form matrix of the reduced form VAR without the mean and/or trend
% VC_eps is the covariance matrix of the reduced form (not in companion form - length(VC_eps) = nvars = Nbig
% num_auto_corr is the number of auto-correlations to compute for each variable
function [auto_corr] = get_auto_corr(Bcomp,VC_eps,num_auto_corr);
%%
Nbig = length(VC_eps);
Nbigcomp = length(Bcomp);
%%
%% 1.) COMPUTE THE LONG RUN COVARIANCE MATRIX
%%
VC_eps_comp = [VC_eps, zeros(Nbig,(Nbigcomp-Nbig)); zeros(Nbigcomp-Nbig,Nbigcomp)];
%%
%% SIGMA_Z = Bcomp*SIGMA_Z*Bcomp' + VC_eps_comp
%% USE VEC PROPERTY (cf. Hamilton, 1994, pg. 3048 footnote)
%% vec(A*B*C) = kron(C',A)*vec(B)
%% vec(SIGMA_Z) 		= kron(Bcomp,Bcomp)*vec(SIGMA_Z) 				+ vec(VC_eps_comp)
%% (Nbigcomp^2 x 1) = (Nbigcomp^2 x Nbigcomp^2)*(Nbigcomp^2 x 1)	+ (Nbigcomp^2 x 1)
%%
%% vec(SIGMA_Z)		= [I - kron(Bcomp,Bcomp)]\vec(VC_eps_comp); 
%% SIGMA_Z = reshape(vec(SIGMA_Z),Nbigcomp,Nbigcomp);
%%
vec_VC_eps_comp 	= reshape(VC_eps_comp,Nbigcomp^2,1);
long_run_MAT		= eye(Nbigcomp^2,Nbigcomp^2) - kron(Bcomp,Bcomp);
vec_SIGMA_Z			= long_run_MAT\vec_VC_eps_comp;
SIGMA_Z				= reshape(vec_SIGMA_Z,Nbigcomp,Nbigcomp);
%%
%% 2.) NOW COMPUTE THE AUTO-COVARIANCES
%%
auto_covariances = zeros(num_auto_corr+1,Nbigcomp);
%%
for ii = 1:1:(num_auto_corr+1);
	auto_covariances(ii,:) = diag((Bcomp^(ii-1))*SIGMA_Z)';
	%% IF YOU WANTED CROSS-AUTO-CORRELATIONS KEEP THE WHOLE MATRIX INSEAD OF THE DIAGONAL
end;
%%
for ii = 1:1:num_auto_corr;
	for jj = 1:1:Nbig;
		auto_corr(ii,jj) = auto_covariances(ii+1,jj)/auto_covariances(1,jj);
	end;
end;