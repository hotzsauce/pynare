% THIS FUNCTION COMPUTES CORRELATIONS AND STD. ERROS OF VARIABLES FROM A VAR WHERE
% Bcomp is the companion form matrix of the reduced form VAR without the mean and/or trend
% VC_eps is the covariance matrix of the reduced form (not in companion form - length(VC_eps) = nvars = Nbig
function [moments] = get_moments(Bcomp,VC_eps);
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
VC_eps_comp			= VC_eps_comp/1000;
vec_VC_eps_comp 	= reshape(VC_eps_comp,Nbigcomp^2,1);
long_run_MAT		= eye(Nbigcomp^2) - kron(Bcomp,Bcomp);
vec_SIGMA_Z			= long_run_MAT\vec_VC_eps_comp;
SIGMA_Z				= reshape(vec_SIGMA_Z,Nbigcomp,Nbigcomp);
%%
%% 2.) NOW COMPUTE THE CORRELATIONS AND STANDARD DEVIATIONS
%%
moments = zeros(Nbig,Nbig);
%%
for ii = 1:1:Nbig;
	
	for jj = ii:1:Nbig;
		if ii == jj;
			moments(ii,jj) = (SIGMA_Z(ii,jj)^.5);
		else
			moments(ii,jj) = SIGMA_Z(ii,jj)/((SIGMA_Z(ii,ii)*SIGMA_Z(jj,jj))^.5);
		end;
	end;
	
end;