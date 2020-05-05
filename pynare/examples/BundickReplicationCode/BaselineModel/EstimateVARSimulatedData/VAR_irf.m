
% EDITED BY LEE SMITH 5/20/2013 TO CHANGE LABELING OF GRAPHS

function [IR] = VAR_irf(Bcomp, VC_eps, IRhoriz, IRtype, sizesho, ispl);


%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This routine computes and plots the impulse response function
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
% ispl:    1 to plot the impulse responses; 0 for no plotting 
%
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% IR:      impulse responses (in level)
%
%--------------------------------------------------------------------------
%
% Author:  Paolo Z., September 2011
%
%--------------------------------------------------------------------------



Nbig       = length(VC_eps);
Nbigcomp   = length(Bcomp);

VC_epschol = chol(VC_eps)'; 
stdepsvec  = diag(VC_eps).^0.5;

if IRtype=='c';
    if sizesho==0; 
        shockvecmat = VC_epschol*eye(Nbig);
    else;
        shockvecmat = VC_epschol*diag(sizesho./stdepsvec);
    end
    
elseif IRtype=='g'; 
    
    if sizesho==0; 
        shockvecmat = VC_eps*diag(stdepsvec.^(-1));
    else;
        shockvecmat = VC_epschol*diag(stdepsvec.^(-2))*diag(sizesho);
    end
end 


[nz,nofIR]         = size(shockvecmat);
aux                = zeros(Nbigcomp, Nbigcomp);
aux(1:Nbig,1:Nbig) = shockvecmat;
shockvecmat        = aux;

IR     = zeros(IRhoriz+1, Nbig, nofIR);
Impmat = eye(length(Bcomp));


for hh=1:IRhoriz+1;
    IRbig  = Impmat*shockvecmat;
    IR(hh, 1:Nbig, 1:Nbig) = IRbig(1:Nbig, 1:Nbig);
    Impmat = Impmat*Bcomp;
end



if ispl==1;
    xaxis     = 0:IRhoriz;
    zerol     = zeros(IRhoriz+1, 1);
    plotcount = 1;
    
%% NOTE TO SELF:
%% ii loops through rows - top to bottom - (which are variables)
%% jj loops through columns - left to right - (which are shocks)

    for ii=1:Nbig;
        
        yaxmin = 1.1*min(min(IR(:, ii, :)));
        yaxmax = 1.1*max(max(IR(:, ii, :)));
        
        for jj=1:nofIR;
            
            subplot( Nbig, nofIR, plotcount );  
            plot( xaxis, IR(:,ii,jj), xaxis, zerol );  
            axis( [0 (IRhoriz) yaxmin yaxmax] ); 
            %title( ['shock' int2str(jj) ' to var ' int2str(ii)] );
            title( ['var' int2str(ii) ' to shock ' int2str(jj)] );
			plotcount = plotcount+1;
            
        end
    end

end

