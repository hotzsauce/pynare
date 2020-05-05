% bsse.m: Maximizes (minimizes -1 times the log likelihood function for the
% New Keynesian model with habit formation in preferences and partial
% indexaton in price setting and computes standard errors for the estimated
% parameters using a parametric bootstrapping algorithm.
%
% When maximizing the log likelihood function, the parameters are
% transformed to satisfy theoretical restrictions. The log likelihood
% function with transformed parameters in contained in llfn.m.
%
% Produces and saves as output:
%
%   tstar: vector of maximum likelihood estimates
%   sevec: vector of bootstrapped standard errors
%   tstarbs: matrix of boostrapped parameter estimates
%
% This program was written for MATLAB 7.10.0 (R2010a) by
%
% Peter Ireland
% Boston College
% Department of Economics
% 140 Commonwealth Avenue
% Chestnut Hill, MA 02467
% irelandp@bc.edu
% http://www2.bc.edu/~irelandp
%
% Copyright (c) 2010 by Peter Ireland. Redistribution is permitted for
% educational and research purposes, so long as no changes are made. All
% copies must be provided free of charge and must include this copyright
% notice.

% STEP 1: ESTIMATE THE MODEL

% load and demean data

global gt pit rt z beta

load gpr.dat;
  
z = exp(mean(gpr(:,1)));
piss = exp(mean(gpr(:,2)));
beta = (z*piss)/(exp(mean(gpr(:,3))));

gt = gpr(:,1) - log(z);
pit = gpr(:,2) - log(piss);
rt = gpr(:,3) - log(z) + log(beta) - log(piss);
  
% set starting values
  
gammatr = sqrt(0.50/0.50);
alphatr = sqrt(0.50/0.50);
psitr = 0.10;
  
rhortr = 1;
rhoptr = 0.25;
rhoxtr = 0;
rhogtr = 0.05;
  
rhoatr = sqrt(0.75/0.25);
rhoetr = sqrt(0.25/0.75);

sigatr = 0.01;
sigetr = 0.001;
sigztr = 0.01;
sigrtr = 0.0025;
  
bigtheto = [ gammatr alphatr ...
             rhoptr rhogtr ...
             rhoatr rhoetr ...
             sigatr sigetr sigztr sigrtr ]';
  
% maximize likelihood

options = optimset('Display','iter','LargeScale','off','MaxFunEvals',10000,'MaxIter',10000);

[ thetstar,llfnstar,exitflag ] = fminunc(@llfn,bigtheto,options);
  
% untransform estimates

thetstar = real(thetstar);

gammatr = thetstar(1);
alphatr = thetstar(2);
psitr = 0.10;
  
rhortr = 1;
rhoptr = thetstar(3);
rhoxtr = 0;
rhogtr = thetstar(4);

rhoatr = thetstar(5);
rhoetr = thetstar(6);
  
sigatr = thetstar(7);
sigetr = thetstar(8);
sigztr = thetstar(9);
sigrtr = thetstar(10);

gamma = gammatr^2/(1+gammatr^2);
alpha = alphatr^2/(1+alphatr^2);
psi = abs(psitr);
  
rhor = abs(rhortr);
rhop = abs(rhoptr);
rhox = abs(rhoxtr);
rhog = abs(rhogtr);
  
rhoa = rhoatr^2/(1+rhoatr^2);
rhoe = rhoetr^2/(1+rhoetr^2);
  
siga = abs(sigatr);
sige = abs(sigetr);
sigz = abs(sigztr);
sigr = abs(sigrtr);
  
tstar = [ gamma alpha rhop rhog rhoa rhoe ...
          siga sige sigz sigr ]';
        
% STEP 2: SOLVE THE ESTIMATED MODEL

% form matrices A, B, C, and P

biga = [ 0 0 0 0 0 0 0 beta*gamma*z 0 0 ; ...
         0 0 1 0 0 0 1 0 -1 0 ; ...
         0 0 0 0 0 0 0 0 beta 0 ; ...
         1 0 0 0 0 0 0 0 0 0 ; ...
         0 0 0 0 0 0 0 0 0 beta*gamma*z ; ...
         1 0 0 -1 0 0 0 0 0 0 ; ...
         0 0 1 0 0 0 0 0 0 0 ; ...
         1 0 0 0 0 0 0 0 0 0 ; ...
         0 1 0 0 0 0 0 0 0 0 ; ...
         0 0 0 1 0 0 0 0 0 0 ];
       
bigb = [ -gamma*z 0 0 0 0 0 (z-beta*gamma)*(z-gamma) z^2+beta*gamma^2 0 0 ; ....
         0 0 0 0 0 0 1 0 0 0 ; ...
         0 -alpha 0 0 0 0 psi 0 1+beta*alpha 0 ; ...
         1 0 0 0 0 1 0 0 0 0 ; ...
         0 0 0 -gamma*z 0 0 0 0 0 z^2+beta*gamma^2 ; ...
         0 0 0 0 1 0 0 0 0 0 ; ...
         0 0 rhor 0 rhox rhog 0 0 rhop 0 ; ...
         0 0 0 0 0 0 0 1 0 0 ; ...
         0 0 0 0 0 0 0 0 1 0 ; ...
         0 0 0 0 0 0 0 0 0 1 ];
           
bigc = [ -(z-beta*gamma*rhoa)*(z-gamma) 0 gamma*z  0 ; ...
         0 0 0 0 ; ...
         -psi -1 0 0 ; ...
         0 0 -1 0 ; ...
         -beta*gamma*(z-gamma)*(1-rhoa) 0 gamma*z 0 ; ...
         0 0 0 0 ; ...
         0 0 0 1 ; ...
         0 0 0 0 ; ...
         0 0 0 0 ; ...
         0 0 0 0 ];

bigp = [ rhoa 0 0 0 ; ...
         0 rhoe 0 0 ; ...
         0 0 0 0 ; ...
         0 0 0 0 ];

% form matrices Q, Z, S, and T  
  
[bigs,bigt,bigq,bigz] = qz(biga,bigb);
  
lamvec = ordeig(bigs,bigt);
  
[ lamvec2 lamind ] = sort(abs(lamvec));
  
lamord(lamind) = (1:10)';
  
[bigs,bigt,bigq,bigz] = ordqz(bigs,bigt,bigq,bigz,lamord);

bigq1 = bigq(1:4,:);
bigq2 = bigq(5:10,:);
  
bigz11 = bigz(1:4,1:4);
bigz12 = bigz(1:4,5:10);
bigz21 = bigz(5:10,1:4);
bigz22 = bigz(5:10,5:10);
  
bigs11 = bigs(1:4,1:4);
bigs12 = bigs(1:4,5:10);
bigs22 = bigs(5:10,5:10);

bigt11 = bigt(1:4,1:4);
bigt12 = bigt(1:4,5:10);
bigt22 = bigt(5:10,5:10);

if abs(bigt11(4,4)) > abs(bigs11(4,4))
    'error - no solution'
end

if abs(bigt22(1,1)) < abs(bigs22(1,1))
    'error - multiple solutions'
end

% form matrix R

bigra = bigs22/bigt22;
  
bigrb = bigq2*bigc;
  
vecr = (eye(24)-kron(bigp,bigra))\bigrb(:);

bigr = reshape(vecr,6,4);

% form matrices M1, M2, M3, and M4

bigm1 = bigz21/bigz11;
  
bigm2 = -(bigz22-bigz21/bigz11*bigz12)/bigt22*bigr;

bigm3 = bigz11/bigs11*bigt11/bigz11;
  
bigm4a = bigt11/bigz11*bigz12/bigt22*bigr ...
           + bigq1*bigc + bigs12/bigt22*bigr*bigp ...
           - bigt12/bigt22*bigr;
  
bigm4 = bigz11/bigs11*bigm4a - bigz12/bigt22*bigr*bigp;

% form matrices PI, W, and U

bigpi = [ bigm3 bigm4 ; zeros(4,4) bigp ];
  
bigpi = real(bigpi);

bigw = [ zeros(4,4) ; eye(4) ];
  
bigu = [ bigm1 bigm2 ];
  
bigu = real(bigu);
  
sigvec = [ siga ; sige ; sigz ; sigr ];
  
bigcx = [ bigu(2,:) ; bigu(5,:) ; bigpi(3,:) ];
  
% STEP 3: LOOP THROUGH THE BOOTSTRAP PROCEDURE

capt = length(gt);

tstarbs = zeros(1000,10);

for bsrep = 1:1000
      
    bsrep
      
% generate bootstrapped data

    gprbs = zeros(capt,3);

    stvec = zeros(8,1);  

    for t = 1:500
        
        etvec = sigvec.*randn(4,1);
      
        stvec = bigpi*stvec + bigw*etvec;
      
    end
    
    for t = 1:capt
        
        etvec = sigvec.*randn(4,1);
      
        stvec = bigpi*stvec + bigw*etvec;
      
        gprbs(t,:) = (bigcx*stvec)';
      
    end
    
    gt = gprbs(:,1);
    pit = gprbs(:,2);
    rt = gprbs(:,3);
    
% estimate model with bootstrapped data
 
    thetstbs = fminunc(@llfn,thetstar,options);
   
% untransform and save bootstrapped estimates

    thetstbs = real(thetstbs);
    
    gammatr = thetstbs(1);
    alphatr = thetstbs(2);
 
    rhoptr = thetstbs(3);
    rhogtr = thetstbs(4);

    rhoatr = thetstbs(5);
    rhoetr = thetstbs(6);
  
    sigatr = thetstbs(7);
    sigetr = thetstbs(8);
    sigztr = thetstbs(9);
    sigrtr = thetstbs(10);

    gamma = gammatr^2/(1+gammatr^2);
    alpha = alphatr^2/(1+alphatr^2);
  
    rhop = abs(rhoptr);
    rhog = abs(rhogtr);
  
    rhoa = rhoatr^2/(1+rhoatr^2);
    rhoe = rhoetr^2/(1+rhoetr^2);
  
    siga = abs(sigatr);
    sige = abs(sigetr);
    sigz = abs(sigztr);
    sigr = abs(sigrtr);
    
    tstarbs(bsrep,:) = [ gamma alpha rhop rhog rhoa rhoe ...
                         siga sige sigz sigr ]';
                     
end
  
% STEP 4: CALCULATE BOOTSTRAPPED STANDARD ERRORS

sevec = std(tstarbs)'; 

save tstar tstar

save sevec sevec

save tstarbs tstarbs