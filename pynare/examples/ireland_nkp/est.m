% est.m: Maximizes (minimizes -1 times) the log likelihood function for the
% New Keynesian model with habit formation in preferences, and partial
% indexation in price setting.
%
% When maximizing the log likelihood function, the parameters are
% transformed to satisfy theoretical restrictions. The log likelihood
% function with transformed parameters in contained in llfn.m.
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