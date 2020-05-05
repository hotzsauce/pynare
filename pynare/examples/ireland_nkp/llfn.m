function llfn = llfn(bigthet)
% Uses the Kalman filter to evaluate the negative log likelihood function
% for the New Keynesian model with habit formation in preferences and
% partial indexation in price setting. The parameters are transformed to
% satisfy theoretical restrictions.
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

% define variables and parameters

global gt pit rt z beta

capt = length(gt);
  
bigthet = real(bigthet);
  
gammatr = bigthet(1);
alphatr = bigthet(2);
psitr = 0.10;
  
rhortr = 1;
rhoptr = bigthet(3);
rhoxtr = 0;
rhogtr = bigthet(4);

rhoatr = bigthet(5);
rhoetr = bigthet(6);
  
sigatr = bigthet(7);
sigetr = bigthet(8);
sigztr = bigthet(9);
sigrtr = bigthet(10);
  
% untransform parameters

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
  
lamviol = 0;

if abs(bigt11(4,4)) > abs(bigs11(4,4))
    lamviol = 1;
end

if abs(bigt22(1,1)) < abs(bigs22(1,1))
    lamviol = 1;
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

% form matrices AX, BX, CX, VX, and BVBX

bigax = bigpi;

bigbx = bigw;

bigcx = [ bigu(2,:) ; bigu(5,:) ; bigpi(3,:) ];

bigvx = diag([siga^2 sige^2 sigz^2 sigr^2]);

bigbvbx = bigbx*bigvx*bigbx';

% assemble data vector

dthat = [ gt pit rt ];

% evaluate negative log likelihood

st = zeros(8,1);

bigsig1 = (eye(64)-kron(bigax,bigax))\bigbvbx(:);

bigsigt = reshape(bigsig1,8,8);

llfn = (3*capt/2)*log(2*pi);

for t = 1:capt

    ut = dthat(t,:)' - bigcx*st;

    omegt = bigcx*bigsigt*bigcx';
    
    if rcond(omegt) < eps
        
        llfn = 1e12;
       
    else

        llfn = llfn + (1/2)*(log(det(omegt))+ut'/omegt*ut);

        bigkt = bigax*bigsigt*bigcx'/omegt;

        st = bigax*st + bigkt*ut;

        bigsigt = bigbvbx + bigax*bigsigt*bigax' ...
                    - bigax*bigsigt*bigcx'/omegt*bigcx*bigsigt*bigax';

    end
    
end

% penalize eigenvalue constraint violations

if lamviol

    llfn = llfn + 1e12;

end

if imag(llfn) > 0
      
    llfn = llfn + 1e12;
    
end