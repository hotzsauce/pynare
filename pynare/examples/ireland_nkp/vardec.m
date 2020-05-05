% vardec.m: Uses the parameter estimates found by est.m to compute forecast
% error variance decompositions for the New Keynesian model with habit
% formation in preferences, and partial indexation in price setting.
%
% Produces as output the 41x8 matrices varstot, varsa, varse, varsz, and
% varsr that report the forecast error variance in each element of the
% model's state vector at horizons k=1 through k=41 quarters ahead, in
% total and attributable to each of the model's four shocks.
%
% Produces as output the 41x8 matrices varsapct, varsepct, varszpct, and
% varsrpct that decomposre the forecasts error variances in the model's
% state variables into percentages due to each of the models' four shocks.
%
% Produces as output the 41x6 matrices varftot, varfa, varfe, varfz, and
% varfr that report the forecast error variance in each element of the
% model's flow vector at horizons k=1 through k=41 quarters ahead, in
% total and attributable to each of the model's four shocks.
%
% Produces as output the 41x6 matrices varfapct, varfepct, varfzpct, and
% varfrpct that decomposre the forecasts error variances in the model's
% flow variables into percentages due to each of the models' five shocks.
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
  
% form matrices AX, BX, and CX

bigax = bigpi;
  
bigbx = bigw;
  
bigcx = [ bigu(2,:) ; bigu(5,:) ; bigpi(3,:) ];
  
% calculate total k-step ahead forecast error variances
  
bigvx = diag([siga^2 sige^2 sigz^2 sigr^2]);

bigbvbx = bigbx*bigvx*bigbx';
  
bigaxk = eye(8);
  
siglitsk = zeros(8,8);
  
varstot = zeros(41,8);
  
varftot = zeros(41,6);
  
for k = 1:41
      
    siglitsk = siglitsk + bigaxk*bigbvbx*bigaxk';
    
    siglitfk = bigu*siglitsk*bigu';
    
    varstot(k,:) = diag(siglitsk)';
    
    varftot(k,:) = diag(siglitfk)';
    
    bigaxk = bigaxk*bigax;
      
end
  
% calculate k-step ahead forecast error variances due to preference shocks
  
bigvx = diag([siga^2 0 0 0]);

bigbvbx = bigbx*bigvx*bigbx';
  
bigaxk = eye(8);
  
siglitsk = zeros(8,8);
  
varsa = zeros(41,8);
  
varfa = zeros(41,6);
  
for k = 1:41
      
    siglitsk = siglitsk + bigaxk*bigbvbx*bigaxk';
    
    siglitfk = bigu*siglitsk*bigu';
    
    varsa(k,:) = diag(siglitsk)';
    
    varfa(k,:) = diag(siglitfk)';
    
    bigaxk = bigaxk*bigax;
      
end
  
% calculate k-step ahead forecast error variances due to cost push shocks
  
bigvx = diag([0 sige^2 0 0]);

bigbvbx = bigbx*bigvx*bigbx';
  
bigaxk = eye(8);
  
siglitsk = zeros(8,8);
  
varse = zeros(41,8);
  
varfe = zeros(41,6);
  
for k = 1:41
      
    siglitsk = siglitsk + bigaxk*bigbvbx*bigaxk';
    
    siglitfk = bigu*siglitsk*bigu';
    
    varse(k,:) = diag(siglitsk)';
    
    varfe(k,:) = diag(siglitfk)';
    
    bigaxk = bigaxk*bigax;
      
end
  
% calculate k-step ahead forecast error variances due to technology shocks
  
bigvx = diag([0 0 sigz^2 0]);

bigbvbx = bigbx*bigvx*bigbx';
  
bigaxk = eye(8);
  
siglitsk = zeros(8,8);
  
varsz = zeros(41,8);
  
varfz = zeros(41,6);
  
for k = 1:41
      
    siglitsk = siglitsk + bigaxk*bigbvbx*bigaxk';
    
    siglitfk = bigu*siglitsk*bigu';
    
    varsz(k,:) = diag(siglitsk)';
    
    varfz(k,:) = diag(siglitfk)';
    
    bigaxk = bigaxk*bigax;
      
end
  
% calculate k-step ahead forecast error variances due to monetary policy shocks
  
bigvx = diag([0 0 0 sigr^2]);

bigbvbx = bigbx*bigvx*bigbx';
  
bigaxk = eye(8);
  
siglitsk = zeros(8,8);
  
varsr = zeros(41,8);
  
varfr = zeros(41,6);
  
for k = 1:41
      
    siglitsk = siglitsk + bigaxk*bigbvbx*bigaxk';
    
    siglitfk = bigu*siglitsk*bigu';
    
    varsr(k,:) = diag(siglitsk)';
    
    varfr(k,:) = diag(siglitfk)';
    
    bigaxk = bigaxk*bigax;
      
end
  
% compute percentage due to each shock
  
varstotpct = varstot + (varstot==0);
varftotpct = varftot + (varftot==0);
    
varsapct = 100*varsa./varstotpct;
varsepct = 100*varse./varstotpct;
varszpct = 100*varsz./varstotpct;
varsrpct = 100*varsr./varstotpct;
    
varfapct = 100*varfa./varftotpct;
varfepct = 100*varfe./varftotpct;
varfzpct = 100*varfz./varftotpct;
varfrpct = 100*varfr./varftotpct;