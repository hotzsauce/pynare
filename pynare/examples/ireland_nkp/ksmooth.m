% ksmooth.m: Uses the Kalman smoother, together with the untransformed
% estimates found by est.m, to produce a series of smoothed estimates of
% the unobserved state variables from the New Keynesian model with habit
% formation in preferences and partial indexation in price setting.
%
% The results are contained in the 8xT matrix sbigtvec, which keeps track
% of the 8 variables in the model's state equation. Additional results are
% contained in the 6xT matrix fbigtvec, which keeps track of the 6
% variables in the model's flow equation, and the 4xT matrix ebigtvec,
% which keeps track of the model's 4 innovations. Finally, the covariance
% and correlation matrices for the innovations are reported as ecovs and
% ecors.
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
  
% form matrices AX, BX, CX, VX, and BVBX

bigax = bigpi;

bigbx = bigw;

bigcx = [ bigu(2,:) ; bigu(5,:) ; bigpi(3,:) ];

bigvx = diag([siga^2 sige^2 sigz^2 sigr^2]);

bigbvbx = bigbx*bigvx*bigbx';

% assemble data vector

dthat = [ gt pit rt ];
  
capt = length(gt);

% run through the Kalman Filter
  
sttvec = zeros(8,capt);
  
stlvec = zeros(8,capt);
  
sigttvec = zeros(64,capt);
 
sigtlvec = zeros(64,capt);

st = zeros(8,1);

bigsigt = (eye(64)-kron(bigax,bigax))\bigbvbx(:);
  
bigsigt = reshape(bigsigt,8,8);
  
for t = 1:capt
      
    stlvec(:,t) = st;
    
    sigtlvec(:,t) = bigsigt(:);
      
    omegt = bigcx*bigsigt*bigcx';
    
    ut = dthat(t,:)' - bigcx*st;
    
    stt = st + bigsigt*bigcx'/omegt*ut;
    
    bigsigtt = bigsigt - bigsigt*bigcx'/omegt*bigcx*bigsigt;
    
    sttvec(:,t) = stt;
    
    sigttvec(:,t) = bigsigtt(:);
    
    st = bigax*stt;
    
    bigsigt = bigbvbx + bigax*bigsigtt*bigax';
    
end  

% run through Kalman smoother

sbigtvec = zeros(8,capt);
  
sbigtvec(:,capt) = sttvec(:,capt);
  
for j = 1:capt-1
      
    bigsigtt = sigttvec(:,capt-j);
    
    bigsigtt = reshape(bigsigtt,8,8);
    
    bigsigtp = sigtlvec(:,capt-j+1);
    
    bigsigtp = reshape(bigsigtp,8,8);
    
    bigjt = bigsigtt*bigax'*pinv(bigsigtp);
    
    stt = sttvec(:,capt-j);
    
    stpbigt = sbigtvec(:,capt-j+1);
    
    stpt = stlvec(:,capt-j+1);
    
    sbigtvec(:,capt-j)= stt + bigjt*(stpbigt-stpt);
    
end
  
% construct smoothed flow and innovation vectors

fbigtvec = bigu*sbigtvec;
  
ebigtvec = zeros(4,capt);
  
ebigtvec(:,1) = sbigtvec(5:8,1);
  
ebigtvec(:,2:capt) = sbigtvec(5:8,2:capt) - bigp*sbigtvec(5:8,1:capt-1);
                     
ecovs = cov(ebigtvec');
  
ecors = corrcoef(ebigtvec');