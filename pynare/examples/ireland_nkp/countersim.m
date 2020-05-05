% countersim.m: Constructs counterfactual paths for output, inflation, and
% interest rates during the 1991, 2001, and 2008 recessions using the
% maximum likelihood estimates and the smoothed estimates of the shocks
% from the New Keynesian model with habit formation in preferences and
% partial indexation in price setting.
%
% Produces as output for xx=91, xx=01, and xx=08
% 
%   yacxx = actual path for output in levels
%   pacxx = actual path for inflation in annualized terms
%   racxx = actual path for the interest rate in annualized terms
%
%   ycfaxx = counterfactual path for output with preference shocks only
%   ycfexx = counterfactual path for output with cost-push shocks only
%   ycfzxx = counterfactual path for output with technology shocks only
%   ycfrxx = counterfactual path for output with monetary policy shocks only
%
%   pcfaxx = counterfactual path for inflation with preference shocks only
%   pcfexx = counterfactual path for inflation with cost-push shocks only
%   pcfzxx = counterfactual path for inflation with technology shocks only
%   pcfrxx = counterfactual path for inflation with monetary policy shocks only
%
%   pcfaxx = counterfactual path for the interest rate with preference shocks only
%   pcfexx = counterfactual path for the interest rate with cost-push shocks only
%   pcfzxx = counterfactual path for the interest rate with technology shocks only
%   pcfrxx = counterfactual path for the interest rate with monetary policy shocks only
%
% Also produces as output for xx=08
%
%   ycfxrxx = counterfactual path for output without monetary policy shocks
%   pcfxrxx = counterfactual path for inflation without monetary policy shocks
%   rcfxrxx = counterfactual path for the interest rate without monetary policy shocks
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
  
% construct counterfactual paths for 1991 recession

yac91 = [ 0 ; cumsum(gt(31:39)+log(z)) ];
pac91 = 4*(pit(30:39)+log(piss));
rac91 = 4*(rt(30:39)+log(z)-log(beta)+log(piss));

ebigtvec91 = ebigtvec;
  
ebigtvec91(2:4,31:108) = zeros(3,78);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec91(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfa91 = [ 0 ; cumsum(cfftvec(2,31:39)'+log(z)) ];
pcfa91 = 4*(cfftvec(5,30:39)'+log(piss));
rcfa91 = 4*(cfrtvec(30:39)'+log(z)-log(beta)+log(piss));
  
ebigtvec91 = ebigtvec;
  
ebigtvec91(1,31:108) = zeros(1,78);
ebigtvec91(3:4,31:108) = zeros(2,78);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec91(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfe91 = [ 0 ; cumsum(cfftvec(2,31:39)'+log(z)) ];
pcfe91 = 4*(cfftvec(5,30:39)'+log(piss));
rcfe91 = 4*(cfrtvec(30:39)'+log(z)-log(beta)+log(piss));  
  
ebigtvec91 = ebigtvec;
 
ebigtvec91(1:2,31:108) = zeros(2,78);
ebigtvec91(4,31:108) = zeros(1,78);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec91(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfz91 = [ 0 ; cumsum(cfftvec(2,31:39)'+log(z)) ];
pcfz91 = 4*(cfftvec(5,30:39)'+log(piss));
rcfz91 = 4*(cfrtvec(30:39)'+log(z)-log(beta)+log(piss)); 
  
ebigtvec91 = ebigtvec;
  
ebigtvec91(1:3,31:108) = zeros(3,78);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec91(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfr91 = [ 0 ; cumsum(cfftvec(2,31:39)'+log(z)) ];
pcfr91 = 4*(cfftvec(5,30:39)'+log(piss));
rcfr91 = 4*(cfrtvec(30:39)'+log(z)-log(beta)+log(piss)); 
  
% construct counterfactual paths for 2001 recession

yac01 = [ 0 ; cumsum(gt(73:81)+log(z)) ];
pac01 = 4*(pit(72:81)+log(piss));
rac01 = 4*(rt(72:81)+log(z)-log(beta)+log(piss));

ebigtvec01 = ebigtvec;
  
ebigtvec01(2:4,73:108) = zeros(3,36);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec01(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfa01 = [ 0 ; cumsum(cfftvec(2,73:81)'+log(z)) ];
pcfa01 = 4*(cfftvec(5,72:81)'+log(piss));
rcfa01 = 4*(cfrtvec(72:81)'+log(z)-log(beta)+log(piss));
  
ebigtvec01 = ebigtvec;
  
ebigtvec01(1,73:108) = zeros(1,36);
ebigtvec01(3:4,73:108) = zeros(2,36);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec01(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfe01 = [ 0 ; cumsum(cfftvec(2,73:81)'+log(z)) ];
pcfe01 = 4*(cfftvec(5,72:81)'+log(piss));
rcfe01 = 4*(cfrtvec(72:81)'+log(z)-log(beta)+log(piss));  
  
ebigtvec01 = ebigtvec;
  
ebigtvec01(1:2,73:108) = zeros(2,36);
ebigtvec01(4,73:108) = zeros(1,36);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec01(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfz01 = [ 0 ; cumsum(cfftvec(2,73:81)'+log(z)) ];
pcfz01 = 4*(cfftvec(5,72:81)'+log(piss));
rcfz01 = 4*(cfrtvec(72:81)'+log(z)-log(beta)+log(piss)); 
  
ebigtvec01 = ebigtvec;
  
ebigtvec01(1:3,73:108) = zeros(3,36);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec01(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfr01 = [ 0 ; cumsum(cfftvec(2,73:81)'+log(z)) ];
pcfr01 = 4*(cfftvec(5,72:81)'+log(piss));
rcfr01 = 4*(cfrtvec(72:81)'+log(z)-log(beta)+log(piss)); 
  
% construct counterfactual paths for 2008 recession

yac08 = [ 0 ; cumsum(gt(100:108)+log(z)) ];
pac08 = 4*(pit(99:108)+log(piss));
rac08 = 4*(rt(99:108)+log(z)-log(beta)+log(piss));

ebigtvec08 = ebigtvec;
  
ebigtvec08(2:4,100:108) = zeros(3,9);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec08(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfa08 = [ 0 ; cumsum(cfftvec(2,100:108)'+log(z)) ];
pcfa08 = 4*(cfftvec(5,99:108)'+log(piss));
rcfa08 = 4*(cfrtvec(99:108)'+log(z)-log(beta)+log(piss));
  
ebigtvec08 = ebigtvec;

ebigtvec08(1,100:108) = zeros(1,9);
ebigtvec08(3:4,100:108) = zeros(2,9);

cfstvec = zeros(8,capt);
 
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec08(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfe08 = [ 0 ; cumsum(cfftvec(2,100:108)'+log(z)) ];
pcfe08 = 4*(cfftvec(5,99:108)'+log(piss));
rcfe08 = 4*(cfrtvec(99:108)'+log(z)-log(beta)+log(piss));
  
ebigtvec08 = ebigtvec;
  
ebigtvec08(1:2,100:108) = zeros(2,9);
ebigtvec08(4,100:108) = zeros(1,9);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec08(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfz08 = [ 0 ; cumsum(cfftvec(2,100:108)'+log(z)) ];
pcfz08 = 4*(cfftvec(5,99:108)'+log(piss));
rcfz08 = 4*(cfrtvec(99:108)'+log(z)-log(beta)+log(piss));
  
ebigtvec08 = ebigtvec;
  
ebigtvec08(1:3,100:108) = zeros(3,9);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec08(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfr08 = [ 0 ; cumsum(cfftvec(2,100:108)'+log(z)) ];
pcfr08 = 4*(cfftvec(5,99:108)'+log(piss));
rcfr08 = 4*(cfrtvec(99:108)'+log(z)-log(beta)+log(piss));
  
ebigtvec08 = ebigtvec;
  
ebigtvec08(4,100:108) = zeros(1,9);

cfstvec = zeros(8,capt);
  
cfftvec = zeros(6,capt);
  
cfrtvec = zeros(1,capt);
  
st = sbigtvec(:,1);
  
ft = bigu*st;
  
cfstvec(:,1) = st';
    
cfftvec(:,1) = ft';
  
cfrtvec(:,1) = (bigpi(3,:)*st)';

for t = 2:capt
      
    et = ebigtvec08(:,t);
    
    st = bigpi*st + bigw*et;
    
    ft = bigu*st;
    
    cfstvec(:,t) = st';
    
    cfftvec(:,t) = ft';
    
    cfrtvec(:,t) = (bigpi(3,:)*st)';
    
end
  
ycfxr08 = [ 0 ; cumsum(cfftvec(2,100:108)'+log(z)) ];
pcfxr08 = 4*(cfftvec(5,99:108)'+log(piss));
rcfxr08 = 4*(cfrtvec(99:108)'+log(z)-log(beta)+log(piss));
 