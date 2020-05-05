% imp.m: Uses the solution found by solve.m to compute impulse responses
% for the New Keynesian model with habit formation in preferences and
% partial indexation in price setting.
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

% IS shock

sa = zeros(52,8);
fa = zeros(52,6);

evec = [ 100*siga 0 0 0 ]';

sa(2,:) = (bigw*evec)';
fa(2,:) = (bigu*sa(2,:)')';

for t = 3:52

    sa(t,:) = (bigpi*sa(t-1,:)')';
    fa(t,:) = (bigu*sa(t,:)')';

end
  
% cost-push shock

se = zeros(52,8);
fe = zeros(52,6);

evec = [ 0 100*sige 0 0 ]';

se(2,:) = (bigw*evec)';
fe(2,:) = (bigu*se(2,:)')';

for t = 3:52

    se(t,:) = (bigpi*se(t-1,:)')';
    fe(t,:) = (bigu*se(t,:)')';

end
  
% technology shock

sz = zeros(52,8);
fz = zeros(52,6);

evec = [ 0 0 100*sigz 0 ]';

sz(2,:) = (bigw*evec)';
fz(2,:) = (bigu*sz(2,:)')';

for t = 3:52

    sz(t,:) = (bigpi*sz(t-1,:)')';
    fz(t,:) = (bigu*sz(t,:)')';

end
  
% monetary policy shock

sr = zeros(52,8);
fr = zeros(52,6);
  
evec = [ 0 0 0 100*sigr ]';

sr(2,:) = (bigw*evec)';
fr(2,:) = (bigu*sr(2,:)')';

for t = 3:52

    sr(t,:) = (bigpi*sr(t-1,:)')';
    fr(t,:) = (bigu*sr(t,:)')';

end  

% create output vectors

ya = sa(2:52,1);
pa = sa(2:52,2);
ra = sa(2:52,3);
qa = sa(2:52,4);
aa = sa(1:51,5);
ea = sa(1:51,6);
za = sa(1:51,7);
  
xa = fa(1:51,1);
ga = fa(1:51,2);
lama = fa(1:51,3);
  
bigya = cumsum(ga);
  
ye = se(2:52,1);
pe = se(2:52,2);
re = se(2:52,3);
qe = se(2:52,4);
ae = se(1:51,5);
ee = se(1:51,6);
ze = se(1:51,7);
  
xe = fe(1:51,1);
ge = fe(1:51,2);
lame = fe(1:51,3);
  
bigye = cumsum(ge);
  
yz = sz(2:52,1);
pz = sz(2:52,2);
rz = sz(2:52,3);
qz = sz(2:52,4);
az = sz(1:51,5);
ez = sz(1:51,6);
zz = sz(1:51,7);
  
xz = fz(1:51,1);
gz = fz(1:51,2);
lamz = fz(1:51,3);
  
bigyz = cumsum(gz);
  
yr = sr(2:52,1);
pr = sr(2:52,2);
rr = sr(2:52,3);
qr = sr(2:52,4);
ar = sr(1:51,5);
er = sr(1:51,6);
zr = sr(1:51,7);
  
xr = fr(1:51,1);
gr = fr(1:51,2);  
lamr = fr(1:51,3);
  
bigyr = cumsum(gr);