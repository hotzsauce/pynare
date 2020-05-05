function [g,c] = hpfilter(y,lambda)
% Constructs the growth and cyclical components of the observed series y
% using the Hodrick-Prescott filter with smoothing parameter lambda.
%
% All series are Tx1 vectors with the t=1 component at the top and the t=T
% component at the bottom.
%
% Program written for MATLAB 7.10.0 (R1010a) by Peter Ireland, Department
% of Economics, Boston College (irelandp@bc.edu).

bigt = length(y);

biglam = zeros(bigt,bigt);

biglam(1,1:3) = [ 1+lambda -2*lambda lambda ];
biglam(2,1:4) = [ -2*lambda 1+5*lambda -4*lambda lambda ];

for t = 3:bigt-2
    
    biglam(t,t-2:t+2) = [ lambda -4*lambda 1+6*lambda -4*lambda lambda ];
    
end

biglam(bigt-1,bigt-3:bigt) = [ lambda -4*lambda 1+5*lambda -2*lambda ];
biglam(bigt,bigt-2:bigt) = [ lambda -2*lambda 1+lambda ];

g = biglam\y;

c = y - g;
