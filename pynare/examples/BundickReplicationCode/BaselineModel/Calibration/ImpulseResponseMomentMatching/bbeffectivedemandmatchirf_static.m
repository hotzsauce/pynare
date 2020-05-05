function [residual, g1, g2] = bbeffectivedemandmatchirf_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 27, 1);

%
% Model equations
%

T18 = params(23)*(y(27)*y(12))^(1-params(1));
T22 = (y(21)*y(10))^params(1);
T50 = y(2)^params(6)*(1-y(12))^(1-params(6));
T55 = T50^((1-params(9))/params(12));
T62 = params(24)*y(1)*T55+params(2)*y(8)^(1/params(12));
T66 = y(23)^(1-params(9));
T89 = y(9)/y(10);
T94 = 1-y(4)-params(8)/2*(T89-params(3))^2;
T122 = y(14)^(-1);
T137 = y(14)/params(38);
T152 = params(7)/2*(T137-1)^2;
T160 = params(8)*(T89-params(3));
T192 = y(13)^2;
T233 = getPowerDeriv(T62,params(12)/(1-params(9)),1);
T240 = getPowerDeriv(T50,(1-params(9))/params(12),1);
T268 = getPowerDeriv(T66/y(8),1-1/params(12),1);
T276 = (-(params(8)/2*1/y(10)*2*(T89-params(3))));
T301 = (-(params(8)/2*2*(T89-params(3))*(-y(9))/(y(10)*y(10))));
T365 = 1/params(38);
lhs =y(26)+params(15);
rhs =T18*T22;
residual(1)= lhs-rhs;
lhs =y(2)+y(10)*params(14)/y(18);
rhs =y(10)*params(14)+y(12)*y(25)+y(3);
residual(2)= lhs-rhs;
lhs =y(25);
rhs =y(2)*(1-params(6))/params(6)/(1-y(12));
residual(3)= lhs-rhs;
lhs =y(23);
rhs =T62^(params(12)/(1-params(9)));
residual(4)= lhs-rhs;
lhs =y(8);
rhs =T66;
residual(5)= lhs-rhs;
lhs =y(12)*y(25);
rhs =(y(26)+params(15))*(1-params(1))/y(11);
residual(6)= lhs-rhs;
lhs =y(10)*y(21)*y(19);
rhs =(y(26)+params(15))*params(1)/y(11);
residual(7)= lhs-rhs;
lhs =y(10)*y(21)*y(16)*y(5);
rhs =(y(26)+params(15))*params(1)/y(11);
residual(8)= lhs-rhs;
lhs =y(10);
rhs =y(9)+y(10)*T94;
residual(9)= lhs-rhs;
lhs =y(4);
rhs =params(3)+params(4)*(y(21)-1)+params(5)/2*(y(21)-1)^2;
residual(10)= lhs-rhs;
lhs =y(5);
rhs =params(4)+(y(21)-1)*params(5);
residual(11)= lhs-rhs;
lhs =y(20);
rhs =params(2)*(T66/y(8))^(1-1/params(12));
residual(12)= lhs-rhs;
lhs =1;
rhs =y(18)*y(20);
residual(13)= lhs-rhs;
lhs =1;
rhs =y(20)*y(17)*T122;
residual(14)= lhs-rhs;
lhs =1;
rhs =y(20)*(y(3)+y(13))/y(13);
residual(15)= lhs-rhs;
lhs =log(y(17));
rhs =(1-params(16))*(log(params(41))+params(17)*log(T137))+log(y(17))*params(16);
residual(16)= lhs-rhs;
lhs =y(3);
rhs =y(26)-y(12)*y(25)-y(9)-y(26)*T152-params(14)*(y(10)-y(10)/y(18));
residual(17)= lhs-rhs;
lhs =1;
rhs =y(20)*(y(21)*y(19)+y(16)*(T94+T89*T160))/y(16);
residual(18)= lhs-rhs;
lhs =1/y(16);
rhs =1-T160;
residual(19)= lhs-rhs;
lhs =T137*params(7)*(T137-1);
rhs =1-params(13)+params(13)/y(11)+T137*(T137-1)*y(20)*params(7);
residual(20)= lhs-rhs;
lhs =y(15);
rhs =y(26)*(y(11)-1)-params(15);
residual(21)= lhs-rhs;
lhs =y(6);
rhs =(y(3)+y(13))/y(13);
residual(22)= lhs-rhs;
lhs =y(7);
rhs =(y(3)+y(13))^2/T192;
residual(23)= lhs-rhs;
lhs =y(22);
rhs =y(7)-y(6)^2;
residual(24)= lhs-rhs;
lhs =y(1);
rhs =(1-params(19))*params(25)+y(1)*params(19)+y(24)*x(1);
residual(25)= lhs-rhs;
lhs =y(24);
rhs =y(24)*params(20)+(1-params(20))*params(48)+params(21)*x(2);
residual(26)= lhs-rhs;
lhs =y(27);
rhs =(1-params(22))*params(52)+y(27)*params(22)+params(49)*x(3);
residual(27)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(27, 27);

  %
  % Jacobian matrix
  %

  g1(1,10)=(-(T18*y(21)*getPowerDeriv(y(21)*y(10),params(1),1)));
  g1(1,12)=(-(T22*params(23)*y(27)*getPowerDeriv(y(27)*y(12),1-params(1),1)));
  g1(1,21)=(-(T18*y(10)*getPowerDeriv(y(21)*y(10),params(1),1)));
  g1(1,26)=1;
  g1(1,27)=(-(T22*params(23)*y(12)*getPowerDeriv(y(27)*y(12),1-params(1),1)));
  g1(2,2)=1;
  g1(2,3)=(-1);
  g1(2,10)=params(14)/y(18)-params(14);
  g1(2,12)=(-y(25));
  g1(2,18)=(-(y(10)*params(14)))/(y(18)*y(18));
  g1(2,25)=(-y(12));
  g1(3,2)=(-((1-params(6))/params(6)/(1-y(12))));
  g1(3,12)=(-(y(2)*(1-params(6))/params(6)/((1-y(12))*(1-y(12)))));
  g1(3,25)=1;
  g1(4,1)=(-(params(24)*T55*T233));
  g1(4,2)=(-(T233*params(24)*y(1)*(1-y(12))^(1-params(6))*getPowerDeriv(y(2),params(6),1)*T240));
  g1(4,8)=(-(T233*params(2)*getPowerDeriv(y(8),1/params(12),1)));
  g1(4,12)=(-(T233*params(24)*y(1)*T240*y(2)^params(6)*(-(getPowerDeriv(1-y(12),1-params(6),1)))));
  g1(4,23)=1;
  g1(5,8)=1;
  g1(5,23)=(-(getPowerDeriv(y(23),1-params(9),1)));
  g1(6,11)=(-((-((y(26)+params(15))*(1-params(1))))/(y(11)*y(11))));
  g1(6,12)=y(25);
  g1(6,25)=y(12);
  g1(6,26)=(-((1-params(1))/y(11)));
  g1(7,10)=y(21)*y(19);
  g1(7,11)=(-((-((y(26)+params(15))*params(1)))/(y(11)*y(11))));
  g1(7,19)=y(21)*y(10);
  g1(7,21)=y(10)*y(19);
  g1(7,26)=(-(params(1)/y(11)));
  g1(8,5)=y(10)*y(21)*y(16);
  g1(8,10)=y(21)*y(16)*y(5);
  g1(8,11)=(-((-((y(26)+params(15))*params(1)))/(y(11)*y(11))));
  g1(8,16)=y(10)*y(21)*y(5);
  g1(8,21)=y(10)*y(16)*y(5);
  g1(8,26)=(-(params(1)/y(11)));
  g1(9,4)=y(10);
  g1(9,9)=(-(1+y(10)*T276));
  g1(9,10)=1-(T94+y(10)*T301);
  g1(10,4)=1;
  g1(10,21)=(-(params(4)+params(5)/2*2*(y(21)-1)));
  g1(11,5)=1;
  g1(11,21)=(-params(5));
  g1(12,8)=(-(params(2)*(-T66)/(y(8)*y(8))*T268));
  g1(12,20)=1;
  g1(12,23)=(-(params(2)*T268*getPowerDeriv(y(23),1-params(9),1)/y(8)));
  g1(13,18)=(-y(20));
  g1(13,20)=(-y(18));
  g1(14,14)=(-(y(20)*y(17)*getPowerDeriv(y(14),(-1),1)));
  g1(14,17)=(-(y(20)*T122));
  g1(14,20)=(-(y(17)*T122));
  g1(15,3)=(-(y(20)/y(13)));
  g1(15,13)=(-((y(20)*y(13)-y(20)*(y(3)+y(13)))/(y(13)*y(13))));
  g1(15,20)=(-((y(3)+y(13))/y(13)));
  g1(16,14)=(-((1-params(16))*params(17)*T365/T137));
  g1(16,17)=1/y(17)-params(16)*1/y(17);
  g1(17,3)=1;
  g1(17,9)=1;
  g1(17,10)=params(14)*(1-1/y(18));
  g1(17,12)=y(25);
  g1(17,14)=y(26)*params(7)/2*T365*2*(T137-1);
  g1(17,18)=params(14)*(-((-y(10))/(y(18)*y(18))));
  g1(17,25)=y(12);
  g1(17,26)=(-(1-T152));
  g1(18,4)=(-(y(20)*(-y(16))/y(16)));
  g1(18,9)=(-(y(20)*y(16)*(T276+T160*1/y(10)+T89*params(8)*1/y(10))/y(16)));
  g1(18,10)=(-(y(20)*y(16)*(T301+T160*(-y(9))/(y(10)*y(10))+T89*params(8)*(-y(9))/(y(10)*y(10)))/y(16)));
  g1(18,16)=(-((y(16)*y(20)*(T94+T89*T160)-y(20)*(y(21)*y(19)+y(16)*(T94+T89*T160)))/(y(16)*y(16))));
  g1(18,19)=(-(y(21)*y(20)/y(16)));
  g1(18,20)=(-((y(21)*y(19)+y(16)*(T94+T89*T160))/y(16)));
  g1(18,21)=(-(y(19)*y(20)/y(16)));
  g1(19,9)=params(8)*1/y(10);
  g1(19,10)=params(8)*(-y(9))/(y(10)*y(10));
  g1(19,16)=(-1)/(y(16)*y(16));
  g1(20,11)=(-((-params(13))/(y(11)*y(11))));
  g1(20,14)=params(7)*(T137-1)*T365+T137*params(7)*T365-((T137-1)*y(20)*params(7)*T365+T137*y(20)*params(7)*T365);
  g1(20,20)=(-(T137*params(7)*(T137-1)));
  g1(21,11)=(-y(26));
  g1(21,15)=1;
  g1(21,26)=(-(y(11)-1));
  g1(22,3)=(-(1/y(13)));
  g1(22,6)=1;
  g1(22,13)=(-((y(13)-(y(3)+y(13)))/(y(13)*y(13))));
  g1(23,3)=(-(2*(y(3)+y(13))/T192));
  g1(23,7)=1;
  g1(23,13)=(-((T192*2*(y(3)+y(13))-(y(3)+y(13))^2*2*y(13))/(T192*T192)));
  g1(24,6)=2*y(6);
  g1(24,7)=(-1);
  g1(24,22)=1;
  g1(25,1)=1-params(19);
  g1(25,24)=(-x(1));
  g1(26,24)=1-params(20);
  g1(27,27)=1-params(22);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],27,729);
end
end
