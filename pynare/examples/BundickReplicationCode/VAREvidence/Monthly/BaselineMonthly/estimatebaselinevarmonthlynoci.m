%% STRUCTURAL VAR
%% SET PRELIMINARY CHOICES
ispl = 0;			% 0 DOESN'T PLOT THE IRFS
save_output = 0;	% 0 DOESN'T SAVE IT
save_graphs   = 0; 	% 0 DOESN'T SAVE THEM
save_marker   = 'Uncertainty'; % THIS CAN BE ANY STRING
%% FIRST LOAD THE DATA
data = csvread('VARData.csv',1,0);

% Data Description
% 1     Date
% 2     CBOE Market Volatility Index, VXO (Avg) 
% 3 	Spliced Monthly GDP (SAAR, Bil.Chn.2009$)  
% 4     Personal Consumption Expenditures (SAAR, Bil.Chn.2009.$)  
% 5     Personal Consumption Expenditures: Durable Goods (SAAR, Bil.Chn.2009.$)  
% 6     Personal Consumption Expenditures: Nondurable Goods (SAAR, Bil.Chn.2009.$)  
% 7     Personal Consumption Expenditures: Services (SAAR, Bil.Chn.2009.$)  
% 8     Federal Funds [Effective] Rate (Avg, % p.a.) 
% 9     5-Year Treasury Note Yield at Constant Maturity (Avg, % p.a.) 
% 10	Agg Wkly Hrs Production & Nonsupervisory Employees: Total Private(SA, Bils. Hrs)  
% 11	Money Stock: M2 (SA, Bil.$) 
% 12	Civilian Noninstitutional Population: 16 Years and Over (NSA, Thous)  
% 13	PCE: Chain Price Index (SA, 2009=100)  
% 14	PCE less Food & Energy: Chain Price Index (SA, 2009=100)  
% 15	Wu-Xia Shadow Rate

date = data(:,1);
vxo = data(:,2);
logvxo = log(data(:,2));

% Convert to per-capita and take logs
% Consumption = Nondurables & Services 
% Investment = Durables & Nonresidential Investment

population = data(:,12)/1000000;
logoutput = log(data(:,3)./population);
logconsumption = log((data(:,6)+data(:,7))./population);
logdurableconsumption = log(data(:,5)./population);
loghours = log(data(:,10)./population);
logprices = log(data(:,14));
logm2 = log(data(:,11));
shadowrate = data(:,15);

startdate = find(date == 1986);
enddate = length(date);
%% DEFINE VARIABLES
vardata = [logvxo logoutput logconsumption logdurableconsumption loghours logprices logm2 shadowrate];
vardata = vardata(startdate:enddate,:);
%% LOAD INTO 'Y' (nperiods x nvars) (HERE YOU CAN CHANGE THE ORDERING OF THE VARIABLES, ORDERINGS MATTER FOR IDENTIFYING STRUCTURAL SHOCKS)
Y = vardata;
[Tbig,Nbig] = size(Y);

%% CREATE LABELS
var_labels = [];
shock_labels = [];
%% SET ESTIMATION OPTIONS
%% SET LAG LENGTH
vlag = 12;
%% ESTIMATE WITH CONSTANT BUT NO TIME TREND
iscon = 1;
istr  = 0;
%% NOW ESTIMATE
[Bcomp, cvec, dvec, Bpl_ev, VC_eps, Residmat] = estim(Y, vlag, iscon, istr);

%% SET NUMBER OF PERIODS (adds first period so 15 implies 16 periods will be computed)
IRhoriz = 59;
%% SET TYPE OF SHOCK (1 PERCENT or 1 STD DEV)
IRtype = 'c';
sizesho = 0;

IRpoint = VAR_irf(Bcomp, VC_eps, IRhoriz, IRtype, sizesho, sizesho);

monthlyobs = IRhoriz+1;
quarterlyobs = monthlyobs/3;

IRpointmonthly = IRpoint(:,:,1);

IRpointquarterly = zeros(quarterlyobs,Nbig);

for i = 1:quarterlyobs
    IRpointquarterly(i,:) = (IRpointmonthly(3*(i-1)+1,:) + IRpointmonthly(3*(i-1)+2,:) + IRpointmonthly(3*(i-1)+3,:))./3;
end

irflength = quarterlyobs;
t = 1:irflength;

figure(1)
subplot(2,3,1)
hold on
plot(t,100*IRpointquarterly(:,2),'s','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
title('Output','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);

subplot(2,3,2)
hold on
plot(t,100*IRpointquarterly(:,3),'s','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
title('Consumption','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);

subplot(2,3,4)
hold on
plot(t,100*IRpointquarterly(:,5),'s','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
title('Hours','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);


subplot(2,3,5)
hold on
plot(t,100*IRpointquarterly(:,1),'s','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
title('Implied Stock Market Volatility','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);










