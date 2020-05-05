%% STRUCTURAL VAR
%% SET PRELIMINARY CHOICES
ispl = 0;			% 0 DOESN'T PLOT THE IRFS
save_output = 0;	% 0 DOESN'T SAVE IT
save_graphs   = 0; 	% 0 DOESN'T SAVE THEM
save_marker   = 'Uncertainty'; % THIS CAN BE ANY STRING
%% FIRST LOAD THE DATA
data = csvread('VARData.csv',1,0);
% Data		Description
% 1		Date
% 2		CBOE Market Volatility Index, VXO (Avg) 
% 3		Real Gross Domestic Product (SAAR, Bil.Chn.2009$)  
% 4		Real Personal Consumption Expenditures (SAAR, Bil.Chn.2009.$)  
% 5		Real Personal Consumption Expenditures: Durable Goods (SAAR, Bil.Chn.2009.$)  
% 6		Real Personal Consumption Expenditures: Nondurable Goods (SAAR, Bil.Chn.2009.$)  
% 7		Real Personal Consumption Expenditures: Services (SAAR, Bil.Chn.2009.$)  
% 8		Real Private Fixed Investment (SAAR, Bil.Chn.2009$)  
% 9		Real Private Nonresidential Fixed Investment (SAAR, Bil.Chn.2009$)  
% 10		Gross Domestic Product: Chain Price Index (SA, 2009=100)  
% 11		Federal Funds [Effective] Rate (Avg, % p.a.) 
% 12		5-Year Treasury Note Yield at Constant Maturity (Avg, % p.a.) 
% 13		Aggregate Hours: Nonfarm Payrolls, Total (SAAR, Bil.Hrs)  
% 14		Money Stock: M2 (SA, Bil.$) 
% 15		Civilian Noninstitutional Population: 16 Years and Over (NSA, Thous)  
% 16		CBOE Market Volatility Index, VIX (Avg)
% 17		Stock Price Index: Standard & Poor's 100 (Avg, Close, Jan-2-76=100)
% 18		Standard & Poor's 500 Stock Price Index (Avg, 1941-43=10)
% 19		Wu-Xia Shadow Rate
% 20		Bloom (2009) Spliced Actual/Implied Volatility
% 21		GDP Deflator Inflation Rate (Annualized)
% 22		Ex-Post Real Rate (Annualized)
% 23		Ex-Post S&P 500 Equity Return (Annualized)
% 24		Ex-Post Equity Premium (Annualized)

date = data(:,1);
vxo = data(:,2);
logvxo = log(data(:,2));

% Convert to per-capita and take logs
% Consumption = Nondurables & Services 
% Investment = Durables & Nonresidential Investment

population = data(:,15)/1000000;
logoutput = log(data(:,3)./population);
logconsumption = log((data(:,6)+data(:,7))./population);
loginvestment = log((data(:,5)+data(:,8))./population);
loghours = log(data(:,13)./population);
logprices = log(data(:,10));
logm2 = log(data(:,14));
shadowrate = data(:,19);

startdate = find(date == 1986);
enddate = find(date == 2014+3/4);
%% DEFINE VARIABLES
vardata = [logvxo logoutput logconsumption loginvestment loghours logprices logm2 shadowrate];
vardata = vardata(startdate:enddate,:);
%% LOAD INTO 'Y' (nperiods x nvars) (HERE YOU CAN CHANGE THE ORDERING OF THE VARIABLES, ORDERINGS MATTER FOR IDENTIFYING STRUCTURAL SHOCKS)
Y = vardata;
%% CREATE LABELS
var_labels = [];
shock_labels = [];
%% SET ESTIMATION OPTIONS
%% SET LAG LENGTH
vlag = 4;
%% ESTIMATE WITH CONSTANT BUT NO TIME TREND
iscon = 1;
istr  = 0;
%% NOW ESTIMATE
[Bcomp, cvec, dvec, Bpl_ev, VC_eps, Residmat] = ...
									estim(Y, vlag, iscon, istr);
%% SET HISTORICAL DECOMP OPTIONS
[Tbig,Nbig] = size(Y);
hist_decomp_opt = zeros(2,max([4,Nbig]));
hist_decomp_opt(1,1) = 0; 			    % 0 START FORECAST IN FIRST PERIOD
hist_decomp_opt(1,2) = 0; 				% SUBTRACT BASE FORECAST (0 to leave it)
hist_decomp_opt(1,3) = 0; 				% PLOT THE OUTPUT (0 to not plot it)
hist_decomp_opt(1,4) = 1; 				% COMPUTE HIST DECOMPS (1 not to compute it)
hist_decomp_opt(2,:) = zeros(1,Nbig);	% COMBINE THE SHOCKS WITH A 1 (0 to not include)
[hist_decomp,start_date,struct_shocks] = get_hist_decomp(Bcomp, VC_eps, Y, cvec, hist_decomp_opt);

close all

figure(1)
subplot(2,1,1)
hold on
plot(date(startdate+vlag:enddate),vxo(startdate+vlag:enddate),'-','Color','b','LineWidth',3)
plot(date(startdate+vlag:enddate),mean(vxo(startdate+vlag:enddate))*ones(length(startdate+vlag:enddate),1),'--','Color','k','LineWidth',1)
title('Implied Stock Market Volatility','FontSize',18);xlim([1987 2014+3/4]);set(gca,'XTick',[1987:2:2013],'FontSize',12);
ylim([0 70]);set(gca,'YTick',[0:10:70],'FontSize',14);ylabel('Annualized Percent','FontSize',14);

subplot(2,1,2)
hold on
plot(date(startdate+vlag:enddate),struct_shocks(:,1),'-','Color','r','LineWidth',3)
plot(date(startdate+vlag:enddate),zeros(length(startdate+vlag:enddate),1),'--','Color','k','LineWidth',1)
title('Estimated Uncertainty Shocks','FontSize',18);xlim([1987 2014+3/4]);set(gca,'XTick',[1987:2:2013],'FontSize',14);ylabel('Standard Deviation','FontSize',14)
ylim([-3 5]);set(gca,'YTick',[-3:1:5],'FontSize',14);

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [-0.75 -.65 12.0 9.25]);
print('-dpdf',['VXOLevelShocks.pdf'])










