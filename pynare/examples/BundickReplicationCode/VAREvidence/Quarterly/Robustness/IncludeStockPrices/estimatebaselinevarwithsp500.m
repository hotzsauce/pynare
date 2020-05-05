clear

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
logsp500 = log(data(:,18));

startdate = find(date == 1986);
enddate = find(date == 2014+3/4);
%% DEFINE VARIABLES
vardata = [logvxo logoutput logconsumption loginvestment loghours logprices logsp500 shadowrate];
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

%% SET NUMBER OF PERIODS (adds first period so 15 implies 16 periods will be computed)
IRhoriz = 19;
%% SET TYPE OF SHOCK (1 PERCENT or 1 STD DEV)
IRtype = 'c';
sizesho = 0;
%% SET OPTIONS FOR CONFIDENCE BOUNDS (MONTE CARLO INTEGRATION OF CONFIDENCE BANDS, SEE KOOP'S CHAPTER IN THE HANDBOOK OF BAYESIAN ECONOMETRICS)
NMC    = 5000;							% NUMBER OF DRAWS FROM THE POSTERIOR DISTRIBUTION
CIperc = [.025 , .16, .84,  .975];		% CHOOSE WHICH PERCENTILES TO PLOT (CAN CHOOSE 2 or 4)
num_auto_corr = 5;						% HOW MANY AUTO-CORRELATIONS TO COMPUTE
no_moments = 0; 						% COMPUTE MOMENTS
var_decomp_opt = 1:1:40;  				% COMPUTE FEVD AT THESE HORIZONS
var_decomp_opt = var_decomp_opt';		
%% SET PLOT OPTIONS
% USE plot_options TO CHANGE PLOTTING LAYOUT	(COLUMNS REPRESENTS SHOCKS, ROWS REPRESENT VARIABLES - NEGATIVES PLOT NEGATIVE SHOCKS)
plot_options = zeros(Nbig,Nbig);
%%
[IRpoint, IRup_2, IRup_1, IRlo_1, IRlo_2, IRmed, IRstd, Bevmean, Bevstd, auto_corr, moments, var_decomp] = ...
confid_int(Bcomp,VC_eps,cvec,dvec,IRhoriz,IRtype,sizesho,NMC,Tbig,CIperc,ispl,num_auto_corr,no_moments,var_decomp_opt,hist_decomp_opt,Y,plot_options,save_graphs,save_marker,var_labels,shock_labels);

irflength = IRhoriz+1;
t = 1:irflength;

close all

save baselinevarestimationwithsp500.mat

figure(1)
subplot(3,3,1)
hold on
plot(t,100*IRpoint(:,2,1),'--','Color','r','LineWidth',3)
plot(t,100*IRlo_1(:,2,1),'-','Color','b','LineWidth',3)
plot(t,100*IRup_1(:,2,1),'-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Output','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.6 0.6]);set(gca,'YTick',[-0.6:0.2:0.6],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,2)
hold on
plot(t,100*IRpoint(:,3,1),'--','Color','r','LineWidth',3)
plot(t,100*IRlo_1(:,3,1),'-','Color','b','LineWidth',3)
plot(t,100*IRup_1(:,3,1),'-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Consumption','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.4 0.4]);set(gca,'YTick',[-0.6:0.2:0.6],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,3)
hold on
plot(t,100*IRpoint(:,4,1),'--','Color','r','LineWidth',3)
plot(t,100*IRlo_1(:,4,1),'-','Color','b','LineWidth',3)
plot(t,100*IRup_1(:,4,1),'-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Investment','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-1.5 2]);set(gca,'YTick',[-1.5:0.5:2],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,4)
hold on
plot(t,100*IRpoint(:,5,1),'--','Color','r','LineWidth',3)
plot(t,100*IRlo_1(:,5,1),'-','Color','b','LineWidth',3)
plot(t,100*IRup_1(:,5,1),'-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Hours','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.6 0.6]);set(gca,'YTick',[-0.6:0.2:0.6],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,5)
hold on
plot(t,100*IRpoint(:,6,1),'--','Color','r','LineWidth',3)
plot(t,100*IRlo_1(:,6,1),'-','Color','b','LineWidth',3)
plot(t,100*IRup_1(:,6,1),'-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Price Level','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.3 0.3]);set(gca,'YTick',[-0.3:0.1:0.3],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,6)
hold on
plot(t,IRpoint(:,8,1),'--','Color','r','LineWidth',3)
plot(t,IRlo_1(:,8,1),'-','Color','b','LineWidth',3)
plot(t,IRup_1(:,8,1),'-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Policy Rate','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.4 0.4]);set(gca,'YTick',[-0.4:0.2:0.4],'FontSize',12);ylabel('Percentage Points','FontSize',12)


subplot(3,3,7)
hold on
plot(t,100*IRpoint(:,7,1),'--','Color','r','LineWidth',3)
plot(t,100*IRlo_1(:,7,1),'-','Color','b','LineWidth',3)
plot(t,100*IRup_1(:,7,1),'-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Stock Prices','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
%ylim([-10 20]);set(gca,'YTick',[-10:5:20],'FontSize',12);ylabel('Percent','FontSize',12)


subplot(3,3,8)
hold on
plot(t,100*IRpoint(:,1,1),'--','Color','r','LineWidth',3)
plot(t,100*IRlo_1(:,1,1),'-','Color','b','LineWidth',3)
plot(t,100*IRup_1(:,1,1),'-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Implied Stock Market Volatility','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-10 20]);set(gca,'YTick',[-10:5:20],'FontSize',12);ylabel('Percent','FontSize',12)
legend('Impulse Response','95% Confidence Interval','Location','NorthEast')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'portrait');
set(gcf, 'PaperPosition', [-0.5 -0.75 9.5 12.0]);
print('-dpdf',['BaselineVARWithStockPrices.pdf'])










