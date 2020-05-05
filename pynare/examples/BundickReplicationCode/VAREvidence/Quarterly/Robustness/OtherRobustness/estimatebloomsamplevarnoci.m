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
vxo = data(:,20);
logvxo = log(data(:,20));

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
fedfundsrate = data(:,11);

startdate = find(date == 1962+1/2);
enddate = find(date == 2008+1/4);
%% DEFINE VARIABLES
vardata = [logvxo logoutput logconsumption loginvestment loghours logprices fedfundsrate];
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
[Bcomp, cvec, dvec, Bpl_ev, VC_eps, Residmat] = estim(Y, vlag, iscon, istr);

%% SET NUMBER OF PERIODS (adds first period so 15 implies 16 periods will be computed)
IRhoriz = 19;
%% SET TYPE OF SHOCK (1 PERCENT or 1 STD DEV)
IRtype = 'c';
sizesho = 0;

IRpoint = VAR_irf(Bcomp, VC_eps, IRhoriz, IRtype, sizesho, sizesho);

irflength = IRhoriz+1;
t = 1:irflength;

figure(1)
subplot(2,3,1)
hold on
plot(t,100*IRpoint(:,2,1),'d','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1,'HandleVisibility','off')
title('Output','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.4 0.4]);set(gca,'YTick',[-0.4:0.2:0.4],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,2)
hold on
plot(t,100*IRpoint(:,3,1),'d','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1,'HandleVisibility','off')
title('Consumption','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.4 0.4]);set(gca,'YTick',[-0.6:0.2:0.6],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,3)
hold on
plot(t,100*IRpoint(:,4,1),'d','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1,'HandleVisibility','off')
title('Investment','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-1.5 1]);set(gca,'YTick',[-1.5:0.5:1],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,4)
hold on
plot(t,100*IRpoint(:,5,1),'d','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1,'HandleVisibility','off')
title('Hours','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.4 0.4]);set(gca,'YTick',[-0.4:0.2:0.4],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,5)
hold on
plot(t,100*IRpoint(:,1,1),'d','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0],'MarkerSize',7)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1,'HandleVisibility','off')
title('Implied Stock Market Volatility','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-5 20]);set(gca,'YTick',[-5:5:20],'FontSize',12);ylabel('Percent','FontSize',12)

