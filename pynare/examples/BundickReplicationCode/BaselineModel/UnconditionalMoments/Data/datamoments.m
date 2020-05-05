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
% 19		Wu-Xia ShadowRate
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
expostrealrate = data(:,22);
expostequitypremium = data(:,24);

startdate = find(date == 1986);
enddate = find(date == 2014+3/4);


[logoutputtrend,logoutputcycle] = hpfilter(logoutput,1600);
[logconsumptiontrend,logconsumptioncycle] = hpfilter(logconsumption,1600);
[loginvestmenttrend,loginvestmentcycle] = hpfilter(loginvestment,1600);
[loghourstrend,loghourscycle] = hpfilter(loghours,1600);


actualdata       = [100*logoutputcycle 100*logconsumptioncycle 100*loginvestmentcycle 100*loghourscycle expostrealrate expostequitypremium vxo];

datameans  = mean(actualdata(startdate:enddate,:))';

datastandarddeviations  = std(actualdata(startdate:enddate,:))';

clc

disp('Empirical Moments')

disp('  ')

disp([char('Variable','Y:   ','C:    ','I: ','N:   ','Ex-Post RR:   ','Ex-Post Equity Premium:   ','VXO:   ') strvcat('    Mean   ',num2str(datameans))])

disp('  ')

disp([char('Variable','Y:   ','C:    ','I: ','N:   ','Ex-Post RR:   ','Ex-Post Equity Premium:   ','VXO:   ') strvcat('Standard Deviation',num2str(datastandarddeviations))])
     

disp('  ')

datalength        = size(actualdata(startdate:enddate,:),1);

rollingwindowlength     = 20;
rollingstddeviations    = zeros(datalength,size(actualdata,2)-1); %Don't calculate rolling for VXO due to short sample
% 
% 
for idata = 1:size(actualdata,2)-1;
    for irolling = 1:datalength
            rollingstddeviations(irolling,idata) = std(actualdata(startdate-rollingwindowlength+irolling:startdate+irolling-1,idata));
    end
end

stddeviationrollingstddeviations = std(rollingstddeviations)';

disp([char('Variable','Y:   ','C:    ','I: ','N:   ','Ex-Post RR:   ','Ex-Post Equity Premium:   ') strvcat('Std Dev(Rolling Std Dev)',num2str(stddeviationrollingstddeviations))])
     
disp('  ')

