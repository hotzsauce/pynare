
clear
clc
close all

dynare bbeffectivedemandsimulation.mod noclearall

samplelength = 120;
nsample      = 10000;

estimatedirfy              = zeros(20,nsample);
estimatedirfmu             = zeros(20,nsample);
estimatedirfvxo            = zeros(20,nsample);

for isample = 1:nsample

    samplestart = 250+(isample-1)*250;
    sampleend   = samplestart + samplelength - 1;
    
    vardata = [log(vxo(samplestart:sampleend)) log(y(samplestart:sampleend)) log(c(samplestart:sampleend)) log(inv(samplestart:sampleend))...
                log(n(samplestart:sampleend)) log(cumprod(pie(samplestart:sampleend))) 400*log(r(samplestart:sampleend)) log(mu(samplestart:sampleend))];
     
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
    
    
    
   estimatedirfy(:,isample)     = IRpoint(:,2,1);
   estimatedirfmu(:,isample)    = IRpoint(:,8,1);
   estimatedirfvxo(:,isample)   = IRpoint(:,1,1);

   
end



irflength = IRhoriz+1;
t = 1:irflength;

figure(1)
subplot(2,3,1)
hold on
plot(t,100*median(estimatedirfy,2)','--','Color','r','LineWidth',3)
plot(t,100*prctile(estimatedirfy',2.5)','-','Color','b','LineWidth',3)
plot(t,100*prctile(estimatedirfy',97.5)','-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Output','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.4 0.4]);set(gca,'YTick',[-0.4:0.2:0.4],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,2)
hold on
plot(t,100*median(estimatedirfmu,2),'--','Color','r','LineWidth',3)
plot(t,100*prctile(estimatedirfmu',2.5)','-','Color','b','LineWidth',3)
plot(t,100*prctile(estimatedirfmu',97.5)','-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Markup','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.1 0.2]);set(gca,'YTick',[-0.1:0.1:0.2],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,3)
hold on
plot(t,100*median(estimatedirfvxo,2),'--','Color','r','LineWidth',3)
plot(t,100*prctile(estimatedirfvxo',2.5)','-','Color','b','LineWidth',3)
plot(t,100*prctile(estimatedirfvxo',97.5)','-','Color','b','LineWidth',3)
title('Implied Stock Market Volatility','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-20 60]);set(gca,'YTick',[-20:20:60],'FontSize',12);ylabel('Percent','FontSize',12)
legend('Model - True Response',sprintf('Estimated Response \nUsing SVAR on Simulated Data'),'95% Probability Interval','Location','NorthWest')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [-0.75 -.65 12.0 9.25]);
print('-dpdf',['ModelvsSimulatedVARMarkup.pdf'])







% 
    