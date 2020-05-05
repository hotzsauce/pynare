
clear
clc
close all

dynare bbeffectivedemandsimulation.mod noclearall

samplelength = 120;
nsample     = 10000;

shockcorrelation = zeros(nsample,1);

estimatedirfy    = zeros(20,nsample);
estimatedirfc    = zeros(20,nsample);
estimatedirfinv  = zeros(20,nsample);
estimatedirfn    = zeros(20,nsample);
estimatedirfvxo  = zeros(20,nsample);

comovementsamples = 0;

for isample = 1:nsample

    samplestart = 250+(isample-1)*250;
    sampleend   = samplestart + samplelength - 1;
    
    vardata = [log(vxo(samplestart:sampleend)) log(y(samplestart:sampleend)) log(c(samplestart:sampleend)) log(inv(samplestart:sampleend))...
                log(n(samplestart:sampleend)) log(cumprod(pie(samplestart:sampleend))) 400*log(r(samplestart:sampleend))];
     
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
    
    

    %% SET HISTORICAL DECOMP OPTIONS
    [Tbig,Nbig] = size(Y);
    hist_decomp_opt = zeros(2,max([4,Nbig]));
    hist_decomp_opt(1,1) = 0; 			    % 0 START FORECAST IN FIRST PERIOD
    hist_decomp_opt(1,2) = 0; 				% SUBTRACT BASE FORECAST (0 to leave it)
    hist_decomp_opt(1,3) = 0; 				% PLOT THE OUTPUT (0 to not plot it)
    hist_decomp_opt(1,4) = 1; 				% COMPUTE HIST DECOMPS (1 not to compute it)
    hist_decomp_opt(2,:) = zeros(1,Nbig);	% COMBINE THE SHOCKS WITH A 1 (0 to not include)
    [hist_decomp,start_date,struct_shocks] = get_hist_decomp(Bcomp, VC_eps, Y, cvec, hist_decomp_opt);
               
    evola   = oo_.exo_simul(samplestart + vlag:sampleend,2);        

    shockcorrelation(isample) = corr(struct_shocks(:,1),evola);
 
   estimatedirfy(:,isample)     = IRpoint(:,2,1);
   estimatedirfc(:,isample)     = IRpoint(:,3,1);
   estimatedirfinv(:,isample)   = IRpoint(:,4,1);
   estimatedirfn(:,isample)     = IRpoint(:,5,1);
   estimatedirfvxo(:,isample)   = IRpoint(:,1,1);
   
      if (max(IRpoint(1,2:5,1)) < 0) 
       comovementsamples = comovementsamples + 1;
      end
       

end

disp(' ')
disp(['Shock Correlation: ',num2str(mean(shockcorrelation))])
disp(['95% CI:   (',num2str(prctile(shockcorrelation,2.5)),',',num2str(prctile(shockcorrelation,97.5)),')'])
disp(' ')

disp(' ')
disp(['% of SVARs With Comovement at Impact: ',num2str(comovementsamples/nsample)])
disp(' ')

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
plot(t,100*median(estimatedirfc,2),'--','Color','r','LineWidth',3)
plot(t,100*prctile(estimatedirfc',2.5)','-','Color','b','LineWidth',3)
plot(t,100*prctile(estimatedirfc',97.5)','-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Consumption','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.2 0.2]);set(gca,'YTick',[-0.3:0.1:0.3],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,3)
hold on
plot(t,100*median(estimatedirfinv,2),'--','Color','r','LineWidth',3)
plot(t,100*prctile(estimatedirfinv',2.5)','-','Color','b','LineWidth',3)
plot(t,100*prctile(estimatedirfinv',97.5)','-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Investment','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-1.5 1.5]);set(gca,'YTick',[-1.5:0.5:1.5],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,4)
hold on
plot(t,100*median(estimatedirfn,2),'--','Color','r','LineWidth',3)
plot(t,100*prctile(estimatedirfn',2.5)','-','Color','b','LineWidth',3)
plot(t,100*prctile(estimatedirfn',97.5)','-','Color','b','LineWidth',3)
plot(t,zeros(irflength,1),'--','Color','k','LineWidth',1)
title('Hours','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.3 0.3]);set(gca,'YTick',[-0.3:0.1:0.3],'FontSize',12);ylabel('Percent','FontSize',12)

subplot(2,3,5)
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
print('-dpdf',['ModelvsSimulatedVAR.pdf'])







% 
    