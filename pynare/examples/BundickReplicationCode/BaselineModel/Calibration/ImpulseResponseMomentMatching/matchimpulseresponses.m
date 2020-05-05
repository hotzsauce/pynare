
clear
clc

% Setup Dynare

cd OriginalModFiles

dynare bbeffectivedemandmatchirf.mod noclearall
clear
clc

cd ..

% Load empirical impulse responses
% Copied from VAREvidence/Quarterly/Baseline Folder
load irfdata.mat

phikguess               = 2.2009;
rhovolaguess            = 0.7533;
volvolaguess			= 0.0025;
volassguess             = 0.0033;
volzssguess             = 0.0013;
rhoaguess               = 0.9131;
rhozguess               = 0.9903;

% Scale the parameters to be between (0,1)
% Makes numerical optimization easier

parameterguess(1)      = log(phikguess/10);
parameterguess(2)      = log(rhovolaguess);
parameterguess(3)      = log(volvolaguess*100);
parameterguess(4)      = log(volassguess*100);
parameterguess(5)      = log(volzssguess*100);
parameterguess(6)      = log(rhoaguess);
parameterguess(7)      = log(rhozguess);

options=optimset('MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-6,'TolFun',1e-6,'Display','off');

% Impulse response matching

[parametersestimation,distance] = fminsearch(@irfdistance,parameterguess,options,irfdata,weights);


% Transform estimated parameters

phikestimate            = 10*exp(parametersestimation(1));
rhovolaestimate 		= exp(parametersestimation(2));
volvolaestimate			= exp(parametersestimation(3))/100;
volassestimate			= exp(parametersestimation(4))/100;
volzssestimate			= exp(parametersestimation(5))/100;
rhoaestimate            = exp(parametersestimation(6));
rhozestimate            = exp(parametersestimation(7));

finalparameters         = [phikestimate;rhovolaestimate;volvolaestimate;volassestimate;volzssestimate;rhoaestimate;rhozestimate];

disp([' '])

disp([char('Distance: ') num2str(distance)])

disp([' '])

disp([char('Phik:    ',...
           'Rhovola: ',...
           'Volvola: ',...
           'Volass:  ',...
           'Volzss:  ',...
           'Rhoa:    ',...
           'Rhoz:    ') num2str(finalparameters)])

disp([' '])

% Generate IRFs and simulate using final parameter estimates

bbeffectivedemandmatchirf;

% Plot estimated impulse responses versus data

t               = 1:irflength; % irflength is set by computeirfs.m, called by bbeffectivedemandmatchirf.m

vxo_vxodata     = 100*IRpoint(1:irflength,1,1);
y_vxodata       = 100*IRpoint(1:irflength,2,1);
c_vxodata       = 100*IRpoint(1:irflength,3,1);
inv_vxodata     = 100*IRpoint(1:irflength,4,1);
n_vxodata       = 100*IRpoint(1:irflength,5,1);
pl_vxodata      = 100*IRpoint(1:irflength,6,1);
r_vxodata       = IRpoint(1:irflength,8,1);

vxocilo_vxodata     = 100*IRlo_1(1:irflength,1,1);
ycilo_vxodata       = 100*IRlo_1(1:irflength,2,1);
ccilo_vxodata       = 100*IRlo_1(1:irflength,3,1);
invcilo_vxodata     = 100*IRlo_1(1:irflength,4,1);
ncilo_vxodata       = 100*IRlo_1(1:irflength,5,1);
plcilo_vxodata      = 100*IRlo_1(1:irflength,6,1);
rcilo_vxodata       = IRlo_1(1:irflength,8,1);

vxociup_vxodata     = 100*IRup_1(1:irflength,1,1);
yciup_vxodata       = 100*IRup_1(1:irflength,2,1);
cciup_vxodata       = 100*IRup_1(1:irflength,3,1);
invciup_vxodata     = 100*IRup_1(1:irflength,4,1);
nciup_vxodata       = 100*IRup_1(1:irflength,5,1);
plciup_vxodata      = 100*IRup_1(1:irflength,6,1);
rciup_vxodata       = IRup_1(1:irflength,8,1);



figure(1)
subplot(3,3,1)
hold on
plot(t,y_vxodata,'--','Color','r','LineWidth',3)
plot(t,ycilo_vxodata,'-','Color','b','LineWidth',3)
plot(t,yciup_vxodata,'-','Color','b','LineWidth',3,'HandleVisibility','off')
plot(t,100*y_evola,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
ylim([-0.6 0.6]);set(gca,'YTick',[-0.6:0.2:0.6],'FontSize',12);ylabel('Percent','FontSize',12)
title('Output','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12)

subplot(3,3,2)
hold on
plot(t,c_vxodata,'--','Color','r','LineWidth',3)
plot(t,ccilo_vxodata,'-','Color','b','LineWidth',3)
plot(t,cciup_vxodata,'-','Color','b','LineWidth',3,'HandleVisibility','off')
plot(t,100*c_evola,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
ylim([-0.4 0.4]);set(gca,'YTick',[-0.4:0.2:0.4],'FontSize',12);ylabel('Percent','FontSize',12)
title('Consumption','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12)

subplot(3,3,3)
hold on
plot(t,inv_vxodata,'--','Color','r','LineWidth',3)
plot(t,invcilo_vxodata,'-','Color','b','LineWidth',3)
plot(t,invciup_vxodata,'-','Color','b','LineWidth',3,'HandleVisibility','off')
plot(t,100*inv_evola,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
ylim([-2 2]);set(gca,'YTick',[-2:0.5:2],'FontSize',12);ylabel('Percent','FontSize',12)
title('Investment','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12)

subplot(3,3,4)
hold on
plot(t,n_vxodata,'--','Color','r','LineWidth',3)
plot(t,ncilo_vxodata,'-','Color','b','LineWidth',3)
plot(t,nciup_vxodata,'-','Color','b','LineWidth',3,'HandleVisibility','off')
plot(t,100*n_evola,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
ylim([-0.6 0.6]);set(gca,'YTick',[-0.6:0.2:0.6],'FontSize',12);ylabel('Percent','FontSize',12)
title('Hours Worked','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12)


subplot(3,3,5)
hold on
plot(t,pl_vxodata,'--','Color','r','LineWidth',3)
plot(t,plcilo_vxodata,'-','Color','b','LineWidth',3)
plot(t,plciup_vxodata,'-','Color','b','LineWidth',3,'HandleVisibility','off')
plot(t,100*pl_evola,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
ylim([-0.3 0.3]);set(gca,'YTick',[-0.3:0.1:0.3],'FontSize',12);ylabel('Percent','FontSize',12)
title('Price Level','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);


subplot(3,3,6)
hold on
plot(t,r_vxodata,'--','Color','r','LineWidth',3)
plot(t,rcilo_vxodata,'-','Color','b','LineWidth',3)
plot(t,rciup_vxodata,'-','Color','b','LineWidth',3,'HandleVisibility','off')
plot(t,r_evola,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
ylim([-0.4 0.4]);set(gca,'YTick',[-0.4:0.2:0.4],'FontSize',12);ylabel('Percentage Points','FontSize',12)
title('Policy Rate','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);

subplot(3,3,8)
hold on
plot(t,100*vola_evola,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
title('Preference Shock Volatility','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,7)
hold on
plot(t,vxo_vxodata,'--','Color','r','LineWidth',3)
plot(t,vxocilo_vxodata,'-','Color','b','LineWidth',3)
plot(t,vxociup_vxodata,'-','Color','b','LineWidth',3,'HandleVisibility','off')
plot(t,100*logvxo_evola,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
title('Implied Stock Market Volatility','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Percent','FontSize',12)
legend('Data - Estimated Response','Data - 95% Confidence Interval','Model','Location','NorthWest')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'portrait');
set(gcf, 'PaperPosition', [-0.5 -0.75 9.5 12.0]);
print('-dpdf',['BaselineVARvsBaselineModel.pdf'])

save calibrationresults.mat
delete *.jnl





























