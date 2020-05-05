// Basu & Bundick (2016)
// "Uncertainty Shocks in a Model of Effective Demand"

var a, c, de, deltau, deltauprime, expre, expre2, expvf1sigma, inv, k, mu, n,
	pe, pie, profit, q, r, rr, rrk, sdf, u, varexpre, vf, vola, w, y, z;

varexo ea, evola, ez;

parameters 	alpha, beta, delta0, delta1, delta2, eta, phip, phik, 
			sigma, ies, frischlse, thetavf, thetamu, leverageratio, fixedcost,
			rhor, rhopie, rhoy,	rhoa, rhovola, volvola, rhoz, productionconstant, utilityconstant, 
            ass, css, dess, deltauss, deltauprimess, express, expre2ss, expvf1sigmass, invss, kss, nss, muss, 
            pess, piess, profitss, qss, rss, rrss, rrkss, sdfss, uss, varexpress, vfss, volass, volzss, wss, yss, zss;

alpha           = 0.333;
beta            = 0.994;
delta0          = 0.025;
delta1          = (1/beta - 1 + delta0);
delta2          = 0.00031;
phip            = 100;
phik            = 2.0901;
sigma           = 80;
ies             = 0.95;
frischlse       = 2.0;
thetavf         = (1 - sigma)/(1 - 1/ies);
thetamu         = 6.0;
leverageratio   = 0.90;

piess           = 1.005;
rss             = piess/beta;
rhor            = 0.000;
rhopie          = 1.5;
rhoy            = 0.2;

ass             = 1;
rhoa            = 0.93564;
volass          = 0.0026251;
rhovola         = 0.74227;
volvola         = 0.0025022;

zss             = 1;
rhoz            = 0.98793;
volzss          = 0.0012857;

// Steady State Calibration

yss         = 1.0;
rrkss       = 1/beta - 1 + delta0;
muss        = thetamu/(thetamu - 1);
uss			= 1;
fixedcost   = (muss - 1)*yss;
kss         = alpha*(yss + fixedcost)/(rrkss*muss);
invss       = delta0*kss;
css         = yss - invss;

nsssolution = fsolve(@(nssguess) 1+(-1+1/ies)/(1+((-1+alpha)*(-1+nssguess)*(fixedcost+yss))/(css*nssguess*muss)) -(frischlse/ies)*nssguess*(1-nssguess)^(-1),0.5,optimoptions('fsolve','Display','off'));
nss         = nsssolution;

wss                 = (1 - alpha)*(yss + fixedcost)/(nss*muss);
eta                 = (wss*(1 - nss)/css + 1)^(-1);
qss                 = 1;
productionconstant 	= (yss + fixedcost)/(kss^(alpha)*(nss)^(1 - alpha));
vfss                = 1;
utilityconstant     = vfss^((1 - sigma)/thetavf)*(1 - beta)/(css^(eta)*(1 - nss)^(1 - eta))^((1 - sigma)/thetavf);
dess                = css - wss*nss + leverageratio*kss*(beta - 1);
pess                = beta*(1-beta)^(-1)*dess;
profitss            = (muss - 1)*yss - fixedcost;
express             = 1/beta;
expre2ss            = 1/beta^(2);
varexpress          = 0;
expvf1sigmass       = vfss^(1-sigma);
qss                 = 1;
rrss                = 1/beta;
sdfss				= beta;
deltauss			= delta0;
deltauprimess		= delta1;

model;

y + fixedcost 			= productionconstant*(z*n)^(1 - alpha)*(u*k(-1))^(alpha);
c + leverageratio*k/rr 	= w*n  + de + leverageratio*k(-1);
w/wss				=   (1-nss)/(1-n);
vf 					= (utilityconstant*a*(c^(eta)*(1 - n)^(1 - eta))^((1 - sigma)/thetavf)  + beta*expvf1sigma^(1/thetavf))^(thetavf/(1 - sigma));
expvf1sigma 			= vf(+1)^(1 - sigma);
w*n 				= (1 - alpha)*(y + fixedcost)/mu;
rrk*u*k(-1) 		= alpha*(y + fixedcost)/mu;
q*deltauprime*u*k(-1) = alpha*(y + fixedcost)/mu;
k 			= ((1 - deltau) - (phik/2)*(inv/k(-1) - delta0)^(2))*k(-1) + inv;
deltau		= delta0 + delta1*(u-1) + (delta2/2)*(u-1)^(2);
deltauprime	= delta1 + delta2*(u-1);
sdf = beta*(a/a(-1))*((c^(eta)*(1 - n)^(1 - eta))/(c(-1)^(eta)*(1 - n(-1))^(1 - eta)))^((1 - sigma)/thetavf)
			*(c(-1)/c)*(vf^(1 - sigma)/expvf1sigma(-1))^(1 - 1/thetavf);
1 = rr*sdf(+1);
1 = r*sdf(+1)*(pie(+1))^(-1);
1 = sdf(+1)*(de(+1) + pe(+1))/pe;
log(r) = (1 - rhor)*(log(rss) + rhopie*log(pie/piess) + rhoy*log(y/y(-1))) + rhor*log(r(-1));
de = y - w*n - inv - (phip/2)*(pie/piess - 1)^(2)*y - leverageratio*(k(-1) - k/rr);
1 = sdf(+1)*(u(+1)*rrk(+1) + 
     q(+1)*((1 - deltau(+1)) - (phik/2)*(inv(+1)/k - delta0)^(2) + phik*(inv(+1)/k - delta0)*(inv(+1)/k)))/q;
1/q = 1 - phik*(inv/k(-1) - delta0);
phip*(pie/piess - 1)*(pie/piess) = (1 - thetamu) + thetamu/mu + 
  sdf(+1)*phip*(pie(+1)/piess - 1)*(y(+1)/y)*(pie(+1)/piess);
profit = (mu - 1)*y - fixedcost;
expre = (de(+1) + pe(+1))/pe;
expre2 = (de(+1) + pe(+1))^(2)/pe^(2);
varexpre = expre2 - (expre)^(2);
a = (1 - rhoa)*ass + rhoa*a(-1) + vola(-1)*ea;
vola = rhovola*vola(-1) + (1 - rhovola)*volass + volvola*evola;
z = (1 - rhoz)*zss + rhoz*z(-1) + volzss*ez;
end;

initval;
a 	 		= 	ass;
c 	 		= 	css;
de 	 		= 	dess;
deltau 		= 	delta0;
deltauprime = 	delta1;
expre 	 	= 	express;
expre2 	 	= 	expre2ss;
expvf1sigma = 	expvf1sigmass;
inv 	 	= 	invss;
k 	 		= 	kss;
mu 	 		= 	muss;
n 	 		= 	nss;
pe 	 		= 	pess;
pie 	 	= 	piess;
profit 		= 	profitss;
q 	 		= 	qss;
r 	 		= 	rss;
rr 	 		= 	rrss;
rrk 	 	= 	rrkss;
sdf			= 	sdfss;
u			= 	uss;
varexpre 	= 	varexpress;
vf 	 		= 	vfss;
vola 	 	= 	volass;
w 	 		= 	wss;
y	 		= 	yss;
z			= 	zss;
end;

steady;

check;

shocks;
var ea 		= 1;
var evola 	= 1;
var ez 		= 1;
end;

stoch_simul(order=3,periods=0,irf=0,noprint,nograph,nomoments,nofunctions,nocorr,pruning);

irflength           = 20;
irfburninlength     = 200;

oo_.stochastic_steady_state = sss(oo_.dr, irfburninlength, options_.order);

for ivariable = 1:M_.endo_nbr
    assignin('base',[deblank(M_.endo_names(ivariable,:)) '_sss'],...
                     oo_.stochastic_steady_state(ivariable,:)');                     
end

for ishock = 1:M_.exo_nbr    

    variables_irfsss  =  irfsss(oo_.dr,M_.Sigma_e(:,ishock),irflength,irfburninlength,options_.order)';

    for ivariable = 1:M_.endo_nbr
        assignin('base',[deblank(M_.endo_names(ivariable,:)) '_' deblank(M_.exo_names(ishock,:))],...
                          variables_irfsss(:,ivariable));                         
    end

end

t         		= 1:irflength;

figure(1)
subplot(2,3,1)
hold on
plot(t,100*y_evola./y_sss,'--','Color','r','LineWidth',3)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
title('Output','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.3 0.2]);set(gca,'YTick',[-0.3:0.1:0.2],'FontSize',12);ylabel('Percent','FontSize',12)


figure(1)
subplot(2,3,2)
hold on
plot(t,100*c_evola./c_sss,'--','Color','r','LineWidth',3)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
title('Consumption','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.2 0.1]);set(gca,'YTick',[-0.2:0.1:0.1],'FontSize',12);ylabel('Percent','FontSize',12)

figure(1)
subplot(2,3,3)
hold on
plot(t,100*inv_evola./inv_sss,'--','Color','r','LineWidth',3)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
title('Investment','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-1 1]);set(gca,'YTick',[-1:0.5:1],'FontSize',12);ylabel('Percent','FontSize',12)


subplot(2,3,4)
hold on
plot(t,100*mu_evola./mu_sss,'--','Color','r','LineWidth',3)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
title('Markup','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.2 0.2]);set(gca,'YTick',[-0.2:0.1:0.2],'FontSize',12);ylabel('Percent','FontSize',12)

figure(1)
subplot(2,3,5)
hold on
plot(t,100*n_evola./n_sss,'--','Color','r','LineWidth',3)
plot(zeros(irflength,1),'k--','HandleVisibility','off'); xlim([0 irflength]);set(gca,'XTick',[4:4:irflength],'FontSize',12);
title('Hours Worked','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);
ylim([-0.2 0.2]);set(gca,'YTick',[-0.2:0.1:0.2],'FontSize',12);ylabel('Percent','FontSize',12)






