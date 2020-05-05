clear

alpha           = 0.333;
beta            = 0.994;
delta0          = 0.025;
delta1          = (1/beta - 1 + delta0);
delta2          = 0.00031;
phip            = 100;
phik            = 10;
sigma           = 15;
ies             = 0.95;
frischlse       = 2.0;
thetavf         = (1 - sigma)/(1 - 1/ies);
thetamu         = 6.0;

piess           = 1.005;
rss             = piess/beta;
rhopie          = 1.5;
rhoy            = 0.2;

ass             = 1;
rhoa            = 0.93564;
volass          = 0.0026251;
rhovola         = 0.74227;
volvola         = 0.0025022;

yss         = 1.0;
rrkss       = 1/beta - 1 + delta0;
muss        = thetamu/(thetamu - 1);
uss			= 1;
fixedcost   = (muss - 1)*yss;
kss         = alpha*(yss + fixedcost)/(rrkss*muss);
iss         = delta0*kss;
css         = yss - iss;

nsssolution = fsolve(@(nssguess) 1+(-1+1/ies)/(1+((-1+alpha)*(-1+nssguess)*(fixedcost+yss))/(css*nssguess*muss)) -(frischlse/ies)*nssguess*(1-nssguess)^(-1),0.5,optimoptions('fsolve','Display','off'));
nss         = nsssolution;

wss                 = (1 - alpha)*(yss + fixedcost)/(nss*muss);
eta                 = (wss*(1 - nss)/css + 1)^(-1);
qss                 = 1;
productionconstant 	= (yss + fixedcost)/(kss^(alpha)*(nss)^(1 - alpha));
vfss                = 1;
utilityconstant     = vfss^((1 - sigma)/thetavf)*(1 - beta)/(css^(eta)*(1 - nss)^(1 - eta))^((1 - sigma)/thetavf);
dess                = css - wss*nss;
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

na = 14; 
mina = 0.90;
maxa = 1.03;
stepa = (maxa - mina) / (na - 1);
grida = linspace(mina,maxa,na)';

nk = 17; 
mink = 10.4;
maxk = 11.2;
stepk = (maxk - mink) / (nk - 1);
gridk = linspace(mink,maxk,nk)';

nvola = 10; 
minvola = 0.001;
maxvola = 0.00775;
stepvola = (maxvola - minvola) / (nvola - 1);
gridvola = linspace(minvola,maxvola,nvola)';

nlagy = 13; 
minlagy = 0.97;
maxlagy = 1.03;
steplagy = (maxlagy - minlagy) / (nlagy - 1);
gridlagy = linspace(minlagy,maxlagy,nlagy)';

nepsilona = 17;
rangestdepsilona = 4;
minepsilona = -1*rangestdepsilona;
maxepsilona = 1*rangestdepsilona;
stepepsilona = (maxepsilona - minepsilona) / (nepsilona - 1);
gridepsilona = linspace(minepsilona,maxepsilona,nepsilona)';

nepsilonvola = 17;
rangestdepsilonvola = 4;
minepsilonvola = -1*rangestdepsilonvola;
maxepsilonvola = 1*rangestdepsilonvola;
stepepsilonvola = (maxepsilonvola - minepsilonvola) / (nepsilonvola - 1);
gridepsilonvola = linspace(minepsilonvola,maxepsilonvola,nepsilonvola)';
