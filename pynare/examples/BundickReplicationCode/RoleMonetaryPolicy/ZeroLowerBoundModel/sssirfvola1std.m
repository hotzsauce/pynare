
nperiods = 2021;
shockperiod = 2001;

yshock    			= zeros(nperiods,1);
ydatashock 			= zeros(nperiods,1);
pieshock  			= zeros(nperiods,1);
rshock   		    = zeros(nperiods,1);
rdshock   		    = zeros(nperiods,1);
lagyshock          = zeros(nperiods,1);
cshock  	 	    = zeros(nperiods,1);
ishock				= zeros(nperiods,1);
ushock				= zeros(nperiods,1);
deltaushock			= zeros(nperiods,1);
deltauprimeshock	= zeros(nperiods,1);
qshock				= zeros(nperiods,1);
kshock				= zeros(nperiods,1);
k1shock				= zeros(nperiods,1);
nshock	       		= zeros(nperiods,1);
wshock		  		= zeros(nperiods,1);
rrkshock		  	= zeros(nperiods,1);
mushock		    	= zeros(nperiods,1);
vfshock		    	= zeros(nperiods,1);
expvf1sigmashock	= zeros(nperiods,1);
ashock      		= zeros(nperiods,1);
lagvolashock   		= zeros(nperiods,1);
volashock   		= zeros(nperiods,1);
epsilonashock		= zeros(nperiods,1);
epsilonvolashock	= zeros(nperiods,1);

ybaseline    		= zeros(nperiods,1);
ydatabaseline   	= zeros(nperiods,1);
piebaseline  		= zeros(nperiods,1);
rbaseline   	    = zeros(nperiods,1);
rdbaseline   		= zeros(nperiods,1);
lagybaseline       = zeros(nperiods,1);
cbaseline  	 	    = zeros(nperiods,1);
ibaseline			= zeros(nperiods,1);
ubaseline			= zeros(nperiods,1);
deltaubaseline		= zeros(nperiods,1);
deltauprimebaseline	= zeros(nperiods,1);
qbaseline			= zeros(nperiods,1);
kbaseline			= zeros(nperiods,1);
k1baseline			= zeros(nperiods,1);
nbaseline           = zeros(nperiods,1);
wbaseline		  	= zeros(nperiods,1);
rrkbaseline		  	= zeros(nperiods,1);
mubaseline		    = zeros(nperiods,1);
vfbaseline		    = zeros(nperiods,1);
expvf1sigmabaseline	= zeros(nperiods,1);
abaseline      		= zeros(nperiods,1);
lagvolabaseline   	= zeros(nperiods,1);
volabaseline   		= zeros(nperiods,1);
epsilonabaseline	= zeros(nperiods,1);
epsilonvolabaseline	= zeros(nperiods,1);


for i = 1:nperiods

   if i == 1
        
        epsilonashock(i,1) = 0;
        epsilonvolashock(i,1) = 0;
        volashock(i,1) = volass;
        ashock(i,1) = ass;
        lagyshock(i,1) = yss;
        lagvolashock(i,1) = volass;
        kshock(i,1) = kss;
        
        
        epsilonabaseline(i,1) = 0;
        epsilonvolabaseline(i,1) = 0;
        volabaseline(i,1) = volass;
        abaseline(i,1) = ass;
        lagybaseline(i,1) = yss;
        lagvolabaseline(i,1) = volass;
        kbaseline(i,1) = kss;
        
        
    else
        
        if i == shockperiod
			
        	epsilonashock(i,1) = 0;%/(100*lagvolashock(i,1));
            ashock(i,1) = rhoa*ashock(i-1,1) + (1-rhoa)*ass + lagvolashock(i,1)*epsilonashock(i,1);
            
            epsilonvolashock(i,1) = 1;
			volashock(i,1) = rhovola*lagvolashock(i,1) + (1-rhovola)*volass + volvola*epsilonvolashock(i,1);
            
        	epsilonabaseline(i,1) = 0/(100*lagvolashock(i,1));
			abaseline(i,1) = rhoa*abaseline(i-1,1) + (1-rhoa)*ass + lagvolabaseline(i,1)*epsilonabaseline(i,1);
            
            epsilonvolabaseline(i,1) = 0;
			volabaseline(i,1) = rhovola*lagvolabaseline(i,1) + (1-rhovola)*volass + volvola*epsilonvolabaseline(i,1);
            
            volashock(i,1) = max(0.00005,volashock(i,1));
            volabaseline(i,1) = max(0.00005,volabaseline(i,1));
            
        else
            
            epsilonashock(i,1) = 0;
   			ashock(i,1) = rhoa*ashock(i-1,1) + (1-rhoa)*ass + lagvolashock(i,1)*epsilonashock(i,1);
            
            epsilonvolashock(i,1) = 0;
			volashock(i,1) = rhovola*lagvolashock(i,1) + (1-rhovola)*volass + volvola*epsilonvolashock(i,1);
            
        	epsilonabaseline(i,1) = 0;
			abaseline(i,1) = rhoa*abaseline(i-1,1) + (1-rhoa)*ass + lagvolabaseline(i,1)*epsilonabaseline(i,1);
            
            epsilonvolabaseline(i,1) = 0;
			volabaseline(i,1) = rhovola*lagvolabaseline(i,1) + (1-rhovola)*volass + volvola*epsilonvolabaseline(i,1);

            volashock(i,1) = max(0.00005,volashock(i,1));
            volabaseline(i,1) = max(0.00005,volabaseline(i,1));

        end
        
    end
    
   
    nshock(i,1) 			= linearinterp4(na, grida, stepa, ashock(i,1), nk, gridk, stepk, kshock(i,1), nvola, gridvola, stepvola, volashock(i,1), nlagy, gridlagy, steplagy, lagyshock(i,1), rule(:,:,:,:,1));
    ishock(i,1) 			= linearinterp4(na, grida, stepa, ashock(i,1), nk, gridk, stepk, kshock(i,1), nvola, gridvola, stepvola, volashock(i,1), nlagy, gridlagy, steplagy, lagyshock(i,1), rule(:,:,:,:,2));
    pieshock(i,1) 			= linearinterp4(na, grida, stepa, ashock(i,1), nk, gridk, stepk, kshock(i,1), nvola, gridvola, stepvola, volashock(i,1), nlagy, gridlagy, steplagy, lagyshock(i,1), rule(:,:,:,:,3));
    ushock(i,1) 			= linearinterp4(na, grida, stepa, ashock(i,1), nk, gridk, stepk, kshock(i,1), nvola, gridvola, stepvola, volashock(i,1), nlagy, gridlagy, steplagy, lagyshock(i,1), rule(:,:,:,:,4));
    expvf1sigmashock(i,1)	= linearinterp4(na, grida, stepa, ashock(i,1), nk, gridk, stepk, kshock(i,1), nvola, gridvola, stepvola, volashock(i,1), nlagy, gridlagy, steplagy, lagyshock(i,1), rule(:,:,:,:,5));
    
    
    yshock(i,1) 			= productionconstant*(ushock(i,1)*kshock(i,1))^(alpha)*nshock(i,1)^(1-alpha) - fixedcost;
    deltaushock(i,1)		= delta0 + delta1*(ushock(i,1)-1) + (delta2/2)*(ushock(i,1)-1D0)^(2D0);
    deltauprimeshock(i,1)	= delta1 + delta2*(ushock(i,1)-1);
	k1shock(i,1)			= (1 - deltaushock(i,1) - (phik/2)*(ishock(i,1)/kshock(i,1)-delta0)^(2))*kshock(i,1) + ishock(i,1);
	qshock(i,1)				= (1 - phik*(ishock(i,1)/kshock(i,1) - delta0))^(-1);
	cshock(i,1)				= yshock(i,1) - ishock(i,1) - (phip/2)*(pieshock(i,1)/piess-1)^(2)*yshock(i,1);
	wshock(i,1)				= ((1-eta) / eta)*cshock(i,1)/(1-nshock(i,1));
	mushock(i,1)			= (1-alpha)*(yshock(i,1) + fixedcost)/(wshock(i,1)*nshock(i,1));
	rrkshock(i,1)			= alpha*(yshock(i,1) + fixedcost)/(ushock(i,1)*kshock(i,1)*mushock(i,1));
	rdshock(i,1)			= exp(log(rss) + rhopie*log(pieshock(i,1)/piess) + rhoy*log(yshock(i,1)/lagyshock(i,1)));
    rshock(i,1)             = max(1,rdshock(i,1));
	vfshock(i,1)			= (utilityconstant*ashock(i,1)*(cshock(i,1)^(eta)*(1-nshock(i,1))^(1-eta))^((1-sigma)/thetavf) ...
								+ beta * expvf1sigmashock(i,1)^(1/thetavf))^(thetavf/(1-sigma));
    ydatashock(i,1) 		= cshock(i,1) + ishock(i,1);
                            
    
    nbaseline(i,1)              = linearinterp4(na, grida, stepa, abaseline(i,1), nk, gridk, stepk, kbaseline(i,1), nvola, gridvola, stepvola, volabaseline(i,1), nlagy, gridlagy, steplagy, lagybaseline(i,1), rule(:,:,:,:,1));
    ibaseline(i,1)              = linearinterp4(na, grida, stepa, abaseline(i,1), nk, gridk, stepk, kbaseline(i,1), nvola, gridvola, stepvola, volabaseline(i,1), nlagy, gridlagy, steplagy, lagybaseline(i,1), rule(:,:,:,:,2));
    piebaseline(i,1)            = linearinterp4(na, grida, stepa, abaseline(i,1), nk, gridk, stepk, kbaseline(i,1), nvola, gridvola, stepvola, volabaseline(i,1), nlagy, gridlagy, steplagy, lagybaseline(i,1), rule(:,:,:,:,3));
    ubaseline(i,1)              = linearinterp4(na, grida, stepa, abaseline(i,1), nk, gridk, stepk, kbaseline(i,1), nvola, gridvola, stepvola, volabaseline(i,1), nlagy, gridlagy, steplagy, lagybaseline(i,1), rule(:,:,:,:,4));
    expvf1sigmabaseline(i,1)    = linearinterp4(na, grida, stepa, abaseline(i,1), nk, gridk, stepk, kbaseline(i,1), nvola, gridvola, stepvola, volabaseline(i,1), nlagy, gridlagy, steplagy, lagybaseline(i,1), rule(:,:,:,:,5));
   
    ybaseline(i,1)              = productionconstant*(ubaseline(i,1)*kbaseline(i,1))^(alpha)*nbaseline(i,1)^(1-alpha) - fixedcost;
    deltaubaseline(i,1)         = delta0 + delta1*(ubaseline(i,1)-1) + (delta2/2)*(ubaseline(i,1)-1D0)^(2D0);
    deltauprimebaseline(i,1)	= delta1 + delta2*(ubaseline(i,1)-1);
	k1baseline(i,1)             = (1 - deltaubaseline(i,1) - (phik/2)*(ibaseline(i,1)/kbaseline(i,1)-delta0)^(2))*kbaseline(i,1) + ibaseline(i,1);
	qbaseline(i,1)				= (1 - phik*(ibaseline(i,1)/kbaseline(i,1) - delta0))^(-1);
	cbaseline(i,1)				= ybaseline(i,1) - ibaseline(i,1) - (phip/2)*(piebaseline(i,1)/piess-1)^(2)*ybaseline(i,1);
	wbaseline(i,1)				= ((1-eta) / eta)*cbaseline(i,1)/(1-nbaseline(i,1));
	mubaseline(i,1)             = (1-alpha)*(ybaseline(i,1) + fixedcost)/(wbaseline(i,1)*nbaseline(i,1));
	rrkbaseline(i,1)			= alpha*(ybaseline(i,1) + fixedcost)/(ubaseline(i,1)*kbaseline(i,1)*mubaseline(i,1));
	rdbaseline(i,1)             = exp(log(rss) + rhopie*log(piebaseline(i,1)/piess) + rhoy*log(ybaseline(i,1)/lagybaseline(i,1)));
    rbaseline(i,1)              = max(1,rdbaseline(i,1));
	vfbaseline(i,1)             = (utilityconstant*abaseline(i,1)*(cbaseline(i,1)^(eta)*(1-nbaseline(i,1))^(1-eta))^((1-sigma)/thetavf) ...
								+ beta * expvf1sigmabaseline(i,1)^(1/thetavf))^(thetavf/(1-sigma));
    ydatabaseline(i,1)          = cbaseline(i,1) + ibaseline(i,1);
                            
                            
    if i < nperiods

       kshock(i+1,1) = k1shock(i,1);
       kbaseline(i+1,1) = k1baseline(i,1);
 
       lagvolashock(i+1,1) = volashock(i,1);
       lagvolabaseline(i+1,1) = volabaseline(i,1);
       
       lagyshock(i+1,1) = yshock(i,1);
       lagybaseline(i+1,1) = ybaseline(i,1);
       
    end
    
end

ystochasticss            = yshock(shockperiod-1,1);
ydatastochasticss        = ydatashock(shockperiod-1,1);
piestochasticss          = pieshock(shockperiod-1,1);
rdstochasticss           = rdshock(shockperiod-1,1);
rstochasticss            = rshock(shockperiod-1,1);
cstochasticss            = cshock(shockperiod-1,1);
istochasticss            = ishock(shockperiod-1,1);
ustochasticss            = ushock(shockperiod-1,1);
nstochasticss            = nshock(shockperiod-1,1);
wstochasticss            = wshock(shockperiod-1,1);
rrkstochasticss          = rrkshock(shockperiod-1,1);
kstochasticss            = kshock(shockperiod-1,1);
mustochasticss           = mushock(shockperiod-1,1);
vfstochasticss           = vfshock(shockperiod-1,1);
expvf1sigmastochasticss  = expvf1sigmashock(shockperiod-1,1);
astochasticss            = ashock(shockperiod-1,1);
volastochasticss         = volashock(shockperiod-1,1);

disp('-----------------------------------------------------------------------------------')
disp(' ')
disp('Stochastic Steady State')
disp(' ')
disp(['Ouptut                    ',num2str(ystochasticss)])
disp(['Data Consistent Ouptut    ',num2str(ydatastochasticss)])
disp(['Inflation                 ',num2str(piestochasticss)])
disp(['Nominal Interest Rate     ',num2str(rstochasticss)])
disp(['Consumption               ',num2str(cstochasticss)])
disp(['Investment                ',num2str(istochasticss)])
disp(['Utilization               ',num2str(ustochasticss)])
disp(['Capital                   ',num2str(kstochasticss)])
disp(['Hours Worked              ',num2str(nstochasticss)])
disp(['Real Wage                 ',num2str(wstochasticss)])
disp(['Rental Rate               ',num2str(rrkstochasticss)])
disp(['Markup                    ',num2str(mustochasticss)])
disp(['Discount Factor Shock     ',num2str(astochasticss)])
disp(['Discount Factor Volatiliy ',num2str(volastochasticss)])

yshockbaselinedeviation = 100*log(yshock./ybaseline);
ydatashockbaselinedeviation = 100*log(ydatashock./ydatabaseline);
pieshockbaselinedeviation = 400*log(pieshock./piebaseline);
rshockbaselinedeviation = 400*log(rshock./rbaseline);
cshockbaselinedeviation = 100*log(cshock./cbaseline);
ishockbaselinedeviation = 100*log(ishock./ibaseline);
ushockbaselinedeviation = 100*log(ushock./ubaseline);
qshockbaselinedeviation = 100*log(qshock./qbaseline);
kshockbaselinedeviation = 100*log(kshock./kbaseline);
rrkshockbaselinedeviation = 100*log(rrkshock./rrkbaseline);
nshockbaselinedeviation = 100*log(nshock./nbaseline);
wshockbaselinedeviation = 100*log(wshock./wbaseline);
mushockbaselinedeviation = 100*log(mushock./mubaseline);
ashockbaselinedeviation = 100*log(ashock./abaseline);
volashockbaselinedeviation = 100*(volashock./volabaseline-1);

t = 1:(nperiods - shockperiod + 1);

figure(1)
subplot(3,3,1);hold on;
plot(t,yshockbaselinedeviation(shockperiod:nperiods,1),'r--','LineWidth',3)
plot(t,zeros(nperiods-shockperiod+1,1),'k--','HandleVisibility','off');
title('Output','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,2);hold on;
plot(t,cshockbaselinedeviation(shockperiod:nperiods,1),'r--','LineWidth',3)
plot(t,zeros(nperiods-shockperiod+1,1),'k--','HandleVisibility','off');
title('Consumption','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,3);hold on;
plot(t,ishockbaselinedeviation(shockperiod:nperiods,1),'r--','LineWidth',3)
plot(t,zeros(nperiods-shockperiod+1,1),'k--','HandleVisibility','off');
title('Investment','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,4);hold on;
plot(t,mushockbaselinedeviation(shockperiod:nperiods,1),'r--','LineWidth',3)
plot(t,zeros(nperiods-shockperiod+1,1),'k--','HandleVisibility','off');
title('Markup','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,5);hold on;
plot(t,nshockbaselinedeviation(shockperiod:nperiods,1),'r--','LineWidth',3)
title('Hours Worked','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Percent','FontSize',12)

subplot(3,3,8);hold on;
plot(t,pieshockbaselinedeviation(shockperiod:nperiods,1),'r--','LineWidth',3)
title('Inflation','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Percent','FontSize',12)


subplot(3,3,6);hold on;
plot(t,400*log(rshock(shockperiod:nperiods,1)),'r--','LineWidth',3)
title('Nominal Interest Rate','FontSize',16);xlim([t(1) t(end)]);set(gca,'XTick',t(4:4:end),'FontSize',12);ylabel('Annualized Percent - Level','FontSize',12)

