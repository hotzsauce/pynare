% Find indicies for variables of interest

M_.endo_names_cellstr       = cellstr(M_.endo_names);    
yindex                      = find(strcmp(M_.endo_names_cellstr,'y') == 1);
cindex                      = find(strcmp(M_.endo_names_cellstr,'c') == 1);
invindex                    = find(strcmp(M_.endo_names_cellstr,'inv') == 1);
nindex                      = find(strcmp(M_.endo_names_cellstr,'n') == 1);
varexpreindex               = find(strcmp(M_.endo_names_cellstr,'varexpre') == 1);
volaindex                   = find(strcmp(M_.endo_names_cellstr,'vola') == 1);
pieindex                    = find(strcmp(M_.endo_names_cellstr,'pie') == 1);
rindex                      = find(strcmp(M_.endo_names_cellstr,'r') == 1);
expreindex                  = find(strcmp(M_.endo_names_cellstr,'expre') == 1);


% Compute stochastic steady state

irflength           = 20;
irfburninlength     = 100;

oo_.stochastic_steady_state = sss(oo_.dr, irfburninlength, options_.order);

y_sss           = oo_.stochastic_steady_state(yindex);
c_sss           = oo_.stochastic_steady_state(cindex);
inv_sss         = oo_.stochastic_steady_state(invindex);
n_sss           = oo_.stochastic_steady_state(nindex);
varexpre_sss    = oo_.stochastic_steady_state(varexpreindex);
vola_sss        = oo_.stochastic_steady_state(volaindex);
pie_sss         = oo_.stochastic_steady_state(pieindex);
r_sss           = oo_.stochastic_steady_state(rindex);
expre_sss       = oo_.stochastic_steady_state(expreindex);


% Compute traditional (no other shock) impulse response from stochastic
%  steady state
 
variables_irfsss    =  irfsss(oo_.dr,M_.Sigma_e(:,2),irflength,irfburninlength,options_.order)';

y_evola           = variables_irfsss(:,yindex);
c_evola           = variables_irfsss(:,cindex);
inv_evola         = variables_irfsss(:,invindex);
n_evola           = variables_irfsss(:,nindex);
varexpre_evola    = variables_irfsss(:,varexpreindex);
vola_evola        = variables_irfsss(:,volaindex);
pie_evola         = variables_irfsss(:,pieindex);
r_evola           = variables_irfsss(:,rindex);
expre_evola       = variables_irfsss(:,expreindex);


y_evola         = y_evola./y_sss;
c_evola         = c_evola./c_sss;
inv_evola       = inv_evola./inv_sss;
n_evola         = n_evola./n_sss;
vola_evola      = vola_evola./vola_sss;
r_evola			= 400*r_evola;
pl_evola        = log( cumprod(pie_evola+pie_sss) ./ cumprod(pie_sss*ones(irflength,1)));
vxo_evola 		= 100*sqrt(4*(max(varexpre_evola+varexpre_sss,1e-16)));
vxo_sss         = 100*sqrt(4*(varexpre_sss));
logvxo_evola	= log(vxo_evola./ vxo_sss);



simulationlength = 25500;
ex_              = randn(simulationlength,3);

variables_simulation    =  simult_(oo_.stochastic_steady_state,oo_.dr,ex_,3);

y           = variables_simulation(yindex,:)';
c           = variables_simulation(cindex,:)';
inv         = variables_simulation(invindex,:)';
n           = variables_simulation(nindex,:)';

nvariables              = 4;
samplelength            = 120;
nsample                 = 100;
rollingwindowlength     = 20;

standarddeviations      = zeros(nvariables,nsample);
stddeviationrollingstddeviations = zeros(nvariables,nsample);


for isample = 1:nsample
    
    samplestart = 250+(isample-1)*250;
    sampleend   = samplestart + samplelength - 1;
    
    ypercentdeviation   = 100*(y(samplestart:sampleend)./mean(y(samplestart:sampleend))-1);
    cpercentdeviation   = 100*(c(samplestart:sampleend)./mean(c(samplestart:sampleend))-1);
    invpercentdeviation = 100*(inv(samplestart:sampleend)./mean(inv(samplestart:sampleend))-1);
    npercentdeviation   = 100*(n(samplestart:sampleend)./mean(n(samplestart:sampleend))-1);

    simulateddata       = [ypercentdeviation cpercentdeviation invpercentdeviation npercentdeviation];
    
%     % Add pre-sample data for rolling standard deviation calculations
%     
    ypercentdeviationrolling   = 100*(y(samplestart-rollingwindowlength-1:sampleend)./mean(y(samplestart-rollingwindowlength-1:sampleend))-1);
    cpercentdeviationrolling   = 100*(c(samplestart-rollingwindowlength-1:sampleend)./mean(c(samplestart-rollingwindowlength-1:sampleend))-1);
    invpercentdeviationrolling = 100*(inv(samplestart-rollingwindowlength-1:sampleend)./mean(inv(samplestart-rollingwindowlength-1:sampleend))-1);
    npercentdeviationrolling   = 100*(n(samplestart-rollingwindowlength-1:sampleend)./mean(n(samplestart-rollingwindowlength-1:sampleend))-1);
%     
    simulateddatarolling       = [ypercentdeviationrolling cpercentdeviationrolling invpercentdeviationrolling npercentdeviationrolling];
%     
    standarddeviations(:,isample)   = std(simulateddata)';
%     
    rollingstddeviations    = zeros(samplelength,nvariables);
% 
    for idata = 1:nvariables;
        for irolling = 1:samplelength
            rollingstddeviations(irolling,idata) = std(simulateddatarolling(irolling:irolling+rollingwindowlength-1,idata));
        end
    end
% 
    stddeviationrollingstddeviations(:,isample) = std(rollingstddeviations)';

end

stdymodel = mean(standarddeviations(1,:));
stdcmodel = mean(standarddeviations(2,:));
stdinvmodel = mean(standarddeviations(3,:));
stdnmodel = mean(standarddeviations(4,:));

svymodel = mean(stddeviationrollingstddeviations(1,:));
svcmodel = mean(stddeviationrollingstddeviations(2,:));
svinvmodel = mean(stddeviationrollingstddeviations(3,:));
svnmodel = mean(stddeviationrollingstddeviations(4,:));







stochasticvolatilitymodel = mean(stddeviationrollingstddeviations')';

