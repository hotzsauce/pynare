
nvariables              = 7;
samplelength            = 120;
nsample                 = 10000;
rollingwindowlength     = 20;

means                   = zeros(nvariables,nsample);
standarddeviations      = zeros(nvariables,nsample);
stddeviationrollingstddeviations = zeros(nvariables,nsample);


for isample = 1:nsample

    samplestart = 250+(isample-1)*250;
    sampleend   = samplestart + samplelength - 1;
    
    vxolevel            = 100*sqrt(4*(max(varexpre(samplestart:sampleend),1e-16)));

    ypercentdeviation   = 100*(y(samplestart:sampleend)./mean(y(samplestart:sampleend))-1);
    cpercentdeviation   = 100*(c(samplestart:sampleend)./mean(c(samplestart:sampleend))-1);
    invpercentdeviation = 100*(inv(samplestart:sampleend)./mean(inv(samplestart:sampleend))-1);
    npercentdeviation   = 100*(n(samplestart:sampleend)./mean(n(samplestart:sampleend))-1);
    rrlevel             = 400*(r(samplestart:sampleend)-pie(samplestart:sampleend));
    equitypremiumlevel  = 400*(expostexpre(samplestart:sampleend) - (1+r(samplestart:sampleend)-pie(samplestart:sampleend)));
    
    simulateddata       = [ypercentdeviation cpercentdeviation invpercentdeviation npercentdeviation rrlevel equitypremiumlevel vxolevel];
    
    % Add pre-sample data for rolling standard deviation calculations
    
    vxolevelrolling            = 100*sqrt(4*(max(varexpre(samplestart-rollingwindowlength-1:sampleend),1e-16)));
    
    ypercentdeviationrolling   = 100*(y(samplestart-rollingwindowlength-1:sampleend)./mean(y(samplestart-rollingwindowlength-1:sampleend))-1);
    cpercentdeviationrolling   = 100*(c(samplestart-rollingwindowlength-1:sampleend)./mean(c(samplestart-rollingwindowlength-1:sampleend))-1);
    invpercentdeviationrolling = 100*(inv(samplestart-rollingwindowlength-1:sampleend)./mean(inv(samplestart-rollingwindowlength-1:sampleend))-1);
    npercentdeviationrolling   = 100*(n(samplestart-rollingwindowlength-1:sampleend)./mean(n(samplestart-rollingwindowlength-1:sampleend))-1);
    rrlevelrolling             = 400*(r(samplestart-rollingwindowlength-1:sampleend)-pie(samplestart-rollingwindowlength-1:sampleend));
    equitypremiumlevelrolling  = 400*(expostexpre(samplestart-rollingwindowlength-1:sampleend) - (1+r(samplestart-rollingwindowlength-1:sampleend)-pie(samplestart-rollingwindowlength-1:sampleend)));
    
    simulateddatarolling       = [ypercentdeviationrolling cpercentdeviationrolling invpercentdeviationrolling npercentdeviationrolling rrlevelrolling equitypremiumlevelrolling vxolevelrolling];
    
    means(:,isample)                = mean(simulateddata)';

    standarddeviations(:,isample)   = std(simulateddata)';
    
    rollingstddeviations    = zeros(samplelength,nvariables);

    for idata = 1:nvariables;
        for irolling = 1:samplelength
            rollingstddeviations(irolling,idata) = std(simulateddatarolling(irolling:irolling+rollingwindowlength-1,idata));
        end
    end

    stddeviationrollingstddeviations(:,isample) = std(rollingstddeviations)';

end

smallsamplemeansmean = mean(means')';
smallsamplemeanslowerci = prctile(means',2.5)';
smallsamplemeansupperci = prctile(means',97.5)';

disp('  ')

disp([char('Variable','Y:   ','C:    ','I: ','N:   ','Ex-Post RR:   ','Ex-Post Equity Premium:   ','VXO:   ') strvcat('    Mean   ',num2str(smallsamplemeansmean)) ...
     char('  ','  ','  ','  ','  ','  ','  ',' ') strvcat('95% Lower CI   ',num2str(smallsamplemeanslowerci)) ...
     char('  ','  ','  ','  ','  ','  ','  ',' ') strvcat(' 95% Upper CI   ',num2str(smallsamplemeansupperci))])

smallsamplestdsmean = mean(standarddeviations')';
smallsamplestdslowerci = prctile(standarddeviations',2.5)';
smallsamplestdsupperci = prctile(standarddeviations',97.5)';

disp('  ')

disp([char('Variable','Y:   ','C:    ','I: ','N:   ','Ex-Post RR:   ','Ex-Post Equity Premium:   ','VXO:   ') strvcat('Standard Deviation',num2str(smallsamplestdsmean)) ...
     char('  ','  ','  ','  ','  ','  ','  ',' ') strvcat('95% Lower CI   ',num2str(smallsamplestdslowerci)) ...
     char('  ','  ','  ','  ','  ','  ','  ',' ') strvcat(' 95% Upper CI   ',num2str(smallsamplestdsupperci))])


disp('  ')

 
smallsamplestdsmean = mean(stddeviationrollingstddeviations')';
smallsamplestdslowerci = prctile(stddeviationrollingstddeviations',2.5)';
smallsamplestdsupperci = prctile(stddeviationrollingstddeviations',97.5)';

disp([char('Variable','Y:   ','C:    ','I: ','N:   ','Ex-Post RR:   ','Ex-Post Equity Premium:   ','VXO:   ') strvcat('Std Dev(Rolling Std Dev)',num2str(smallsamplestdsmean)) ...
     char('  ','  ','  ','  ','  ','  ','  ',' ') strvcat('95% Lower CI   ',num2str(smallsamplestdslowerci)) ...
     char('  ','  ','  ','  ','  ','  ','  ',' ') strvcat(' 95% Upper CI   ',num2str(smallsamplestdsupperci))])

disp('  ')
 
 