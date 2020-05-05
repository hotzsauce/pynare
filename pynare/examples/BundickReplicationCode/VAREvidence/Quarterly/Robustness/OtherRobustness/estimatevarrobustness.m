clear
close all

estimatebaselinevarnoci


estimatenipacivarnoci

cd ../../../Monthly/BaselineMonthly/

estimatebaselinevarmonthlynoci

cd ../RobustnessMonthly/

estimatebaselinevarmonthlynocivxolast

cd ../../Quarterly/Robustness/OtherRobustness

subplot(2,3,5)
legend(sprintf('Baseline VAR Model'),sprintf('NIPA Consumption \n& Investment Definitions'),...
    sprintf('Monthly Frequency Estimation\nQuarterly Aggregated IRF\nSpliced Output Series'), ...
    sprintf('Monthly Frequency Estimation\nQuarterly Aggregated IRF\nVXO Ordered Last'),'Location','NorthEast')

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [-0.75 -.65 12.0 9.25]);
print('-dpdf',['BaselineVARRobustnessNIPAMonthly.pdf'])

clear
close all

estimatebaselinevarnoci

estimate5yearvarnoci

estimateprezlbvarnoci

estimatebloomsamplevarnoci

subplot(2,3,5)
legend(sprintf('Baseline VAR Model'),sprintf('5-Year Treasury'),...
    sprintf('1986Q1 - 2009Q4 Sample\nFederal Funds Rate'),sprintf('Bloom (2009)\n1962Q3 - 2008Q2 Sample\nSpliced Volatility Series'),'Location','NorthEast')

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [-0.75 -.65 12.0 9.25]);
print('-dpdf',['BaselineVARRobustnessPolicy.pdf'])




