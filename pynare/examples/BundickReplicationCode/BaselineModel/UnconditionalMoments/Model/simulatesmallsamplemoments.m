
clear
clc
close all

dynare bbeffectivedemand.mod noclearall

disp('Baseline Model')

calculatemoments

dynare bbeffectivedemandhigherleverage.mod noclearall

disp('Higher Leverage')

calculatemoments

