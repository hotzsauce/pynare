

clear
close all
clc


dynare bbeffectivedemandoptimalpolicy.mod noclearall

!rm bbeffectivedemandoptimalpolicy_*
!rm bbeffectivedemandoptimalpolicy.log
!rm bbeffectivedemandoptimalpolicy.jnl
!rm bbeffectivedemandoptimalpolicy.m

cd ZeroLowerBoundModel

clear; postsolution;

sssirfvola; sssirfvolazlb;

cd ../



