

clear
close all
clc

cd StickyvsFlexiblePrices

dynare bbeffectivedemand.mod noclearall

!rm bbeffectivedemand_*
!rm bbeffectivedemand.log
!rm bbeffectivedemand.jnl
!rm bbeffectivedemand.m

dynare bbeffectivedemandflexibleprices.mod noclearall

!rm bbeffectivedemandflexibleprices_*
!rm bbeffectivedemandflexibleprices.log
!rm bbeffectivedemandflexibleprices.jnl
!rm bbeffectivedemandflexibleprices.m

cd ../

clc
close all
clear

cd ModelSupportVARIdentification

dynare bbeffectivedemandmodelsupportvar.mod noclearall

!rm bbeffectivedemandmodelsupportvar_*
!rm bbeffectivedemandmodelsupportvar.log
!rm bbeffectivedemandmodelsupportvar.jnl
!rm bbeffectivedemandmodelsupportvar.m

cd ../

clc
close all
clear

cd EstimateVARSimulatedData

simulateestimatevarsmallsamples

cd ../

cd ModelFeatures

cd TechnologyUncertaintyShock&RealMUs

dynare bbeffectivedemand.mod noclearall

!rm bbeffectivedemand_*
!rm bbeffectivedemand.log
!rm bbeffectivedemand.jnl
!rm bbeffectivedemand.m

dynare bbeffectivedemandrealmus.mod noclearall

!rm bbeffectivedemandrealmus_*
!rm bbeffectivedemandrealmus.log
!rm bbeffectivedemandrealmus.jnl
!rm bbeffectivedemandrealmus.m


dynare bbeffectivedemandtechnologyuncertainty.mod noclearall

!rm bbeffectivedemandtechnologyuncertainty_*
!rm bbeffectivedemandtechnologyuncertainty.log
!rm bbeffectivedemandtechnologyuncertainty.jnl
!rm bbeffectivedemandtechnologyuncertainty.m


cd ../../

clc
close all
clear

cd ModelFeatures

cd EpsteinZin

dynare bbeffectivedemand.mod noclearall

!rm bbeffectivedemand_*
!rm bbeffectivedemand.log
!rm bbeffectivedemand.jnl
!rm bbeffectivedemand.m

dynare bbeffectivedemandlowersigma.mod noclearall

!rm bbeffectivedemandlowersigma_*
!rm bbeffectivedemandlowersigma.log
!rm bbeffectivedemandlowersigma.jnl
!rm bbeffectivedemandlowersigma.m

dynare bbeffectivedemandlowersigmalargershocks.mod noclearall

!rm bbeffectivedemandlowersigmalargershocks_*
!rm bbeffectivedemandlowersigmalargershocks.log
!rm bbeffectivedemandlowersigmalargershocks.jnl
!rm bbeffectivedemandlowersigmalargershocks.m

cd ../../

clear
clc
close all

cd ModelFeatures

cd InvestmentAdjustmentCosts

dynare bbeffectivedemand.mod noclearall

!rm bbeffectivedemand_*
!rm bbeffectivedemand.log
!rm bbeffectivedemand.jnl
!rm bbeffectivedemand.m

dynare bbeffectivedemandnoinvcosts.mod noclearall

!rm bbeffectivedemandnoinvcosts_*
!rm bbeffectivedemandnoinvcosts.log
!rm bbeffectivedemandnoinvcosts.jnl
!rm bbeffectivedemandnoinvcosts.m

dynare bbeffectivedemandinfiniteinvcosts.mod noclearall

!rm bbeffectivedemandinfiniteinvcosts_*
!rm bbeffectivedemandinfiniteinvcosts.log
!rm bbeffectivedemandinfiniteinvcosts.jnl
!rm bbeffectivedemandinfiniteinvcosts.m

dynare bbeffectivedemandinfiniteinvcostsflex.mod noclearall

!rm bbeffectivedemandinfiniteinvcostsflex_*
!rm bbeffectivedemandinfiniteinvcostsflex.log
!rm bbeffectivedemandinfiniteinvcostsflex.jnl
!rm bbeffectivedemandinfiniteinvcostsflex.m

cd ../../

clear
clc
close all

cd ModelFeatures

cd LeverageOtherShocksVolaPersistence

dynare bbeffectivedemandhigherleverage.mod noclearall

!rm bbeffectivedemandhigherleverage_*
!rm bbeffectivedemandhigherleverage.log
!rm bbeffectivedemandhigherleverage.jnl
!rm bbeffectivedemandhigherleverage.m

dynare bbeffectivedemandhighervolzss.mod noclearall

!rm bbeffectivedemandhighervolzss_*
!rm bbeffectivedemandhighervolzss.log
!rm bbeffectivedemandhighervolzss.jnl
!rm bbeffectivedemandhighervolzss.m

dynare bbeffectivedemandiidvola.mod noclearall

!rm bbeffectivedemandiidvola_*
!rm bbeffectivedemandiidvola.log
!rm bbeffectivedemandiidvola.jnl
!rm bbeffectivedemandiidvola.m

dynare bbeffectivedemand.mod noclearall

!rm bbeffectivedemand_*
!rm bbeffectivedemand.log
!rm bbeffectivedemand.jnl
!rm bbeffectivedemand.m


cd ../../

clc
close all
clear

cd ModelFeatures

cd Utilization

dynare bbeffectivedemand.mod noclearall

!rm bbeffectivedemand_*
!rm bbeffectivedemand.log
!rm bbeffectivedemand.jnl
!rm bbeffectivedemand.m


dynare bbeffectivedemandhigherdelta2.mod noclearall

!rm bbeffectivedemandhigherdelta2_*
!rm bbeffectivedemandhigherdelta2.log
!rm bbeffectivedemandhigherdelta2.jnl
!rm bbeffectivedemandhigherdelta2.m

cd ../../



clc
close all
clear

cd ModelFeatures

cd GHHLaborSupplyElasticity


dynare bbeffectivedemandghhflexibleprices.mod noclearall

!rm bbeffectivedemandghhflexibleprices_*
!rm bbeffectivedemandghhflexibleprices.log
!rm bbeffectivedemandghhflexibleprices.jnl
!rm bbeffectivedemandghhflexibleprices.m

dynare bbeffectivedemand.mod noclearall

!rm bbeffectivedemand_*
!rm bbeffectivedemand.log
!rm bbeffectivedemand.jnl
!rm bbeffectivedemand.m


dynare bbeffectivedemandghh.mod noclearall

!rm bbeffectivedemandghh_*
!rm bbeffectivedemandghh.log
!rm bbeffectivedemandghh.jnl
!rm bbeffectivedemandghh.m

dynare bbeffectivedemandlowerlse.mod noclearall

!rm bbeffectivedemandlowerlse_*
!rm bbeffectivedemandlowerlse.log
!rm bbeffectivedemandlowerlse.jnl
!rm bbeffectivedemandlowerlse.m

cd ../../

clear
clc
close all


clear
clc
close all

cd SSSvsGIRF

dynare bbeffectivedemand.mod noclearall

!rm bbeffectivedemand_*
!rm bbeffectivedemand.log
!rm bbeffectivedemand.jnl
!rm bbeffectivedemand.m

dynare bbeffectivedemandsimulation.mod noclearall

!rm bbeffectivedemandsimulation_*
!rm bbeffectivedemandsimulation.log
!rm bbeffectivedemandsimulation.jnl
!rm bbeffectivedemandsimulation.m









