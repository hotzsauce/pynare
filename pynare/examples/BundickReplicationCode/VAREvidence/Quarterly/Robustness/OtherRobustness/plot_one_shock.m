%% FIRST PLOT IRFS  
function [] = plot_one_shock(l_CI,IRhoriz,plot_options,shock,IRmed,IRlo_2,IRup_2,variables,shocks,vn,sn,IRlo_1,IRup_1,save_graphs,save_marker);
		
		figure;
		
		xaxis=0:IRhoriz;
		zerol=zeros(IRhoriz+1,1);
		plotcount=1;
		
		vars = plot_options(:,shock);
		nvars = size(vars,1);
		pvars = nonzeros(vars);
		Nvarplots = size(pvars,1);
		plot_vars = zeros(Nvarplots,1);
		
		for ii = 1:1:nvars;
			if vars(ii,1) == 1;
				plot_vars(1,1) = ii;
			elseif vars(ii,1) == 2;
				plot_vars(2,1) = ii;
			elseif vars(ii,1) == 3;
				plot_vars(3,1) = ii;
			elseif vars(ii,1) == 4;
				plot_vars(4,1) = ii;
			elseif vars(ii,1) == 5;
				plot_vars(5,1) = ii;
			elseif vars(ii,1) == 6;
				plot_vars(6,1) = ii;
			end;
		end;
		
		xaxis=0:IRhoriz;
		zerol=zeros(IRhoriz+1,1);
		plotcount=1;

		%% NOTE TO SELF:
		%% plot_this loops through rows - top to bottom - (which are variables)
		%% jj loops through columns - left to right - (which are shocks)
		if l_CI == 2;
			for ii = 1:1:Nvarplots;
				plot_this = plot_vars(ii);
				for jj = shock:1:shock;
					%%
					yaxmin=min([min(IRlo_2(:,plot_this,jj)), 0]);
					yaxmin=yaxmin-0.3*abs(yaxmin);
					yaxmax=max([max(IRup_2(:,plot_this,jj)), 0]);
					yaxmax=yaxmax+0.3*abs(yaxmax);
					%%
					subplot(Nvarplots,1,plotcount);
					hold on;
					jbfill(xaxis,IRup_2(:,plot_this,jj)',IRlo_2(:,plot_this,jj)',[.75, .75, .75],[1 1 1],0,1);
					plot(xaxis, zerol, '-k', 'LineWidth', 1.0);
					plot(xaxis, IRlo_2(:,plot_this,jj), '-b','LineWidth', 2.0);
					plot(xaxis, IRup_2(:,plot_this,jj), '-b','LineWidth', 2.0);
					plot(xaxis, IRmed(:,plot_this,jj), '-r', 'LineWidth', 1.5);	 
					%% THE FOLLOWING BOXES IN THE PLOTS (%% COMMENT OUT IF UNDERSIRED)
					%plot([0, 0], [yaxmin, yaxmax], 'k');
					%plot([IRhoriz, IRhoriz], [yaxmin, yaxmax], 'k');
					%plot([0, IRhoriz], [yaxmin, yaxmin], 'k');
					%plot([0, IRhoriz], [yaxmax, yaxmax], 'k');
					%%
					xlim([0 IRhoriz]);
					ylim([yaxmin yaxmax]);
					%%
					if jj == 1;
						ylabel(['\fontsize{10}\bf{% Change}']);
					end;
					%%
					if plotcount == Nvarplots;
						xlabel(['\fontsize{10}\bf{Months}']);
					end;
					axis([0 IRhoriz yaxmin yaxmax]);
					title_variable = variables.(vn{plot_this});
					title_shock    = shocks.(sn{jj});
					if ii == 1;
						title_stack	   = {['\fontsize{11}\bf{' title_shock ' Shock }'] ,'\fontsize{4} ', ['\fontsize{11}\rm{' title_variable '}']};   
					else
						title_stack	   = ['\fontsize{11}\rm{' title_variable '}'];
					end;
					title(title_stack);
					hold off;
					plotcount=plotcount+1;
				end; %% jj
			end %% ii
		elseif l_CI == 4;
			for ii = 1:1:Nvarplots;
				plot_this = plot_vars(ii);
				for jj = shock:1:shock;
					%%
					yaxmin=min([min(IRlo_1(:,plot_this,jj)), 0]);
					yaxmin=yaxmin-0.3*abs(yaxmin);
					yaxmax=max([max(IRup_1(:,plot_this,jj)), 0]);
					yaxmax=yaxmax+0.3*abs(yaxmax);
					%%
					subplot(Nvarplots,1,plotcount);
					hold on;
					jbfill(xaxis,IRup_1(:,plot_this,jj)',IRlo_1(:,plot_this,jj)',[.75, .75, .75],[1 1 1],0,1);
					plot(xaxis, zerol, '-k', 'LineWidth', 1.0);
					plot(xaxis, IRlo_2(:,plot_this,jj), '--b','LineWidth', 2.0);
					plot(xaxis, IRup_2(:,plot_this,jj), '--b','LineWidth', 2.0);
					plot(xaxis, IRlo_1(:,plot_this,jj), '--b','LineWidth',1.00);
					p2 = plot(xaxis, IRup_1(:,plot_this,jj), '--b','LineWidth',1.00);
					p1 = plot(xaxis, IRmed(:,plot_this,jj), '-r', 'LineWidth', 1.5);	 
					%% THE FOLLOWING BOXES IN THE PLOTS (%% COMMENT OUT IF UNDERSIRED)
					%plot([0, 0], [yaxmin, yaxmax], 'k');
					%plot([IRhoriz, IRhoriz], [yaxmin, yaxmax], 'k');
					%plot([0, IRhoriz], [yaxmin, yaxmin], 'k');
					%plot([0, IRhoriz], [yaxmax, yaxmax], 'k');
					%%
					xlim([0 IRhoriz]);
					ylim([yaxmin yaxmax]);
					%%
					%if jj == 1;
						ylabel(['\fontsize{10}\bf{% Change}']);
					%end;
					%%
					if plotcount == Nvarplots;
						xlabel(['\fontsize{10}\bf{Months}']);
					end;
					axis([0 IRhoriz yaxmin yaxmax]);
					title_variable = variables.(vn{plot_this});
					title_shock    = shocks.(sn{jj});
					if ii == 1;
						title_stack	   = {['\fontsize{11}\bf{' title_shock ' Shock }'] ,'\fontsize{4} ', ['\fontsize{11}\rm{' title_variable '}']};   
					else
						title_stack	   = ['\fontsize{11}\rm{' title_variable '}'];
					end;
					title(title_stack);
					hold off;
					plotcount=plotcount+1;
				end; %% jj
			end %% ii
		end; %% if/elseif CI
		%%
		if save_graphs == 1;
			title_graph = regexprep(title_shock,' ','_');
			title_graph = [save_marker,'var',title_graph,'.eps'];
			print(gcf,'-depsc',title_graph);
		end;