%% FIRST PLOT IRFS  
function [] = plot_all_shocks(l_CI,IRhoriz,IRmed,IRlo_2,IRup_2,variables,shocks,vn,sn,IRlo_1,IRup_1,save_graphs,save_marker);
	figure;
	
	Nbig = size(IRmed,2);
	Nshocks = size(IRmed,3);
	xaxis=0:IRhoriz;
    zerol=zeros(IRhoriz+1,1);
    plotcount=1;
	
	%% NOTE TO SELF:
	%% ii loops through rows - top to bottom - (which are variables)
	%% jj loops through columns - left to right - (which are shocks)
	if l_CI == 2;
		for ii = 1:1:Nbig;
			for jj = 1:1:Nshocks; %jj = 1; 
				%%
				yaxmin=min([min(IRlo_2(:,ii,jj)), 0]);
				yaxmin=yaxmin-0.3*abs(yaxmin);
				yaxmax=max([max(IRup_2(:,ii,jj)), 0]);
				yaxmax=yaxmax+0.3*abs(yaxmax);
				%%
				subplot(Nbig,Nshocks,plotcount);
				hold on;
				jbfill(xaxis,IRup_2(:,ii,jj)',IRlo_2(:,ii,jj)',[.75, .75, .75],[1 1 1],0,1);
				plot(xaxis, zerol, '-k', 'LineWidth', 1.0);
				plot(xaxis, IRlo_2(:,ii,jj), '-b','LineWidth', 2.0);
				plot(xaxis, IRup_2(:,ii,jj), '-b','LineWidth', 2.0);
				plot(xaxis, IRmed(:,ii,jj), '-r', 'LineWidth', 1.5);	 
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
				if ii == Nbig;
					xlabel(['\fontsize{10}\bf{Months}']);
				end;
				axis([0 IRhoriz yaxmin yaxmax]);
				title_variable = variables.(vn{ii});
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
		for ii = 1:1:Nbig;
			for jj = 1:1:Nshocks; %jj = 1; 
				%%
				yaxmin=min([min(IRlo_1(:,ii,jj)), 0]);
				yaxmin=yaxmin-0.3*abs(yaxmin);
				yaxmax=max([max(IRup_1(:,ii,jj)), 0]);
				yaxmax=yaxmax+0.3*abs(yaxmax);
				%%
				subplot(Nbig,Nshocks,plotcount);
				hold on;
				jbfill(xaxis,IRup_1(:,ii,jj)',IRlo_1(:,ii,jj)',[.75, .75, .75],[1 1 1],0,1);
				plot(xaxis, zerol, '-k', 'LineWidth', 1.0);
				plot(xaxis, IRlo_2(:,ii,jj), '-b','LineWidth', 2.0);
				plot(xaxis, IRup_2(:,ii,jj), '-b','LineWidth', 2.0);
				plot(xaxis, IRlo_1(:,ii,jj), '--b','LineWidth',1.00);
				plot(xaxis, IRup_1(:,ii,jj), '--b','LineWidth',1.00);
				plot(xaxis, IRmed(:,ii,jj), '-r', 'LineWidth', 1.5);	 
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
				if ii == Nbig;
					xlabel(['\fontsize{10}\bf{Months}']);
				end;
				axis([0 IRhoriz yaxmin yaxmax]);
				title_variable = variables.(vn{ii});
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
	if save_graphs == 1;
		title_graph = [save_marker,'var_all_shocks'];
		title_graph = [title_graph,'.eps'];
		print(gcf,'-deps',title_graph);
	end;