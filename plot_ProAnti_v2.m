function plot_ProAnti_v2(units,CS_modulating_cells,CS_rates, plotType, recArea)

% December 2018
% Input
%           load file: CS_modulating cells.
%           load file: units  - output from extractWholeNeuronResults.m
%           load file: CS_rates
%           plotType (see below for plot list)
%           recArea (vermis or lateral)


% Plot list
%
% 'raster_instr': instruction aligned raster plot for the chosen cell
% 'psth': psth for the chosen cell, aligned to saccade and instruction
% 'delta_change'
% 'bino_dist'

switch plotType
    case 'psth_instr_signif'
        %%
        % get area
        for cellNum = 1:length(units)
            indx(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx = find(indx);
        nunits_area = 1:length(indx);
        
        % get signif
        for i = 1:length(indx)
            indx_instr(i) = logical(units(indx(i)).stats.instr.flags.proVsAnti_instr);
        end
        instr_signif = units(indx_instr); nunits_signif = sum(indx_instr); units_signif= indx(indx_instr);
        p=numSubplots(length(instr_signif));
        
        %extract and plot
        for j=1:nunits_signif
            t= instr_signif(j).pro.neural.instr.ts_pst; % time
            r_pro= instr_signif(j).pro.neural.instr.rate_pst; % psth
            std_pro = std(instr_signif(j).pro.neural.instr.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            sem_pro = std(instr_signif(j).pro.neural.instr.rate_pst)/sqrt(length(units(indx(i)).pro.neural.trial));
            sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
            r_anti = instr_signif(j).anti.neural.instr.rate_pst;
            std_anti = std(instr_signif(j).anti.neural.instr.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= instr_signif(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro = shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti = shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([0 0.4]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            %           xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(units_signif(j))]);
            
            % plot zstat bino dist
            figure(2); subplot(p(1),p(2),j); hold on;
            stat_instr(i,:) = instr_signif(j).stats.instr.pval.pbDist_testStat;
            t_instr = instr_signif(j).pro.neural.instr.ts_pst;
            plot(t_instr, stat_instr(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[0 0.4],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            %           ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(units_signif(j))]);
            
            %             waitforbuttonpress; close all;
        end
        
    case 'psth_instr_CS_mod'
        %%
        % get cells
        indx = CS_modulating_cells.(recArea).either.instruction;
        nunits_area = 1:length(indx);
        % get signif
        for i = 1:length(indx)
            indx_signif(i) = logical(units(indx(i)).stats.instr.flags.proVsAnti_instr);
        end
        cs_units = units(indx);
        space_width = 0.15;
        p=numSubplots(length(indx));
        
        %extract
        for j=1:length(cs_units)
            t= cs_units(j).pro.neural.instr.ts_pst; % time
            r_pro= cs_units(j).pro.neural.instr.rate_pst; % psth
            std_pro = std(cs_units(j).pro.neural.instr.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = cs_units(j).anti.neural.instr.rate_pst;
            std_anti = std(cs_units(j).anti.neural.instr.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= cs_units(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([0 0.4]), 'TickDir', 'out', 'FontSize',12); % analysis window size
            %xlabel('Time (s)'); ylabel ('FR (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            %         annotation('textbox',...
            %             [0.148482014388489 0.147388059701493 0.0829328537170263 0.0671641791044775],...
            %             'String',['ProVSAnti = ' num2str(indx_signif(j))],...
            %             'FitBoxToText','on');
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            %             if j ==1
            %                 xlabel('Time (s)'); ylabel ('FR (spk/s)');
            %             end
            
            
            % plot zstat bino dist
            figure(2); subplot(p(1),p(2),j); hold on;
            stat_instr(j,:) = cs_units(j).stats.instr.pval.pbDist_testStat;
            t_instr = cs_units(j).pro.neural.instr.ts_pst;
            position_stat_sign_instr(j,:) = abs(stat_instr(j,:))>=1.96;
            %figure('Position',[646 836 560 172]);
            plot(t_instr, stat_instr(j,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[0 0.4],'ylim',[-4 4], 'TickDir', 'out', 'FontSize', 18); %ylabel('z-stat');
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]); box off;
            %             if j ==1
            %                 xlabel('Time (s)'); ylabel ('z-stat');
            %             end
            
            %       waitforbuttonpress; close all;
        end
        
        shared_cells = indx(indx_signif);
        p=numSubplots(length(shared_cells));
        % SS
        for i = 1:length(shared_cells)
            r_pro = units(shared_cells(i)).pro.neural.instr.rate_pst; % psth
            std_pro = std(units(shared_cells(i)).pro.neural.instr.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = units(shared_cells(i)).anti.neural.instr.rate_pst; % psth
            std_anti = std(units(shared_cells(i)).anti.neural.instr.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            
            figure(1);subplot(p(1),p(2),i); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([0 0.4]), 'TickDir', 'out', 'FontSize',12); % analysis window size
            set(s_pro.mainLine,'LineWidth', 3), set(s_anti.mainLine,'LineWidth', 3);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(shared_cells(i))]);
        end
        % CS
        for i = 1:length(shared_cells)
            t = CS_rates.lateral.timepoints_instruction;
            r_pro = CS_rates.lateral.rate_pro_instr(shared_cells(i),:); % psth
            std_pro = std(r_pro);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.lateral.rate_anti_instr(shared_cells(i),:); % psth
            std_anti = std(r_anti);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            
            figure(2);subplot(p(1),p(2),i); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([0 0.4]), 'TickDir', 'out', 'FontSize',12); % analysis window size
            set(s_pro.mainLine,'LineWidth', 3), set(s_anti.mainLine,'LineWidth', 3);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(shared_cells(i))]);
        end
        
        
    case 'psth_sacc_CS_mod'
        %%
        % get cells
        indx = CS_modulating_cells.(recArea).either.sacc;
        nunits_area = 1:length(indx);
        % get signif
        for i = 1:length(indx)
            indx_signif(i) = logical(units(indx(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        cs_units = units(indx);
        p=numSubplots(length(indx));
        
        %extract
        for j=1:length(indx)
            t= cs_units(j).pro.neural.sacc.ts_pst; % time
            r_pro= cs_units(j).pro.neural.sacc.rate_pst; % psth
            std_pro = std(cs_units(j).pro.neural.sacc.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = cs_units(j).anti.neural.sacc.rate_pst;
            std_anti = std(cs_units(j).anti.neural.sacc.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= cs_units(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            %             annotation('textbox',...
            %                 [0.148482014388489 0.147388059701493 0.0829328537170263 0.0671641791044775],...
            %                 'String',['ProVSAnti = ' num2str(indx_signif(j))],...
            %                 'FitBoxToText','on');
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(2); subplot(p(1),p(2),j); hold on;
            stat_sacc(i,:) = cs_units(j).stats.sacc.pval.pbDist_testStat;
            t_sacc = cs_units(j).pro.neural.sacc.ts_pst;
            plot(t_sacc, stat_sacc(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.15 0.15],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            %             waitforbuttonpress; close all;
        end
        
        shared_cells = indx(indx_signif);
        p=numSubplots(length(shared_cells));
        % SS
        for i = 1:length(shared_cells)
            t = units(shared_cells(i)).pro.neural.sacc.ts_pst;
            r_pro = units(shared_cells(i)).pro.neural.sacc.rate_pst; % psth
            std_pro = std(units(shared_cells(i)).pro.neural.sacc.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = units(shared_cells(i)).anti.neural.sacc.rate_pst; % psth
            std_anti = std(units(shared_cells(i)).anti.neural.sacc.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            
            figure(1);subplot(p(1),p(2),i); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([-0.15 0.151]), 'TickDir', 'out', 'FontSize',12); % analysis window size
            set(s_pro.mainLine,'LineWidth', 3), set(s_anti.mainLine,'LineWidth', 3);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(shared_cells(i))]);
        end
        % CS
        for i = 1:length(shared_cells)
            t = CS_rates.lateral.timepoints_saccade;
            r_pro = CS_rates.lateral.rate_pro_sacc(shared_cells(i),:); % psth
            std_pro = std(r_pro);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.lateral.rate_anti_sacc(shared_cells(i),:); % psth
            std_anti = std(r_anti);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            
            figure(2);subplot(p(1),p(2),i); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([-0.15 0.151]), 'TickDir', 'out', 'FontSize',12); % analysis window size
            set(s_pro.mainLine,'LineWidth', 3), set(s_anti.mainLine,'LineWidth', 3);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(shared_cells(i))]);
        end
        
        
    case 'psth_sacc_cs_ss_vermis'
        % get indx from CS
        indx_cs = CS_rates.vermis.indx_area;
        ss_units = units(indx_cs);
        p=numSubplots(length(indx_cs));
        for i = 1:length(indx_cs)
            indx_signif(i) = logical(ss_units(i).stats.sacc.flags.proVsAnti_sacc);
        end
        
        % CS
        for j=1:length(indx_cs)
            t_cs= CS_rates.vermis.timepoints_saccade; % time
            r_pro= CS_rates.vermis.rate_pro_sacc(indx_cs(j),:); % psth
            std_pro = CS_rates.vermis.std_pro(indx_cs(j),:);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.vermis.rate_anti_sacc(indx_cs(j),:);
            std_anti = CS_rates.vermis.std_anti(indx_cs(j),:);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t_cs,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t_cs,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([-0.15 0.151]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx_cs(j)) ' sig ' num2str(indx_cs(j))]);
        end
        
        for j=1:length(indx_cs)
            t= ss_units(j).pro.neural.sacc.ts_pst; % time
            r_pro= ss_units(j).pro.neural.sacc.rate_pst; % psth
            std_pro = std(ss_units(j).pro.neural.sacc.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = ss_units(j).anti.neural.sacc.rate_pst;
            std_anti = std(ss_units(j).anti.neural.sacc.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= ss_units(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(2);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(3); subplot(p(1),p(2),j); hold on;
            stat_sacc(i,:) = cs_units(j).stats.sacc.pval.pbDist_testStat;
            t_sacc = cs_units(j).pro.neural.sacc.ts_pst;
            position_stat_sign_instr(i,:) = abs(stat_sacc(i,:))>=1.96;
            plot(t_sacc, stat_sacc(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.15 0.15],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
        end
        
    case 'psth_sacc_cs_ss_lateral'
        % get indx from CS
        indx_cs = CS_rates.lateral.indx_area;
        ss_units = units(indx_cs);
        p=numSubplots(length(indx_cs));
        for i = 1:length(indx_cs)
            indx_signif(i) = logical(ss_units.stats.sacc.flags.proVsAnti_sacc);
        end
        
        % CS
        for j=1:length(indx_cs)
            t_cs= CS_rates.lateral.timepoints_saccade; % time
            r_pro= CS_rates.lateral.rate_pro_sacc(indx_cs(j),:); % psth
            std_pro = CS_rates.lateral.std_pro(indx_cs(j),:);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.lateral.rate_anti_sacc(indx_cs(j),:);
            std_anti = CS_rates.lateral.std_anti(indx_cs(j),:);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t_cs,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t_cs,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx_cs(j)) ' sig ' num2str(indx_cs(j))]);
        end
        
        for j=1:length(indx_cs)
            t= ss_units(j).pro.neural.sacc.ts_pst; % time
            r_pro= ss_units(j).pro.neural.sacc.rate_pst; % psth
            std_pro = std(ss_units(j).pro.neural.sacc.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = ss_units(j).anti.neural.sacc.rate_pst;
            std_anti = std(ss_units(j).anti.neural.sacc.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= ss_units(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(2);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(3); subplot(p(1),p(2),j); hold on;
            stat_sacc(i,:) = cs_units(j).stats.sacc.pval.pbDist_testStat;
            t_sacc = cs_units(j).pro.neural.sacc.ts_pst;
            plot(t_sacc, stat_sacc(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.15 0.15],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
        end
    case 'psth_instr_cs_ss_vermis'
        % get indx from CS
        indx_cs = CS_rates.vermis.indx_area;
        ss_units = units(indx_cs);
        p=numSubplots(length(indx_cs));
        for i = 1:length(indx_cs)
            indx_signif(i) = logical(ss_units(i).stats.instr.flags.proVsAnti_sacc);
        end
        
        % CS
        for j=1:length(indx_cs)
            t_cs= CS_rates.vermis.timepoints_instruction; % time
            r_pro= CS_rates.vermis.rate_pro_instr(indx_cs(j),:); % psth
            std_pro = std(r_pro);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.vermis.rate_anti_instr(indx_cs(j),:);
            std_anti = std(r_anti);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t_cs,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t_cs,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx_cs(j)) ' sig ' num2str(indx_cs(j))]);
        end
        
        for j=1:length(indx_cs)
            t= ss_units(j).pro.neural.instr.ts_pst; % time
            r_pro= ss_units(j).pro.neural.instr.rate_pst; % psth
            std_pro = std(ss_units(j).pro.neural.instr.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = ss_units(j).anti.neural.instr.rate_pst;
            std_anti = std(ss_units(j).anti.neural.instr.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= ss_units(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(2);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(3); subplot(p(1),p(2),j); hold on;
            stat_instr(i,:) = cs_units(j).stats.instr.pval.pbDist_testStat;
            t_instr = cs_units(j).pro.neural.instr.ts_pst;
            plot(t_instr, stat_instr(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.15 0.15],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
        end
        
    case 'psth_instr_cs_ss_lateral'
        % get indx from CS
        indx_cs = CS_rates.lateral.indx_area;
        ss_units = units(indx_cs);
        p=numSubplots(length(indx_cs));
        for i = 1:length(indx_cs)
            indx_signif(i) = logical(ss_units(i).stats.sacc.flags.proVsAnti_sacc);
        end
        
        % CS
        for j=1:length(indx_cs)
            t_cs= CS_rates.lateral.timepoints_instruction; % time
            r_pro= CS_rates.lateral.rate_pro_instr(indx_cs(j),:); % psth
            std_pro = CS_rates.lateral.std_pro(indx_cs(j),:);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.lateral.rate_anti_instr(indx_cs(j),:);
            std_anti = CS_rates.lateral.std_anti(indx_cs(j),:);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t_cs,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t_cs,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx_cs(j)) ' sig ' num2str(indx_cs(j))]);
        end
        
        for j=1:length(indx_cs)
            t= ss_units(j).pro.neural.instr.ts_pst; % time
            r_pro= ss_units(j).pro.neural.instr.rate_pst; % psth
            std_pro = std(ss_units(j).pro.neural.instr.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = ss_units(j).anti.neural.instr.rate_pst;
            std_anti = std(ss_units(j).anti.neural.instr.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= ss_units(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(2);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(3); subplot(p(1),p(2),j); hold on;
            stat_instr(i,:) = cs_units(j).stats.instr.pval.pbDist_testStat;
            t_instr = cs_units(j).pro.neural.instr.ts_pst;
            plot(t_instr, stat_instr(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.15 0.15],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
        end
        
    case 'delta_rate'
        %%
        % get cells
        indx = CS_modulating_cells.(recArea).either.sacc;
        nunits_area = 1:length(indx);
        % get signif
        for i = 1:length(indx)
            indx_signif(i) = logical(units(indx(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        cs_units = units(indx);
        
        for j=1:length(indx)
            delta_pro_base(j,:)=cs_units(j).pro.neural.sacc.delta_rate_base;
            delta_anti_base(j,:)=cs_units(j).anti.neural.sacc.delta_rate_base;
            max_delta_pro_base(j) = max(abs(cs_units(j).pro.neural.sacc.delta_rate_base));
            max_delta_anti_base(j) = max(abs(cs_units(j).anti.neural.sacc.delta_rate_base));
        end
        t = cs_units(1).pro.neural.sacc.ts_pst_win;
        mean_delta_pro = mean(delta_pro_base,1);
        sem_delta_pro = std(delta_pro_base,0,1)/sqrt(numel(indx));
        mean_delta_anti = mean(delta_anti_base,1);
        sem_delta_anti = std(delta_anti_base,0,1)/sqrt(numel(indx));
        
        figure; hold on;
        s_pro=shadedErrorBar(t, mean_delta_pro, sem_delta_pro, 'lineprops','r');
        s_anti=shadedErrorBar(t, mean_delta_anti, sem_delta_anti, 'lineprops','g');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca,'TickDir', 'out', 'FontSize', 18,'ylim',[-2 10],'ytick',[-2 0 10]);
        xlabel('Time (s)'); ylabel('Change in firing rate (Hz)')
        title('CS mod cells')
        
        %stat over time
        [h_change,p_change] = ttest(abs(delta_anti_base), abs(delta_pro_base));
        h_change(h_change==0)=NaN;
        % diff anti-pro change in FR sacc
        figure; hold on;
        plot(t,nanmean(abs(delta_anti_base)-abs(delta_pro_base)),'Color','k', 'LineWidth', 2);
        plot(t,h_change*0.5, '*c')
        set(gca,'TickDir','out','ylim',[-2 2], 'ytick',[-2 0 2], 'FontSize', 18)
        xlabel('Time (s)'); ylabel('Abs change in FR anti-pro')
        title('CS mod cells')
        
        
        % plot max change in FR
        % plot scatter for significantly diff cells from mean
        figure; hold on;
        plot(max_delta_pro_base,max_delta_anti_base, '.k','MarkerSize', 18);
        %set(gca,'XScale','Log','YScale','Log' ,'FontSize', 18, 'TickDir', 'out');axis ([1e0 1e2 1e0 1e2]);
        set(gca,'FontSize', 18, 'TickDir', 'out', 'ylim',[0 30],'xlim',[0 30],'yTick', [0 30],'xTick', [0 30]); %axis ([1e0 1e2 1e0 1e2]);
        % plot([1e0 1e2],[1e0 1e2]);
        plot([0 30], [0 30]);
        xlabel('Max change pro'); ylabel('Max change anti');
        title(['Max change in firing rate >> ' recArea])
        axis square
        [h,p] = ttest(max_delta_pro_base,max_delta_anti_base)
        
    case 'bino_dist_CS_mod_sacc'
        % get cells
        indx = CS_modulating_cells.(recArea).either.sacc;
        nunits_area = 1:length(indx);
        % get signif
        for i = 1:length(indx)
            indx_signif(i) = logical(units(indx(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        cs_units = units(indx);
        
        for j=1:length(indx)
            stat_sacc(j,:) = cs_units(j).stats.sacc.pval.pbDist_testStat;
            position_stat_sign_sacc(j,:) = abs(stat_sacc(j,:))>=1.96;
            t= cs_units(j).pro.neural.instr.ts_pst; % time
            r_pro(j,:)= cs_units(j).pro.neural.sacc.rate_pst; % psth
            r_anti(j,:)= cs_units(j).anti.neural.sacc.rate_pst; % psth
        end
        t_sacc = cs_units(1).pro.neural.sacc.ts_pst;
        
        % plot
        figure; hold on;
        plot(t_sacc,abs(stat_sacc), '.', 'MarkerSize', 18);
        set(gca, 'xlim',[-0.15 0.151], 'TickDir', 'out', 'FontSize', 18);
        hline(1.96); vline(0);xlabel('time(s)'); ylabel('Z-stat')
        %histogram
        figure;histogram(abs(stat_sacc),25);
        set(gca, 'XDir','reverse','xlim',[0 5],'TickDir', 'out', 'FontSize', 18);
        vline(1.96); box off;
        
        figure; hold on;
        plot(t_sacc,smooth(nanmean(abs(stat_sacc)),3),'LineWidth', 2,'Color','k'); % smoothed
        hline(1.96,'k')
        set(gca, 'xlim',[-0.15 0.151],'ylim',[0 2], 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Sacc => ' recArea ' n= ' num2str(numel(indx))]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        % colormap
        abs_stat_sacc = abs(stat_sacc);
        t_sacc_win = t_sacc(t_sacc>-0.15 & t_sacc<0.15);
        stat_sacc_win = abs_stat_sacc(:,t_sacc>-0.15 & t_sacc<0.15);
        figure; colormap(bone);
        imagesc(t_sacc_win,1:size(stat_sacc_win,1),stat_sacc_win, [0 1.96]);
        set(gca,'xlim',[-0.15 0.151],'YTickLabel',[], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off;
        signif = false(size(stat_sacc_win,1),1); count=1;
        for i=1:size(stat_sacc_win,1)
            if sum(stat_sacc_win(i,:)>1.96)>=2
                is = find(stat_sacc_win(i,:) > 1.96);
                indx_stat(count) = is(1); count=count+1;
                signif(i) = 1;
            end
        end
        signif = find(signif);
        figure('Position',[1160 678 80 420]); plot(1,signif, '*k'); box off;
        set(gca, 'xlim',[0.9 1.1], 'ylim',[0 25], 'YDir', 'reverse', 'xTick', []);
        %
        %         % sort
        %         [~,sorted_indx] = sort(indx_stat);
        %         sorted_stat = stat_sacc_win(sorted_indx,:);
        
        % sort by max
        [y,I] = max(stat_sacc_win,[],2);
        
        [~,sorted_indx] = sort(I);
        sorted_stat = stat_sacc_win(sorted_indx,:);
        
        % plot
        figure; colormap(bone); imagesc(t_sacc_win,1:size(sorted_stat,1),sorted_stat,[0 1.96]);
        set(gca,'xlim',[-0.15 0.151], 'YTickLabel',[], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; axis off;
        
        % waterfall plot
        figure; hold on; colormap(winter);
        for j=1:size(sorted_stat,1)
            waterfall(t_sacc_win,j,sorted_stat(j,:));
        end
        set(gca,'xlim', [-0.15 0.151],'zlim', [0 6], 'CameraViewAngle', 9, 'zTick', [], 'ytick', []);
        view(gca,[-3.2000 86.8000]);
        grid off;
        title(['Sorted Z-stat'])
        
        %% plot FR according to this order and normalize FR
        for i=1:size(r_pro,1)
            r_pro_win(i,:) = r_pro(i, t_sacc>0.05 & t_sacc<0.350);
            r_anti_win(i,:) = r_anti(i, t_sacc>0.05 & t_sacc<0.350);
            r_subtracted(i,:) = r_anti_win(i,:) - r_pro_win(i,:); % subtract anti-pro for each cell and plot.
        end
        % sort based on index above
        r_pro_sorted = r_pro_win(sorted_indx,:);
        r_anti_sorted = r_anti_win(sorted_indx,:);
        r_subtracted_sorted = r_subtracted(sorted_indx,:);
        
        % plot FR subtracted
        cmap = colormap(bone);
        colormap(cmap);
        imagesc(t_sacc_win,1:size(r_subtracted_sorted,1),r_subtracted_sorted, [-15 15]);
        set(gca,'xlim',[-0.15 0.151],'YTickLabel',[],'xTick',[-0.15 0.15],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title('Anti-Pro')
        
        % normalize
        [maxR_pro,~]=max(r_pro_sorted,[],2); [maxR_anti,~]=max(r_anti_sorted,[],2);
        for i=1:size(r_pro_sorted,1)
            r_pro_norm(i,:) = r_pro_sorted(i,:)./repmat(maxR_pro(i),[1 size(r_pro_sorted(i,:),2)]);
            r_anti_norm(i,:) = r_anti_sorted(i,:)./repmat(maxR_anti(i),[1 size(r_anti_sorted(i,:),2)]);
        end
        
        % firing rate
        cmap = colormap(bone); % cmap = flipud(cmap);
        colormap(cmap);
        imagesc(t_sacc_win,1:size(r_pro_norm,1),r_pro_norm, [0 1]);
        set(gca,'xlim',[-0.15 0.151],'YTickLabel',[],'xTick',[-0.15 0.15],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title('Pro');
        
        figure; colormap(cmap);
        imagesc(t_sacc_win,1:size(r_anti_norm,1),r_anti_norm, [0 1]);
        set(gca,'xlim',[-0.15 0.151],'YTickLabel',[],'xTick',[-0.15 0.15],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title('Anti')
        
        figure; colormap(cmap);
        imagesc(t_sacc_win,1:size(r_anti_norm,1),r_anti_norm-r_pro_norm, [-0.1 0.1]);
        set(gca,'xlim',[-0.15 0.151],'YTickLabel',[],'xTick',[-0.15 0.15],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title('Anti-Pro')
        
        
    case 'bino_dist_CS_mod_instr'
        % get cells
        % get cells
        indx = CS_modulating_cells.(recArea).either.instruction;
        nunits_area = 1:length(indx);
        % get signif
        for i = 1:length(indx)
            indx_signif(i) = logical(units(indx(i)).stats.instr.flags.proVsAnti_instr);
        end
        cs_units = units(indx);
        
        for j=1:length(indx)
            stat_instr(j,:) = cs_units(j).stats.instr.pval.pbDist_testStat;
            position_stat_sign_instr(j,:) = abs(stat_instr(j,:))>=1.96;
            t= cs_units(j).pro.neural.instr.ts_pst; % time
            r_pro(j,:)= cs_units(j).pro.neural.instr.rate_pst; % psth
            r_anti(j,:)= cs_units(j).anti.neural.instr.rate_pst; % psth
        end
        t_instr = cs_units(1).pro.neural.instr.ts_pst;
        
        % plot
        figure; hold on;
        plot(t_instr,abs(stat_instr), '.', 'MarkerSize', 18);
        set(gca, 'xlim',[0.05 0.350], 'TickDir', 'out', 'FontSize', 18);
        hline(1.96); vline(0);xlabel('time(s)'); ylabel('Z-stat')
        %histogram
        figure;histogram(abs(stat_instr),25);
        set(gca, 'XDir','reverse','xlim',[0 6],'TickDir', 'out', 'FontSize', 18);
        vline(1.96); box off;
        
        figure; hold on;
        plot(t_instr,smooth(nanmean(abs(stat_instr)),3),'LineWidth', 2,'Color','k'); % smoothed
        hline(1.96,'k')
        set(gca, 'xlim',[0.05 0.350],'ylim',[0 2], 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Sacc => ' recArea ' n= ' num2str(numel(indx))]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        % colormap
        abs_stat_instr = abs(stat_instr);
        t_instr_win = t_instr(t_instr>0.05 & t_instr<0.350);
        stat_instr_win = abs_stat_instr(:,t_instr>0.05 & t_instr<0.350);
        figure; colormap(bone);
        imagesc(t_instr_win,1:size(stat_instr_win,1),stat_instr_win, [0 1.96]);
        set(gca,'xlim',[0.05 0.350],'YTickLabel',[], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off;
        signif = false(size(stat_instr_win,1),1); count=1;
        %         for i=1:size(stat_instr_win,1)
        %             if  sum(stat_instr_win(i,:)>1.96)>=2  % ~isempty(find(stat_instr_win(i,:)>1.96))
        %                 is = find(stat_instr_win(i,:)>1.96);
        %                 indx_stat(count) = is(1); count=count+1;
        %                 signif(i) = 1;
        %             end
        %         end
        %         signif = find(signif);
        %         figure('Position',[1160 678 80 420]); plot(1,signif, '*k'); box off;
        %         set(gca, 'xlim',[0.9 1.1], 'ylim',[0 25], 'YDir', 'reverse', 'xTick', []);
        
        %         % sort
        %         [~,sorted_indx] = sort(indx_stat);
        %         sorted_stat = stat_instr_win(sorted_indx,:);
        % sort by max
        [y,I] = max(stat_instr_win,[],2);
        
        [~,sorted_indx] = sort(I);
        sorted_stat = stat_instr_win(sorted_indx,:);
        
        % plot
        figure; colormap(bone); imagesc(t_instr_win,1:size(sorted_stat,1),sorted_stat,[0 1.96]);
        set(gca,'xlim',[0.05 0.350], 'YTickLabel',[], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; axis off;
        
        % waterfall plot
        figure; hold on; colormap(winter);
        for j=1:size(sorted_stat,1)
            waterfall(t_instr_win,j,sorted_stat(j,:));
        end
        set(gca,'xlim', [0.05 0.350],'zlim', [0 6], 'CameraViewAngle', 9, 'zTick', [], 'ytick', []);
        view(gca,[-3.2000 86.8000]);
        grid off;
        title(['Sorted Z-stat'])
        
        %% plot FR according to this order and normalize FR
        for i=1:size(r_pro,1)
            r_pro_win(i,:) = r_pro(i, t_instr>0.05 & t_instr<0.350);
            r_anti_win(i,:) = r_anti(i, t_instr>0.05 & t_instr<0.350);
            r_subtracted(i,:) = r_anti_win(i,:) - r_pro_win(i,:); % subtract anti-pro for each cell and plot.
        end
        % sort based on index above
        r_pro_sorted = r_pro_win(sorted_indx,:);
        r_anti_sorted = r_anti_win(sorted_indx,:);
        r_subtracted_sorted = r_subtracted(sorted_indx,:);
        
        % plot FR subtracted
        cmap = colormap(bone);
        colormap(cmap);
        imagesc(t_instr_win,1:size(r_subtracted_sorted,1),r_subtracted_sorted, [-10 10]);
        set(gca,'xlim',[0.05 0.350],'YTickLabel',[],'xTick',[0.05 0.350],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title('Anti-Pro')
        
        % normalize
        [maxR_pro,~]=max(r_pro_sorted,[],2); [maxR_anti,~]=max(r_anti_sorted,[],2);
        for i=1:size(r_pro_sorted,1)
            r_pro_norm(i,:) = r_pro_sorted(i,:)./repmat(maxR_pro(i),[1 size(r_pro_sorted(i,:),2)]);
            r_anti_norm(i,:) = r_anti_sorted(i,:)./repmat(maxR_anti(i),[1 size(r_anti_sorted(i,:),2)]);
        end
        
        % firing rate
        cmap = colormap(bone); % cmap = flipud(cmap);
        colormap(cmap);
        imagesc(t_instr_win,1:size(r_pro_norm,1),r_pro_norm, [0 1]);
        set(gca,'xlim',[0.05 0.350],'YTickLabel',[],'xTick',[0.05 0.35],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title('Pro');
        
        figure; colormap(cmap);
        imagesc(t_instr_win,1:size(r_anti_norm,1),r_anti_norm, [0 1]);
        set(gca,'xlim',[0.05 0.350],'YTickLabel',[],'xTick',[0.05 0.350],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title('Anti')
        
        figure; colormap(cmap);
        imagesc(t_instr_win,1:size(r_anti_norm,1),r_anti_norm-r_pro_norm, [-0.1 0.1]);
        set(gca,'xlim',[0.05 0.350],'YTickLabel',[],'xTick',[0.05 0.350],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title('Anti-Pro')
        
    case 'cs_mod_all_dir_instr'
        indx=CS_modulating_cells.(recArea).all_directions.either.instruction;
        ss = units(indx);
        p=numSubplots(length(indx));
        for i = 1:length(indx)
            indx_signif(i) = logical(ss(i).stats.instr.flags.proVsAnti_instr);
        end
        
        % CS
        for j=1:length(indx)
            t_cs = CS_rates.(recArea).timepoints_instruction(1,:); % time
            r_pro= CS_rates.(recArea).rate_pro_instr(indx(j),:); % psth
            std_pro = CS_rates.(recArea).std_pro(indx(j),:);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.(recArea).rate_anti_instr(indx(j),:);
            std_anti = CS_rates.(recArea).std_anti(indx(j),:);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t_cs,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t_cs,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([0 0.4]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j))]);
        end
        
        for j=1:length(indx)
            t= ss(j).pro.neural.instr.ts_pst; % time
            r_pro= ss(j).pro.neural.instr.rate_pst; % psth
            std_pro = std(ss(j).pro.neural.instr.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = ss(j).anti.neural.instr.rate_pst;
            std_anti = std(ss(j).anti.neural.instr.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= ss(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(2);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([0 0.4]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(3); subplot(p(1),p(2),j); hold on;
            stat_instr(i,:) = ss(j).stats.instr.pval.pbDist_testStat;
            t_instr = ss(j).pro.neural.instr.ts_pst;
            plot(t_instr, stat_instr(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[0 0.4],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
        end
        
    case 'cs_mod_all_dir_sacc'
        indx=CS_modulating_cells.(recArea).all_directions.either.sacc;
        ss = units(indx);
        p=numSubplots(length(indx));
        for i = 1:length(indx)
            indx_signif(i) = logical(ss(i).stats.sacc.flags.proVsAnti_sacc);
        end
        
        % CS
        for j=1:length(indx)
            t_cs = CS_rates.(recArea).timepoints_saccade(1,:); % time
            r_pro= CS_rates.(recArea).rate_pro_sacc(indx(j),:); % psth
            std_pro = CS_rates.(recArea).std_pro(indx(j),:);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.(recArea).rate_anti_sacc(indx(j),:);
            std_anti = CS_rates.(recArea).std_anti(indx(j),:);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t_cs,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t_cs,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j))]);
        end
        
        for j=1:length(indx)
            t= ss(j).pro.neural.sacc.ts_pst; % time
            r_pro= ss(j).pro.neural.sacc.rate_pst; % psth
            std_pro = std(ss(j).pro.neural.sacc.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = ss(j).anti.neural.sacc.rate_pst;
            std_anti = std(ss(j).anti.neural.sacc.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= ss(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(2);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(3); subplot(p(1),p(2),j); hold on;
            stat_sacc(i,:) = ss(j).stats.sacc.pval.pbDist_testStat;
            t_sacc = ss(j).pro.neural.sacc.ts_pst;
            plot(t_sacc, stat_sacc(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.15 0.15],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
        end
        
        %%
        case 'cs_mod_tuned_instr'
        indx=CS_modulating_cells.(recArea).tuning.either.instruction;
        ss = units(indx);
        p=numSubplots(length(indx));
        for i = 1:length(indx)
            indx_signif(i) = logical(ss(i).stats.instr.flags.proVsAnti_instr);
        end
        
        % CS
        for j=1:length(indx)
            t_cs = CS_rates.(recArea).timepoints_instruction(1,:); % time
            r_pro= CS_rates.(recArea).rate_pro_instr(indx(j),:); % psth
            std_pro = CS_rates.(recArea).std_pro(indx(j),:);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.(recArea).rate_anti_instr(indx(j),:);
            std_anti = CS_rates.(recArea).std_anti(indx(j),:);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t_cs,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t_cs,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([0 0.4]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j))]);
        end
        
        for j=1:length(indx)
            t= ss(j).pro.neural.instr.ts_pst; % time
            r_pro= ss(j).pro.neural.instr.rate_pst; % psth
            std_pro = std(ss(j).pro.neural.instr.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = ss(j).anti.neural.instr.rate_pst;
            std_anti = std(ss(j).anti.neural.instr.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= ss(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(2);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([0 0.4]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(3); subplot(p(1),p(2),j); hold on;
            stat_instr(i,:) = ss(j).stats.instr.pval.pbDist_testStat;
            t_instr = ss(j).pro.neural.instr.ts_pst;
            plot(t_instr, stat_instr(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[0 0.4],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
        end
        
    case 'cs_mod_tuned_sacc'
        indx=CS_modulating_cells.(recArea).all_directions.either.sacc;
        ss = units(indx);
        p=numSubplots(length(indx));
        for i = 1:length(indx)
            indx_signif(i) = logical(ss(i).stats.sacc.flags.proVsAnti_sacc);
        end
        
        % CS
        for j=1:length(indx)
            t_cs = CS_rates.(recArea).timepoints_saccade(1,:); % time
            r_pro= CS_rates.(recArea).rate_pro_sacc(indx(j),:); % psth
            std_pro = CS_rates.(recArea).std_pro(indx(j),:);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = CS_rates.(recArea).rate_anti_sacc(indx(j),:);
            std_anti = CS_rates.(recArea).std_anti(indx(j),:);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t_cs,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t_cs,r_anti,std_anti,'lineprops','g');
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j))]);
        end
        
        for j=1:length(indx)
            t= ss(j).pro.neural.sacc.ts_pst; % time
            r_pro= ss(j).pro.neural.sacc.rate_pst; % psth
            std_pro = std(ss(j).pro.neural.sacc.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            r_anti = ss(j).anti.neural.sacc.rate_pst;
            std_anti = std(ss(j).anti.neural.sacc.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= ss(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(2);subplot(p(1),p(2),j); hold on;
            s_pro=shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            s_anti=shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(3); subplot(p(1),p(2),j); hold on;
            stat_sacc(i,:) = ss(j).stats.sacc.pval.pbDist_testStat;
            t_sacc = ss(j).pro.neural.sacc.ts_pst;
            plot(t_sacc, stat_sacc(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.15 0.15],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
        end
        
    case 'sorted_colormap_sacc'
        % get area
        for cellNum = 1:length(units)
            indx(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx = find(indx);
        nunits_area = 1:length(indx);
        
        t = units(1).pro.neural.sacc.ts_pst_win; % time
        for j=1:length(indx)
            r_pro(j,:) = units(indx(j)).pro.neural.sacc.rate_pst_win; % psth
            r_anti(j,:) = units(indx(j)).anti.neural.sacc.rate_pst_win; % psth
        end
        % pro
        [maxRates,pos_max] = max(r_pro, [], 2);  
        [~,indx_max] = sort(pos_max);
        r_pro_norm = r_pro./repmat(maxRates,[1 size(r_pro,2)]); 
        r_pro_sorted = r_pro_norm(indx_max,:); 
        B = goodcolormap('wr');
        
        % colormap sorted
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flipud(B'));
        imagesc(t,1:size(r_pro_sorted,1),r_pro_sorted, [0 1]);
        set(gca,'xlim',[-0.15 0.151], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted Pro ' recArea])
        % colormap unsorted
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flipud(B'));
        imagesc(t,1:size(r_pro_norm,1),r_pro_norm, [0 1]);
        set(gca,'xlim',[-0.15 0.151],'YTickLabel', [],'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Pro ' recArea])
        
        % anti
        [maxRates_anti,pos_max_anti] = max(r_anti, [], 2);  
        [~,indx_max_anti] = sort(pos_max_anti);
        r_anti_norm = r_anti./repmat(maxRates_anti,[1 size(r_anti,2)]); 
        r_anti_sorted = r_anti_norm(indx_max_anti,:); 
        B = goodcolormap('wr');
        % colormap sorted
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flipud(B'));
        imagesc(t,1:size(r_anti_sorted,1),r_anti_sorted, [0 1]);
        set(gca,'xlim',[-0.15 0.151], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc sorted Anti ' recArea])
        % colormap unsorted
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flipud(B'));
        imagesc(t,1:size(r_anti_norm,1),r_anti_norm, [0 1]);
        set(gca,'xlim',[-0.15 0.151],'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Anti ' recArea])
        
    case 'sorted_colormap_instr'
        % get area
        for cellNum = 1:length(units)
            indx(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx = find(indx);
        nunits_area = 1:length(indx);
        
        t = units(1).pro.neural.instr.ts_pst_win; % time
        for j=1:length(indx)
            r_pro(j,:) = units(indx(j)).pro.neural.instr.rate_pst_win; % psth
            r_anti(j,:) = units(indx(j)).anti.neural.instr.rate_pst_win; % psth
        end
        % pro
        [maxRates,pos_max] = max(r_pro, [], 2);  
        [~,indx_max] = sort(pos_max);
        r_pro_norm = r_pro./repmat(maxRates,[1 size(r_pro,2)]); 
        r_pro_sorted = r_pro_norm(indx_max,:); 
        B = goodcolormap('wr');
        
        % colormap sorted
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flipud(B'));
        imagesc(t,1:size(r_pro_sorted,1),r_pro_sorted, [0 1]);
        set(gca,'xlim',[0.05 0.350], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title(['Instr Sorted Pro ' recArea])
        % colormap unsorted
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flipud(B'));
        imagesc(t,1:size(r_pro_norm,1),r_pro_norm, [0 1]);
        set(gca,'xlim',[0.05 0.350],'YTickLabel', [],'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off; title(['Instr Pro ' recArea])
        
        % anti
        [maxRates_anti,pos_max_anti] = max(r_anti, [], 2);  
        [~,indx_max_anti] = sort(pos_max_anti);
        r_anti_norm = r_anti./repmat(maxRates_anti,[1 size(r_anti,2)]); 
        r_anti_sorted = r_anti_norm(indx_max_anti,:); 
        B = goodcolormap('wr');
        % colormap sorted
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flipud(B'));
        imagesc(t,1:size(r_anti_sorted,1),r_anti_sorted, [0 1]);
        set(gca,'xlim',[0.05 0.350], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Instr sorted Anti ' recArea])
        % colormap unsorted
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flipud(B'));
        imagesc(t,1:size(r_anti_norm,1),r_anti_norm, [0 1]);
        set(gca,'xlim',[0.05 0.350],'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Instr Anti ' recArea])


end

        
end