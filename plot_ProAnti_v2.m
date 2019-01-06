function plot_ProAnti_v2(units,CS_modulating_cells, plotType, recArea)

% December 2018
% Input
%           load file: CS_modulating cells.
%           load file: units  - output from extractWholeNeuronResults.m
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
            set (gca, 'xlim',([-0.05 0.5]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            %           xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none'); 
            set(s_pro.patch, 'FaceAlpha', 0.05); set(s_anti.patch, 'FaceAlpha', 0.05);
            vline(0, 'k-'); box off;
            title(['Cell ' num2str(units_signif(j))]);
            
            % plot zstat bino dist
            figure(2); subplot(p(1),p(2),j); hold on;
            stat_instr(i,:) = instr_signif(j).stats.instr.pval.pbDist_testStat;
            t_instr = instr_signif(j).pro.neural.instr.ts_pst;
            position_stat_sign_instr(i,:) = abs(stat_instr(i,:))>=1.96;
            plot(t_instr, stat_instr(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.05 0.5],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
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
        for j=1:length(indx)
            t= cs_units(j).pro.neural.instr.ts_pst; % time
            r_pro= cs_units(j).pro.neural.instr.rate_pst; % psth
            std_pro = std(cs_units(j).pro.neural.instr.rate_pst);
            std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            sem_pro = std(cs_units(j).pro.neural.instr.rate_pst)/sqrt(length(cs_units(i).pro.neural.trial));
            sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
            r_anti = cs_units(j).anti.neural.instr.rate_pst;
            std_anti = std(cs_units(j).anti.neural.instr.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= cs_units(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.05 0.5]), 'TickDir', 'out', 'FontSize',12); % analysis window size
            %xlabel('Time (s)'); ylabel ('FR (spk/s)');
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
            set(gca, 'xlim',[-0.05 0.5],'ylim',[-4 4], 'TickDir', 'out', 'FontSize', 18); %ylabel('z-stat');
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]); box off;
            %             if j ==1
            %                 xlabel('Time (s)'); ylabel ('z-stat');
            %             end
            
            %       waitforbuttonpress; close all;
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
            sem_pro = std(cs_units(j).pro.neural.sacc.rate_pst)/sqrt(length(cs_units(i).pro.neural.trial));
            sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
            r_anti = cs_units(j).anti.neural.sacc.rate_pst;
            std_anti = std(cs_units(j).anti.neural.sacc.rate_pst);
            std_anti = repmat(std_anti,[1 size(r_anti,2)]);
            mean_base= cs_units(j).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            % plot psth
            figure(1);subplot(p(1),p(2),j); hold on;
            shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
            shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.1 0.2]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            % xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            vline(0, 'k-'); box off;
            %             annotation('textbox',...
            %                 [0.148482014388489 0.147388059701493 0.0829328537170263 0.0671641791044775],...
            %                 'String',['ProVSAnti = ' num2str(indx_signif(j))],...
            %                 'FitBoxToText','on');
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            
            % plot zstat bino dist
            figure(2); subplot(p(1),p(2),j); hold on;
            stat_sacc(i,:) = cs_units(j).stats.sacc.pval.pbDist_testStat;
            t_instr = cs_units(j).pro.neural.instr.ts_pst;
            position_stat_sign_instr(i,:) = abs(stat_sacc(i,:))>=1.96;
            plot(t_instr, stat_sacc(i,:), 'k','MarkerSize', 15);
            hline(-1.96, '--k');hline(1.96, '--k');vline(0, 'c');
            set(gca, 'xlim',[-0.1 0.2],'ylim',[-6 6],'xTick',[], 'TickDir', 'out', 'FontSize', 18)
            % ylabel('z-stat'); title ('Binomial pb dist'), box off;
            title(['Cell ' num2str(indx(j)) ' sig ' num2str(indx_signif(j))]);
            %             waitforbuttonpress; close all;
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
        shadedErrorBar(t, mean_delta_pro, sem_delta_pro, 'lineprops','r');
        shadedErrorBar(t, mean_delta_anti, sem_delta_anti, 'lineprops','g');
        set(gca,'TickDir', 'out', 'FontSize', 18,'ylim',[-5 15],'ytick',[-6 0 15]);
        xlabel('Time (s)'); ylabel('Change in firing rate (Hz)')
        title('CS mod cells')
        
        %stat over time
        [h_change,p_change] = ttest(abs(delta_anti_base), abs(delta_pro_base));
        h_change(h_change==0)=NaN;
        % diff anti-pro change in FR sacc
        figure; hold on;
        plot(t,nanmean(abs(delta_anti_base)-abs(delta_pro_base)),'Color','k', 'LineWidth', 2);
        plot(t,h_change*0.5, '*c')
        set(gca,'TickDir','out','ylim',[-5 5], 'ytick',[-5 0 5], 'FontSize', 18)
        xlabel('Time (s)'); ylabel('Abs change in FR anti-pro')
        title('CS mod cells')
        
        
        % plot max change in FR
        % plot scatter for significantly diff cells from mean
        figure; hold on;
        plot(max_delta_pro_base,max_delta_anti_base, '.k','MarkerSize', 18);
        %set(gca,'XScale','Log','YScale','Log' ,'FontSize', 18, 'TickDir', 'out');axis ([1e0 1e2 1e0 1e2]);
        set(gca,'FontSize', 18, 'TickDir', 'out', 'ylim',[0 30],'xlim',[0 30]); %axis ([1e0 1e2 1e0 1e2]);
        % plot([1e0 1e2],[1e0 1e2]);
        plot([0 30], [0 30]);
        xlabel('Max change pro'); ylabel('Max change anti');
        title(['Max change in firing rate >> ' recArea])
        axis square
        [h,p] = ttest(max_delta_pro_base,max_delta_anti_base)
        
    case 'bino_dist_CS_mod'
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
        end
        t_sacc = cs_units(1).pro.neural.sacc.ts_pst;
        
        % plot
        figure; hold on;
        plot(t_sacc,abs(stat_sacc), '.', 'MarkerSize', 18);
        set(gca, 'xlim',[-0.1 0.3], 'TickDir', 'out', 'FontSize', 18);
        hline(1.96); vline(0);xlabel('time(s)'); ylabel('Z-stat')
        %histogram
        figure;histogram(abs(stat_sacc),25);
        set(gca, 'XDir','reverse','xlim',[0 5],'TickDir', 'out', 'FontSize', 18);
        vline(1.96); box off;
        
        figure; hold on;
        plot(t_sacc,smooth(nanmean(abs(stat_sacc)),3),'LineWidth', 2,'Color','k'); % smoothed
        hline(1.96,'k')
        set(gca, 'xlim',[0 0.3],'ylim',[0 2], 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Sacc => ' recArea ' n= ' num2str(numel(indx))]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        % colormap
        colormap(bone);
        imagesc(t_sacc,1:size(stat_sacc,1),abs(stat_sacc), [0 1.96]);
        set(gca,'xlim',[-0.1 0.3],'YTickLabel',[], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off;
        signif = false(size(stat_sacc,1),1);
        for i=1:size(stat_sacc,1)
            if sum(stat_sacc(i,:)>1.96)>=2
                signif(i) = 1;
            end
        end
        signif = find(signif);
        figure('Position',[1160 678 80 420]); plot(1,signif, '*k'); box off;
        set(gca, 'xlim',[0.9 1.1], 'ylim',[0 25], 'YDir', 'reverse', 'xTick', []);
        
    case 'bino_dist_CS_mod_instr'
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
            position_stat_sign_sacc(j,:) = abs(stat_instr(j,:))>=1.96;
        end
        t_instr = cs_units(1).pro.neural.instr.ts_pst;
        
        % plot
        figure; hold on;
        plot(t_instr,abs(stat_instr), '.', 'MarkerSize', 18);
        set(gca, 'xlim',[0 0.5], 'TickDir', 'out', 'FontSize', 18);
        hline(1.96); vline(0);xlabel('time(s)'); ylabel('Z-stat')
        %histogram
        figure;histogram(abs(stat_instr),25);
        set(gca, 'XDir','reverse','xlim',[0 5],'TickDir', 'out', 'FontSize', 18);
        vline(1.96); box off;
        
        figure; hold on;
        plot(t_instr,smooth(nanmean(abs(stat_instr)),3),'LineWidth', 2,'Color','k'); % smoothed
        hline(1.96,'k')
        set(gca, 'xlim',[0 0.5],'ylim',[0 2], 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Sacc => ' recArea ' n= ' num2str(numel(indx))]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        % colormap
        colormap(bone);
        imagesc(t_instr,1:size(stat_instr,1),abs(stat_instr), [0 1.96]);
        set(gca,'xlim',[0 0.5],'YTickLabel',[], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('Neuron'); box off;
        signif = false(size(stat_instr,1),1);
        for i=1:size(stat_instr,1)
            if sum(stat_instr(i,:)>1.96)>=2
                signif(i) = 1;
            end
        end
        signif = find(signif);
        figure('Position',[1160 678 80 420]); plot(1,signif, '*k'); box off;
        set(gca, 'xlim',[0.9 1.1], 'ylim',[0 25], 'YDir', 'reverse', 'xTick', []);
        
        
        
end