function plot_ProAnti(units, plotType, cellNum, recArea)

% Needed: units.mat
% Input:    units  - output from extractWholeNeuronResults.m
%           cellNum - cell you want to plot. If all leave empty = []
%           plotType - plot you want (e.g. raster, psth, etc.)
%           recArea - OMV, lateral Cb

% if running population (all neurons), just inpu t [] in cellNum.
% if running single neuron, on recArea input []. Example: plot_ProAnti(units, 'raster_sacc',1,  [])

% List of plots

% 'eyeKin' : plot eye kinematics for all cells and trials
% 'raster_sacc': saccade aligned raster plot for the chosen cell
% 'raster_instr': instruction aligned raster plot for the chosen cell
% 'psth': psth for the chosen cell, aligned to saccade and instruction
% 'pplsth_sacc_all': psth aligned to saccade onset for all cells. Press any key to plot next cell
% 'psth_instr_all':psth aligned to instruction onset for all cells. Press
% any key to plot next cell
% 'delta_rate': plots time-course of change in firing rate (firing rate-baseline) for all cells and then plots max change in a scatter plot for significant cells.
% 'DDI'
% 'indx_change'
% 'colormap_sacc': plot time-course of normalized firing rate for all cells
% 'waterfall_sacc': waterfall plot of time-course of normalized firing rate for all cells
% 'scatter_pro_anti': scatter plot of nspk of pro vs anti for a cell
% 'scatter_pro_anti_pop': scatter plot of mean firing rates of all cells pro vs anti
% 'peak_resp_instr': scatter plot of peak response instr pro vs anti
% 'peak_resp_sacc': scatter plot of peak response sacc pro vs anti
% 'firingVSamp': firing rate vs sacc amplitude for a cell and also per sacc amplitude blocks
% 'firingVSvel': firing rate vs sacc peak vel
% 'firingVSdur': firing rate vs sacc duration
% 'binomial_pb_dist': Binomial probability distribution across time aligned to instr and sacc. >1.96 is significant.
% 'spk_pb_stat': spike probability density for 50 ms
% 'cv'
% 'cv2'
% 'eye move' plot eye traces for all trials for one neuron
% 'mod_ratio': plot modulation ratio for all neurons and for signif diff neurons

%%
%%

switch plotType
    case 'eye_kin'
        %% plot
        % gather
        eyeKin = eyeKinematics_ProAnti(units);
        
        proAmp = vertcat(eyeKin(1,:).proAmp);
        antiAmp = vertcat(eyeKin(1,:).antiAmp);
        
        proDur = vertcat(eyeKin(1,:).proDur);
        antiDur = vertcat(eyeKin(1,:).antiDur);
        
        proPV = vertcat(eyeKin(1,:).proPV);
        antiPV = vertcat(eyeKin(1,:).antiPV);
        
        proRT = vertcat(eyeKin(1,:).proRT)*1000;
        antiRT = vertcat(eyeKin(1,:).antiRT)*1000;
        
        % plot amp
        figure; hold on
        h1 = histfit(proAmp,20,'kernel');
        h2 = histfit(antiAmp,20,'kernel');
        xlabel('Saccade amplitude (deg)')
        ylabel('Number of trials')
        set (gca, 'TickDir', 'out','FontSize', 18);
        alpha(0.25)
        set(h1(1),'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
        set(h2(1),'FaceColor', [0 1 0],'EdgeColor', [0 1 0]);
        set(h1(2),'Color',[1 0 0]);
        set(h2(2),'Color',[0 1 0]);
        vline(mean(proAmp),'r'); % draw line on mean
        vline(mean(antiAmp),'g'); % draw line on mean
        
        %plot dur
        figure; hold on
        h1 = histfit(proDur,20,'kernel');
        h2 = histfit(antiDur,20,'kernel');
        xlabel('Saccade duration (ms)')
        ylabel('Number of trials')
        set (gca, 'TickDir', 'out','FontSize', 18);
        alpha(0.25)
        set(h1(1),'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
        set(h2(1),'FaceColor', [0 1 0],'EdgeColor', [0 1 0]);
        set(h1(2),'Color',[1 0 0]);
        set(h2(2),'Color',[0 1 0]);
        vline(mean(proDur),'r'); % draw line on mean
        vline(mean(antiDur),'g'); % draw line on mean
        
        % plot PV
        figure; hold on
        h1 = histfit(proPV,20,'kernel');
        h2 = histfit(antiPV,20,'kernel');
        xlabel('Saccade peak velocity (deg/s)')
        ylabel('Number of trials')
        set (gca, 'TickDir', 'out','FontSize', 18);
        alpha(0.25)
        set(h1(1),'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
        set(h2(1),'FaceColor', [0 1 0],'EdgeColor', [0 1 0]);
        set(h1(2),'Color',[1 0 0]);
        set(h2(2),'Color',[0 1 0]);
        vline(mean(proPV),'r'); % draw line on mean
        vline(mean(antiPV),'g'); % draw line on mean
        
        % plot RT
        figure; hold on
        h1 = histfit(proRT,20,'kernel');
        h2 = histfit(antiRT,20,'kernel');
        xlabel('Reaction time (ms)')
        ylabel('Number of trials')
        set (gca, 'TickDir', 'out','FontSize', 18);
        alpha(0.25)
        set(h1(1),'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
        set(h2(1),'FaceColor', [0 1 0],'EdgeColor', [0 1 0]);
        set(h1(2),'Color',[1 0 0]);
        set(h2(2),'Color',[0 1 0]);
        vline(mean(proRT),'r'); % draw line on mean
        vline(mean(antiRT),'g'); % draw line on mean
        
        
    case 'raster_sacc'
        % saccade aligned
        % pro
        [~,indx] = sort([units(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).pro.behav.trial(indx).reactionTime];
        r_pro= units(cellNum).pro.neural.trial;
        recArea = units(cellNum).area;
        
        figure;subplot (2,1,1); hold on;box off
        
        if strcmp(units(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_SS_align_sacc)
                    %plot(sorted_RT(j),j,'.r');
                    %plot(r_pro(indx(j)).tspk_SS_align_sacc(1:2:end),j,'.k'); %plot every n spikes
                    plot(r_pro(indx(j)).tspk_SS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.150 0.151]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            %title(['Pro (aligned to saccade) SS => ' recArea ' cellNum ' num2str(cellNum)]); xlabel('Time (s)');ylabel('Trial Num')
        else
            for j= 1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_CS_align_sacc)
                    %plot(sorted_RT(j),j,'.r');
                    plot(r_pro(indx(j)).tspk_CS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.150 0.151]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            %title(['Pro (aligned to saccade) CS  => ' recArea ' cellNum ' num2str(cellNum) ]); xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
        % anti
        [~,indx] = sort([units(cellNum).anti.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).anti.behav.trial(indx).reactionTime];
        r_anti=units(cellNum).anti.neural.trial;
        subplot (2,1,2); hold on;box off
        
        if strcmp(units(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_SS_align_sacc)
                    %plot(sorted_RT(j),j,'.r');
                    %plot(r_anti(indx(j)).tspk_SS_align_sacc(1:2:end),j,'.k'); %plot every n spikes
                    plot(r_anti(indx(j)).tspk_SS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.150 0.151]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            %title(['Anti (aligned to saccade) SS  => ' recArea]); xlabel('Time (s)');ylabel('Trial Num')
        else
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_CS_align_sacc)
                    %plot(sorted_RT(j),j,'.r');
                    plot(r_anti(indx(j)).tspk_CS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.150 0.151]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            %title(['Anti (aligned to saccade) CS  => ' recArea]); xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
    case 'raster_instr'
        % instr aligned
        % pro
        
        [~,indx] = sort([units(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).pro.behav.trial(indx).reactionTime];
        r_pro= units(cellNum).pro.neural.trial;
        recArea = units(cellNum).area;
        
        figure;subplot (2,1,1); hold on;box off
        
        if strcmp(units(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_SS)
                    plot(r_pro(indx(j)).tspk_SS(1:2:end),j,'.k'); %plot every n spikes
                    %plot(r_pro(indx(j)).tspk_SS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.3]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title(['Pro (aligned to instruction) SS  => ' recArea]); xlabel('Time (s)');ylabel('Trial Num')
        else
            for j=1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_CS)
                    plot(r_pro(indx(j)).tspk_CS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.3]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title(['Pro (aligned to instruction) CS  => ' recArea]); xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
        % anti
        [~,indx] = sort([units(cellNum).anti.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).anti.behav.trial(indx).reactionTime];
        r_anti=units(cellNum).anti.neural.trial;
        subplot (2,1,2); hold on;box off
        
        if strcmp(units(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_SS)
                    plot(r_anti(indx(j)).tspk_SS(1:2:end),j,'.k'); %plot every n spikes
                    %plot(units(cellNum).anti.neural.trial(indx).tspk_SS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.3]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title(['Anti (aligned to instruction) SS  => ' recArea]); xlabel('Time (s)');ylabel('Trial Num')
        else
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_CS)
                    plot(r_anti(indx(j)).tspk_CS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.3]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title(['Anti (aligned to saccade) CS  => ' recArea]); xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
        
    case 'psth'
        
        %gather
        t= units(cellNum).pro.neural.sacc.ts_pst;
        r_pro= units(cellNum).pro.neural.sacc.rate_pst;
        sem_pro = std(units(cellNum).pro.neural.sacc.rate_pst)/sqrt(length(units(cellNum).pro.neural.trial));
        sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
        std_pro = std(units(cellNum).pro.neural.sacc.rate_pst); std_pro = repmat(std_pro,[1 size(r_pro,2)]);
        r_anti = units(cellNum).anti.neural.sacc.rate_pst;
        sem_anti = std(units(cellNum).anti.neural.sacc.rate_pst)/sqrt(length(units(cellNum).anti.neural.trial));
        sem_anti = repmat(sem_anti,[1 size(r_anti,2)]); 
        std_anti = std(units(cellNum).anti.neural.sacc.rate_pst); std_anti = repmat(std_anti,[1 size(r_anti,2)]);
        mean_base= units(cellNum).pro.neural.base.rate_mu;
        mean_base = repmat(mean_base,[1 size(r_pro,2)]);
        recArea = units(cellNum).area;
        proVSanti_sacc = units(cellNum).stats.sacc.flags.proVsAnti_sacc;
        proVSanti_instr = units(cellNum).stats.instr.flags.proVsAnti_instr;
        
        %         %plot
        %         figure; hold on;
        %         plot(t, r_pro, 'r', 'LineWidth', 3);
        %         plot(t, r_anti, 'g', 'LineWidth', 3);
        %         plot(t,mean_base,'--k','LineWidth', 0.3);
        %         set (gca, 'xlim',([-0.5 0.5]), 'TickDir', 'out', 'FontSize',18);
        %         xlabel('Time (s)'); ylabel ('Firing rate (spk/s');
        %         vline(0, 'k-')
        %         box off
        
        % plot w/sem
        figure('Position', [375 403 834 536]); subplot(1,2,1); hold on
        s_pro = shadedErrorBar(t,r_pro,std_pro,'lineprops','r');
        s_anti = shadedErrorBar(t,r_anti,std_anti,'lineprops','g');
        plot(t,mean_base,'--k','LineWidth', 0.3);
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        vline(0, 'k-'); box off;
        
        set (gca, 'xlim',([-0.150 0.151]), 'TickDir', 'out', 'FontSize',18); % analysis window size
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
        vline(0, 'k-');
        box off
        title(['Aligned to saccade => ' recArea]);
        % Create textbox
        annotation('textbox',...
            [0.148482014388489 0.147388059701493 0.0829328537170263 0.0671641791044775],...
            'String',['ProVSAnti = ' num2str(proVSanti_sacc)],...
            'FitBoxToText','on');
        
        
        %% instruction period
        %gather
        
        t= units(cellNum).pro.neural.sacc.ts_pst;
        r_pro= units(cellNum).pro.neural.instr.rate_pst;
        sem_pro = std(units(cellNum).pro.neural.instr.rate_pst)/sqrt(length(units(cellNum).pro.neural.trial));
        sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
        r_anti = units(cellNum).anti.neural.instr.rate_pst; 
        sem_anti = std(units(cellNum).anti.neural.instr.rate_pst)/sqrt(length(units(cellNum).anti.neural.trial));
        sem_anti = repmat(sem_anti,[1 size(r_anti,2)]);
        proVSanti_instr = units(cellNum).stats.instr.flags.proVsAnti_instr;
        
        %         %plot
        %         figure; hold on;
        %         plot(t, r_pro, 'r', 'LineWidth', 3);
        %         plot(t, r_anti, 'g', 'LineWidth', 3);
        %         set (gca, 'xlim',([-0.5 0.5]), 'TickDir', 'out', 'FontSize',18);
        %         xlabel('Time (s)'); ylabel ('Firing rate (spk/s');
        %         vline(0, 'k--')
        %         box off
        
        % plot w/sem
        subplot(1,2,2);
        shadedErrorBar(t, r_pro,sem_pro,'lineprops','r');
        shadedErrorBar(t, r_anti,sem_anti,'lineprops','g');
        set (gca, 'xlim',([-0.1 0.3]), 'TickDir', 'out', 'FontSize',18);
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
        vline(0, 'k--');
        box off
        title(['Aligned to instruction => ' recArea]);
        annotation('textbox',...
            [0.592778202676864 0.128125 0.0859980879541108 0.0484375],...
            'String',['ProVSAnti = ' num2str(proVSanti_instr)],...
            'FitBoxToText','on');
        
    case 'pplsth_sacc_all'   % plot all neurons per area
        % gather indx
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        
        for i=1:length(indx_area)
            t= units(indx_area(i)).pro.neural.sacc.ts_pst; % time
            r_pro= units(indx_area(i)).pro.neural.sacc.rate_pst; % psth
            sem_pro = std(units(indx_area(i)).pro.neural.sacc.rate_pst)/sqrt(length(units(indx_area(i)).pro.neural.trial));
            sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
            r_anti = units(indx_area(i)).anti.neural.sacc.rate_pst;
            sem_anti = std(units(indx_area(i)).anti.neural.sacc.rate_pst)/sqrt(length(units(indx_area(i)).anti.neural.trial));
            sem_anti = repmat(sem_anti,[1 size(r_anti,2)]);
            mean_base= units(indx_area(i)).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            proVSanti_sacc = units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc;
            
            % plot w/sem
            figure; hold on;
            shadedErrorBar(t,r_pro,sem_pro,'lineprops','r');
            shadedErrorBar(t,r_anti,sem_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.1 0.2]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            vline(0, 'k-');
            box off
            title(['Aligned to saccade => ' recArea ' unit= ' num2str(indx_area(i))])
            annotation('textbox',...
                [0.159928571428571 0.154761904761905 0.150785714285714 0.104761904761905],...
                'String',['Pro VS Anti = ' num2str(proVSanti_sacc)],...
                'FitBoxToText','on');
            fname = ['psth_sacc_all_' recArea];
            print(fname,'-append', '-dpsc2')
            waitforbuttonpress; close all;
        end
    case 'psth_instr_all'
        % gather indx
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        for i=1:length(indx_area)
            t= units(indx_area(i)).pro.neural.instr.ts_pst; % time
            r_pro= units(indx_area(i)).pro.neural.instr.rate_pst; % psth
            sem_pro = std(units(indx_area(i)).pro.neural.instr.rate_pst)/sqrt(length(units(indx_area(i)).pro.neural.trial));
            sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
            r_anti = units(indx_area(i)).anti.neural.instr.rate_pst;
            sem_anti = std(units(indx_area(i)).anti.neural.instr.rate_pst)/sqrt(length(units(indx_area(i)).anti.neural.trial));
            sem_anti = repmat(sem_anti,[1 size(r_anti,2)]);
            mean_base= units(indx_area(i)).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            proVSanti_instr = units(indx_area(i)).stats.instr.flags.proVsAnti_instr;
            
            % plot w/sem
            figure; hold on;
            shadedErrorBar(t,r_pro,sem_pro,'lineprops','r');
            shadedErrorBar(t,r_anti,sem_anti,'lineprops','g');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.1 0.3]), 'TickDir', 'out', 'FontSize',18); % analysis window size
            xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
            vline(0, 'k-');
            box off
            title(['Aligned to instruction => ' recArea ' unit= ' num2str(indx_area(i))])
            annotation('textbox',...
                [0.159928571428571 0.154761904761905 0.150785714285714 0.104761904761905],...
                'String',['Pro VS Anti = ' num2str(proVSanti_instr)],...
                'FitBoxToText','on');
            waitforbuttonpress; close all;
        end
        
    case 'delta_rate'
        % Plot change in FR from baseline
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        % gather pro
        t = units(1).pro.neural.sacc.ts_pst_win;
        t_instr = units(1).pro.neural.instr.ts_pst_win;
        for cells = 1:length(indx_area)
            max_delta_pro(cells) = max(abs(units(indx_area(cells)).pro.neural.sacc.delta_rate));
            max_delta_pro_base(cells) = max(abs(units(indx_area(cells)).pro.neural.sacc.delta_rate_base));
            delta_pro(cells,:)=units(indx_area(cells)).pro.neural.sacc.delta_rate;
            delta_pro_base(cells,:)=units(indx_area(cells)).pro.neural.sacc.delta_rate_base;
            delta_pro_base_instr(cells,:)=units(indx_area(cells)).pro.neural.instr.delta_rate_base;
            r_pro(cells,:) = units(indx_area(cells)).pro.neural.sacc.rate_pst_win;
            sig_pro(cells,:) = std(units(indx_area(cells)).pro.neural.sacc.rate_pst_win)/sqrt(sum(indx_area));
            std_pro(cells,:) = std(units(indx_area(cells)).pro.neural.sacc.rate_pst_win);
        end
        % plot(t, mean(abs(delta_pro)'));
        % plot for pro all cells
        
        
        % gather anti
        for cells = 1:length(indx_area)
            max_delta_anti(cells) = max(abs(units(indx_area(cells)).anti.neural.sacc.delta_rate));
            max_delta_anti_base(cells) = max(abs(units(indx_area(cells)).anti.neural.sacc.delta_rate_base));
            delta_anti(cells,:)=units(indx_area(cells)).anti.neural.sacc.delta_rate;
            delta_anti_base(cells,:)=units(indx_area(cells)).anti.neural.sacc.delta_rate_base;
            delta_anti_base_instr(cells,:)=units(indx_area(cells)).anti.neural.instr.delta_rate_base;
            r_anti(cells,:) = units(indx_area(cells)).anti.neural.sacc.rate_pst_win;
            sig_anti(cells,:) = std(units(indx_area(cells)).anti.neural.sacc.rate_pst_win)/sqrt(sum(indx_area));
            std_anti(cells,:) = std(units(indx_area(cells)).anti.neural.sacc.rate_pst_win);
            instr_pro_base(cells,:) = units(indx_area(cells)).pro.neural.instr.delta_rate_base;
            instr_anti_base(cells,:) = units(indx_area(cells)).anti.neural.instr.delta_rate_base;
        end
        %plot(t, delta_anti');
        
        % get significantly different cells
        for i = 1:length(indx_area)
            indx_instr(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr);
            indx_instr_ks_npsk(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr_ks_nspk);
            indx_instr_ks_tpsk(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr_ks_t_spk);
            indx_sign(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
            indx_sacc_ks_nspk(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk);
            indx_sacc_ks_tspk(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc_ks_t_spk);
        end
        nunits = size(abs(delta_pro(indx_sign,:)),2);
        
        % quick plot of psth of mean of all cells
        figure; hold on;
        plot(t,mean(r_pro));
        plot(t,mean(r_anti));
        shadedErrorBar(t, mean(r_pro), repmat(mean(sig_pro),[size(mean(r_pro)) 1]), 'lineprops','r');
        shadedErrorBar(t, mean(r_anti),repmat(mean(sig_anti),[size(mean(r_anti)) 1]), 'lineprops','g');
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [60 70], 'ytick', [60 61 63],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Firing rate (spk/s)'); xlabel('Time (s)')
        
        
        
        %plot change in FR from mean for all cells
        figure; hold on;
        %pro
        plot(t, mean(abs(delta_pro)', 2),'r', 'LineWidth',2);
        set(gca,'TickDir', 'out', 'FontSize', 18)
        % anti
        plot(t, mean(abs(delta_anti)', 2),'g', 'LineWidth',2);
        set(gca,'TickDir', 'out', 'FontSize', 18)
        box off
        xlabel('time from saccade onset') ;ylabel('Change in firing (spks/s)')
        title([recArea ' (all cells)'])
        
        %         % abs vals
        %         mean_delta_pro = mean(abs(delta_pro_base));
        %         sem_delta_pro = std(abs(delta_pro_base))/sqrt(sum(length(indx_area)));
        %         mean_delta_anti = mean(abs(delta_anti_base));
        %         sem_delta_anti = std(abs(delta_anti_base),0,1)/sqrt(sum(length(indx_area)));
        
        % sacc
        mean_delta_pro = mean(abs(delta_pro_base));
        sem_delta_pro = std(delta_pro_base)/sqrt(sum(length(indx_area)));
        mean_delta_anti = mean(abs(delta_anti_base));
        sem_delta_anti = std(delta_anti_base,0,1)/sqrt(sum(length(indx_area)));
        %instr
        mean_delta_pro_instr = mean(delta_pro_base_instr(indx_instr,:),1);
        sem_delta_pro_instr = std(delta_pro_base_instr(indx_instr,:))/sqrt(sum(indx_instr));
        mean_delta_anti_instr = mean(delta_anti_base_instr(indx_instr,:),1);
        sem_delta_anti_instr = std(delta_anti_base_instr(indx_instr,:))/sqrt(sum(indx_instr));
        
        %plot change in FR from baseline for all cells sacc
        figure; hold on;
        s_pro = shadedErrorBar(t, mean_delta_pro, sem_delta_pro, 'lineprops','r');
        s_anti = shadedErrorBar(t, mean_delta_anti, sem_delta_anti, 'lineprops','g');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        vline(0, 'k-'); box off;
        set(gca,'TickDir', 'out', 'FontSize', 18,'ylim',[-3 19], 'ytick', [-3 2], 'xlim', [-0.150 0.151]); box off
        xlabel('Time (s)') ;ylabel('Change in firing from base (spks/s)')
        title([recArea ' (all cells)'])
        
        %plot change in FR from baseline for all cells instr
        figure; hold on;
        shadedErrorBar(t_instr, mean_delta_pro_instr, sem_delta_pro_instr, 'lineprops','r');
        shadedErrorBar(t_instr, mean_delta_anti_instr, sem_delta_anti_instr, 'lineprops','g');
        set(gca,'TickDir', 'out', 'FontSize', 18,'ylim',[-5 20], 'ytick', [0 5]); box off
        xlabel('Time (s)') ;ylabel('Change in firing from base (Instr) (spks/s)')
        title([recArea ' (all cells)'])
        
        
        %plot mean FR for signif cells instr
        figure; hold on;
        plot(t,mean(r_pro(indx_sign)));
        plot(t,mean(r_anti(indx_sign)));
        shadedErrorBar(t, mean(r_pro(indx_sign,:)), repmat(mean(sig_pro(indx_sign,:)),[size(mean(r_pro(indx_sign,:))) 1]), 'lineprops','r');
        shadedErrorBar(t, mean(r_anti(indx_sign,:)),repmat(mean(sig_anti(indx_sign,:)),[size(mean(r_anti(indx_sign,:))) 1]), 'lineprops','g');
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [60 70], 'ytick', [60 61 62 70],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Firing rate (spk/s)'); xlabel('Time (s)')
        title('signif cells')
        
        % plot change in FR for significantly different cells
        % Abs vals
        %         mean_delta_pro = mean(abs(delta_pro_base(indx_sign,:)),1);
        %         sem_delta_pro = std(abs(delta_pro_base(indx_sign,:)),0,1)/sqrt(sum(indx_sign));
        %         mean_delta_anti = mean(abs(delta_anti_base(indx_sign,:)),1);
        %         sem_delta_anti = std(abs(delta_anti_base(indx_sign,:)),0,1)/sqrt(sum(indx_sign));
        
        mean_delta_pro = mean(delta_pro_base(indx_sign,:),1);
        sem_delta_pro = std(delta_pro_base(indx_sign,:),0,1)/sqrt(sum(indx_sign));
        mean_delta_anti = mean(delta_anti_base(indx_sign,:),1);
        sem_delta_anti = std(delta_anti_base(indx_sign,:),0,1)/sqrt(sum(indx_sign));
        
        % all
        figure; hold on;
        plot(t, abs(delta_pro_base(indx_sign,:)));
        set(gca, 'xlim',[-0.150 0.150],'TickDir', 'out', 'FontSize', 18)
        xlabel('Time (s) aligned to saccade onset'); ylabel('Change in FR (Hz)')
        title(['All sign cells - prosaccade ' num2str(sum(indx_sign))]);
        
        figure; hold on;
        plot(t, abs(delta_anti_base(indx_sign,:)));
        set(gca,'xlim',[-0.150 0.150],'TickDir', 'out', 'FontSize', 18)
        xlabel('Time (s) aligned to saccade onset'); ylabel('Change in FR (Hz)')
        title(['All sign cells - antisaccade ' num2str(sum(indx_sign))]);
        
        %stat over time
        [h_change,p_change] = ttest(abs(delta_anti(indx_sign,:)), abs(delta_pro(indx_sign,:)));
        [h_instr,p_instr] = ttest(abs(instr_anti_base(indx_sign,:)), abs(instr_pro_base(indx_sign,:)));
        
        % diff anti-pro change in FR sacc
        figure; hold on;
        plot(t,nanmean(abs(delta_anti_base(indx_sign,:))-abs(delta_pro_base(indx_sign,:))),'Color','k', 'LineWidth', 2);
        % plot(t,h_change*0.5, '*c')
        set(gca,'xlim', [-0.150 0.151],'TickDir','out','ylim',[-1 4], 'ytick',[-1 4], 'FontSize', 18)
        xlabel('Time (s)'); ylabel('Abs change in FR anti-pro')
        title('Signif diff cells')
        
        % diff anti-pro change in FR instr
        figure; hold on;
        plot(t_instr,nanmean(abs(instr_anti_base(indx_sign,:))-abs(instr_pro_base(indx_sign,:))),'Color','k', 'LineWidth', 2);
        %plot(t_instr,h_instr*0.5, '*c')
        set(gca,'TickDir','out','ylim',[-1 2], 'ytick',[-1 2], 'FontSize', 18)
        xlabel('Time (s)'); ylabel('Abs change in FR anti-pro instr')
        title('Signif diff cells')
        
        
        figure; hold on;
        s_pro = shadedErrorBar(t, mean_delta_pro, sem_delta_pro, 'lineprops','r');
        s_anti = shadedErrorBar(t, mean_delta_anti, sem_delta_anti, 'lineprops','g');
        set(gca,'TickDir', 'out', 'FontSize', 18,'ytick',[0 20], 'xlim', [-0.150 0.151]);
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        vline(0, 'k-'); box off;
        xlabel('Time (s)'); ylabel('Change in firing rate (Hz)')
        title('Only sign diff cells')
        
        % plot scatter for significantly diff cells from mean
%         figure; hold on;
%         plot(max_delta_pro,max_delta_anti, '.k','MarkerSize', 18);
%         plot(max_delta_pro(indx_sign),max_delta_anti(indx_sign), '.c','MarkerSize', 18);
%         set(gca,'XScale','Log','YScale','Log' ,'FontSize', 18, 'TickDir', 'out');axis ([1e0 1e2 1e0 1e2]);
%         plot([1e0 1e2],[1e0 1e2]);
%         xlabel('Max change pro'); ylabel('Max change anti');
%         title(['Max change in firing rate >> ' recArea])
%         axis square
%         [h,p] = ttest(max_delta_pro,max_delta_anti)
        
        % plot scatter for significantly diff cells from baseline
        figure; hold on;
        plot(max_delta_pro_base,max_delta_anti_base, '.k','MarkerSize', 18);
        plot(max_delta_pro_base(indx_sign),max_delta_anti_base(indx_sign), '.c','MarkerSize', 18);
        set(gca,'XScale','Log','YScale','Log' ,'FontSize', 18, 'TickDir', 'out');axis ([1e0 1e2 1e0 1e2]);
        plot([1e0 1e2],[1e0 1e2]);
        xlabel('Max change pro'); ylabel('Max change anti');
        title(['Max change in firing rate from base >> ' recArea])
        axis square
        [h,p] = ttest(max_delta_pro_base,max_delta_anti_base)
        
    case 'DDI_sacc'
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        % extract
        for j=1:length(nunits)
            ddi_pro(j,:) = units(indx_area(j)).stats.pro.sacc.DDI;
            ddi_anti(j,:) = units(indx_area(j)).stats.anti.sacc.DDI;
        end
        [y_pro,bin_pro] = hist(ddi_pro,15);
        [y_anti,bin_anti] = hist(ddi_anti,15);
        
        % pro
        figure; hold on;
        plot(bin_pro,y_pro, 'r','LineWidth',2);
        ha1 = area(bin_pro,y_pro,'FaceColor',[1 0 0],'Linewidth',2);
        set(ha1,'FaceAlpha',0.5,'EdgeColor','none'); box off;
        set(gca,'xlim', [0 1],'yTick',[0 5],'xTick',[0 0.5 1], 'TickDir', 'out','FontSize', 18);
        xlabel('DDI'); ylabel('Nr of neurons');vline(mean(ddi_pro));
        
        % anti
        figure; hold on;
        plot(bin_anti,y_anti, 'g','LineWidth',2);
        ha1 = area(bin_anti,y_anti,'FaceColor',[0 1 0],'Linewidth',2);
        set(ha1,'FaceAlpha',0.5,'EdgeColor','none'); box off;
        set(gca,'xlim', [0 1],'yTick',[0 5],'xTick',[0 0.5 1], 'TickDir', 'out','FontSize', 18);
        xlabel('DDI'); ylabel('Nr of neurons');vline(mean(ddi_anti));
        
        % scatter
        figure; hold on;
        plot(ddi_pro,ddi_anti,'.k','MarkerSize', 30);
        plot([0:1],[0:1],'k');
        set (gca,'xlim', [0.5 1],'ylim',[0.5 1],'xTick', [0 0.5 1], 'yTick', [0 0.5 1], 'TickDir', 'out','FontSize', 18);
        axis square; box off;
        xlabel('Pro'); ylabel('Anti'); title('DDI')
        [h,p] = ttest(ddi_pro,ddi_anti)
        
        % diagonal
        figure; hold on;
        histogram(ddi_pro-ddi_anti,25)
        set(gca,'ylim', [0 12],'xlim',[-0.707 0.707] ,'TickDir', 'out','ytick', [0 12],'Fontsize', 18); xlabel('diagonal', 'Fontsize', 18); % 1/sqrt(2) 0.702
        vline(0)
        
    case 'DDI_instr'
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        % extract
        for j=1:length(nunits)
            ddi_pro(j,:) = units(indx_area(j)).stats.pro.instr.DDI;
            ddi_anti(j,:) = units(indx_area(j)).stats.anti.instr.DDI;
        end
        [y_pro,bin_pro] = hist(ddi_pro,15);
        [y_anti,bin_anti] = hist(ddi_anti,15);
        
        % pro sacc
        figure; hold on;
        plot(bin_pro,y_pro, 'r','LineWidth',2);
        ha1 = area(bin_pro,y_pro,'FaceColor',[1 0 0],'Linewidth',2);
        set(ha1,'FaceAlpha',0.5,'EdgeColor','none'); box off;
        set(gca,'xlim', [0 1],'yTick',[0 5],'xTick',[0 0.5 1], 'TickDir', 'out','FontSize', 18);
        xlabel('DDI'); ylabel('Nr of neurons');vline(mean(ddi_pro));
        
        % anti sacc
        figure; hold on;
        plot(bin_anti,y_anti, 'g','LineWidth',2);
        ha1 = area(bin_anti,y_anti,'FaceColor',[0 1 0],'Linewidth',2);
        set(ha1,'FaceAlpha',0.5,'EdgeColor','none'); box off;
        set(gca,'xlim', [0 1],'yTick',[0 5],'xTick',[0 0.5 1], 'TickDir', 'out','FontSize', 18);
        xlabel('DDI'); ylabel('Nr of neurons');vline(mean(ddi_anti));
        
        % scatter
        figure; hold on;
        plot(ddi_pro,ddi_anti,'.k','MarkerSize', 30);
        plot([0:1],[0:1],'k');
        set (gca,'xlim', [0.5 1],'ylim',[0.5 1],'xTick', [0 0.5 1], 'yTick', [0 0.5 1], 'TickDir', 'out','FontSize', 18);
        axis square; box off;
        xlabel('Pro'); ylabel('Anti'); title('DDI')
        [h,p] = ttest(ddi_pro,ddi_anti)
        
        % diagonal
        figure; hold on;
        histogram(ddi_pro-ddi_anti,25)
        set(gca,'ylim', [0 12],'xlim',[-0.707 0.707] ,'TickDir', 'out','ytick', [0 12],'Fontsize', 18); xlabel('diagonal', 'Fontsize', 18); % 1/sqrt(2) 0.702
        vline(0)
        
        
    case 'indx_change'
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        % extract
        for j=1:length(nunits)
            change_instr(j,:) = units(indx_area(j)).stats.instr.change_indx;
            change_sacc(j,:) = units(indx_area(j)).stats.sacc.change_indx;
        end
        
        %plot
        figure; hold on;
        histogram(change_instr,50);
        set(gca,'ylim', [0 20],'xlim',[-0.2 0.2] ,'TickDir', 'out','ytick', [0 25],'Fontsize', 18);
        
        figure; hold on;
        histogram(change_sacc, 25)
        set(gca,'ylim', [0 20],'xlim',[-0.25 0.25] ,'TickDir', 'out','ytick', [0 25],'Fontsize', 18);
        
        
    case 'colormap_sacc'
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        %pro
        t = units(1).pro.neural.sacc.ts_pst; clear r;
        for j=1:length(nunits)
            r_pro(j,:)= units(indx_area(j)).pro.neural.sacc.rate_pst;
            r_anti(j,:)= units(indx_area(j)).anti.neural.sacc.rate_pst;
            r(j,:) = units(indx_area(j)).pro.neural.sacc.norm_rate_pst(1,:);
        end
        
        [r,t] = smooth_colormap(r,t);
        
        
        
        %         % sort
        %                 t_sacc = t>-0.2 & t<=0.255;
        %                 for i=1:length(r(:,1))
        %                     this_r = r(i,:);
        %                     this_r_sacc = this_r(t_sacc);
        %                     [~,max_r(i)]= max(this_r_sacc)
        %                    if max_r(i) == 1
        %                    max_r(i) = max_r(i)+i*0.01
        %                    end
        %                 end
        
        
        % plot colormap
        B = goodcolormap('bwr');
        figure; set(gcf,'Position',[100 200 300 300]);
        % hold on; colormap(B');
        % hold on; colormap(bluewhitered);
        hold on; colormap(bone(4));
        imagesc(t,1:nunits,r,[0,1]);
        %imagesc(t,1:nunits,r_sorted,[0,1]);
        set(gca,'xlim',[-0.2 0.255],'ylim',[1 nunits(end)],...
            'YTickLabel',[],'TickDir','Out','Fontsize',16);
        xlabel('Time (s)'); ylabel('Neuron');
        title('Pro')
        vline(0,'k')
        
        %anti
        t = units(1).anti.neural.sacc.ts_pst; clear r;
        for j=1:length(nunits)
            r(j,:) = units(indx_area(j)).anti.neural.sacc.norm_rate_pst(1,:);
        end
        
        [r,t] = smooth_colormap(r,t);
        % plot colormap
        B = goodcolormap('bwr');
        figure; set(gcf,'Position',[100 200 300 300]);
        %hold on; colormap(B');
        %hold on; colormap(bluewhitered);
        hold on; colormap(bone(4));
        imagesc(t,1:nunits,r,[0,1]);
        set(gca,'xlim',[-0.2 0.255],'ylim',[1 nunits(end)],...
            'YTickLabel',[],'TickDir','Out','Fontsize',16);
        xlabel('Time (s)'); ylabel('Neuron');
        title('Anti')
        vline(0,'k')
        
    case 'waterfall_sacc'  % TODO sort by max?
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        if isempty(recArea)
            disp('no enough input arguments')
        end
        
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        %pro
        t = units(1).pro.neural.sacc.ts_pst; clear r;
        for j=1:length(indx_area)
            r(j,:) = units(indx_area(j)).pro.neural.sacc.norm_rate_pst(1,:);
            indx_neg=find(r(j,:)<0);
            r(j,indx_neg) = 0;
        end
        
        % plot
        figure; hold on; %colormap(d);
        for j=1:length(nunits)
            waterfall(t,j,r(j,:));
        end
        set(gca,'xlim', [-0.255 0.255], 'zlim', [0 2],'CameraViewAngle', 9, 'zTick', [], 'ytick', []);
        view(gca,[-0.399999999999989 59]);
        grid off;
        title(['Pro'])
        
        %anti
        t = units(1).anti.neural.sacc.ts_pst; clear r;
        for j=1:length(indx_area)
            r(j,:) = units(indx_area(j)).anti.neural.sacc.norm_rate_pst(1,:);
            indx_neg=find(r(j,:)<0);
            r(j,indx_neg) = 0;
        end
        
        % plot
        figure; hold on; %colormap(d);
        for j=1:length(nunits)
            waterfall(t,j,r(j,:));
        end
        set(gca,'xlim', [-0.255 0.255], 'zlim', [0 2],'CameraViewAngle', 9, 'zTick', [], 'ytick', []);
        view(gca,[-0.399999999999989 59]);
        grid off;
        title(['Anti'])
        
    case 'scatter_pro_anti'
        %pro
        %gather instr
        % find the min of trials
        pro_trials = size(units(cellNum).pro.neural.instr.nspk,1);anti_trials = size(units(cellNum).anti.neural.instr.nspk,1);
        min_trials = min(pro_trials, anti_trials);
        
        % instr
        r_anti_instr = units(cellNum).anti.neural.instr.nspk(1:min_trials);
        r_pro_instr = units(cellNum).pro.neural.instr.nspk(1:min_trials);
        
        % sacc
        r_anti_sacc = units(cellNum).anti.neural.sacc.nspk(1:min_trials);
        r_pro_sacc = units(cellNum).pro.neural.sacc.nspk(1:min_trials);
        
        %plot
        % instr
        figure; hold on;
        plot(r_pro_instr, r_anti_instr, '.k','MarkerSize', 18);
        xlim([20 110]); ylim([20 110]);
        title('Instruction'); xlabel('Prosaccade'); ylabel('Antisaccade');
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        
        % sacc
        figure; hold on;
        plot(r_pro_sacc, r_anti_sacc, '.k','MarkerSize', 18);
        xlim([20 160]); ylim([20 160]);
        title('Saccade'); xlabel('Prosaccade'); ylabel('Antisaccade');
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        
    case 'scatter_pro_anti_pop'
        % plot avg firing rate pro vs anti for all cells.
        %% instr
        % gather
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        
        for i = 1:length(nunits)
            rate_pro(i) = units(indx_area(i)).pro.neural.instr.rate_mu;
            rate_anti(i) = units(indx_area(i)).anti.neural.instr.rate_mu;
            sig_pro(i) = units(indx_area(i)).pro.neural.instr.rate_sig;
            sig_anti(i) = units(indx_area(i)).anti.neural.instr.rate_sig;
            rate_base_pro(i) = units(indx_area(i)).pro.neural.base.rate_mu;
            sig_base_pro(i) = units(indx_area(i)).pro.neural.base.rate_sig;
            rate_base_anti(i) = units(indx_area(i)).anti.neural.base.rate_mu;
            sig_base_anti(i) = units(indx_area(i)).anti.neural.base.rate_sig;
        end
        % get significantly diff
        for i = 1:length(nunits)
            indx_sign(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr);
        end
        
        %% plot instruction
        figure; hold on;
        errorbar(rate_pro,rate_anti,sig_pro,sig_pro,sig_anti,sig_anti, 'ok','MarkerSize', 4)
        plot(rate_pro(indx_sign),rate_anti(indx_sign), '.c','MarkerSize', 18);
        plot([0 150],[0 150]);
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        title('Instruction'); xlabel('Prosaccade'); ylabel('Antisaccade');
        annotation('textbox',...
            [0.659928571428571 0.152380952380952 0.227571428571429 0.116666666666667],...
            'String',{'n= ' num2str(size(indx_sign,2)),'sign diff = ' num2str(sum(indx_sign))},...
            'FitBoxToText','on');
        axis square;
        
        %% sacc
        for i = 1:length(nunits)
            rate_pro(i) = units(indx_area(i)).pro.neural.sacc.rate_mu;
            rate_anti(i) = units(indx_area(i)).anti.neural.sacc.rate_mu;
            sig_pro(i) = units(indx_area(i)).pro.neural.sacc.rate_sig;
            sig_anti(i) = units(indx_area(i)).anti.neural.sacc.rate_sig;
        end
        % get significantly diff
        clear indx_sign
        for i = 1:length(nunits)
            indx_sign(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        
        % plot
        figure; hold on;
        errorbar(rate_pro,rate_anti,sig_pro,sig_pro,sig_anti,sig_anti, 'ok', 'MarkerSize', 4)
        plot(rate_pro(indx_sign),rate_anti(indx_sign), '.c','MarkerSize', 18);
        plot([0 150],[0 150]);
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        title('Saccade'); xlabel('Prosaccade'); ylabel('Antisaccade');
        annotation('textbox',...
            [0.659928571428571 0.152380952380952 0.227571428571429 0.116666666666667],...
            'String',{'n= ' num2str(size(indx_sign,2)),'sign diff = ' num2str(sum(indx_sign))},...
            'FitBoxToText','on');
        axis square;
        
        %% Plot base vs pro and anti
        % pro
        clear indx_sign
        for i = 1:length(nunits)
            indx_sign(i) = logical(units(indx_area(i)).stats.pro.flags.saccVSbase_nspk);
        end
        
        figure; hold on;
        errorbar(rate_base_pro,rate_pro,sig_base_pro,sig_base_pro,sig_pro,sig_pro, 'ok', 'MarkerSize', 4)
        plot(rate_base_pro(indx_sign),rate_pro(indx_sign), '.k','MarkerSize', 5);
        plot([0 180],[0 180]);
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        title('SaccadeVSBase'); xlabel('Baseline'); ylabel('Prosaccade');
        annotation('textbox',...
            [0.659928571428571 0.152380952380952 0.227571428571429 0.116666666666667],...
            'String',{'n= ' num2str(size(indx_sign,2)),'sign diff = ' num2str(sum(indx_sign))},...
            'FitBoxToText','on');
        axis square;
        
        %anti
        clear indx_sign
        for i = 1:length(nunits)
            indx_sign(i) = logical(units(indx_area(i)).stats.anti.flags.saccVSbase_nspk);
        end
        
        figure; hold on;
        errorbar(rate_base_anti,rate_anti,sig_base_anti,sig_base_anti,sig_anti,sig_anti, 'ok', 'MarkerSize', 4)
        plot(rate_base_anti(indx_sign),rate_anti(indx_sign), '.k','MarkerSize', 18);
        plot([0 180],[0 180]);
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        title('SaccadeVSBase'); xlabel('Baseline'); ylabel('Antisaccade');
        annotation('textbox',...
            [0.659928571428571 0.152380952380952 0.227571428571429 0.116666666666667],...
            'String',{'n= ' num2str(size(indx_sign,2)),'sign diff = ' num2str(sum(indx_sign))},...
            'FitBoxToText','on');
        axis square;
        
        
    case 'peak_resp_instr'
        % gather
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        % indx signif neurons
        for i = 1:length(nunits)
            indx_sign_instr(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr);
        end
        n_instr = sum(indx_sign_instr);
        
        t = units(1).pro.neural.instr.ts_pst_win;
        for i = 1:length(units(indx_area))
            peak_pro(i) = units(indx_area(i)).pro.neural.instr.peak_resp;
            peak_anti(i) = units(indx_area(i)).anti.neural.instr.peak_resp;
            peak_time_pro(i) = units(indx_area(i)).pro.neural.instr.peak_resp_time;
            peak_time_anti(i) = units(indx_area(i)).anti.neural.instr.peak_resp_time;
        end
        % with baseline subtracted
        for i = 1:length(units(indx_area))
            peak_pro_diff(i) = units(indx_area(i)).pro.neural.instr.peak_resp-units(i).pro.neural.base.rate_mu;
            peak_anti_diff(i) = units(indx_area(i)).anti.neural.instr.peak_resp-units(i).anti.neural.base.rate_mu;
        end
        
        
        % plot (baseline subtracted)
        figure; hold on;
        plot(peak_pro_diff,peak_anti_diff, '.k', 'MarkerSize', 16);
        plot(peak_pro_diff(indx_sign_instr), peak_anti_diff(indx_sign_instr), '.c', 'MarkerSize', 16);
        %set(gca,'XScale','Log','YScale','Log' ,'FontSize', 18, 'TickDir', 'out');axis ([1e0 1e2 1e0 1e2]);
        set(gca,'FontSize', 18, 'TickDir', 'out');
        plot([-60 10],[-60 120]);
        title('Instr peak resp (baseline subtracted)'); xlabel('Prosaccade'); ylabel('Antisaccade');
        box off; axis square
        [h,p] = ttest(peak_pro_diff,peak_anti_diff)
        % plot peak time resp
        h_pro = hist(peak_time_pro,t);
        h_anti = hist(peak_time_anti,t);
        figure; hold on;
        plot(t,h_pro, 'r', 'LineWidth', 2);
        plot(t,h_anti, 'g', 'LineWidth', 2);
        set(gca, 'TickDir', 'out'); box off;
        title('Peak resp time'); xlabel('Time (s)'); ylabel('Neuron nr')
        
    case 'peak_resp_sacc'
        % gather
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        % indx signif neurons
        for i = 1:length(nunits)
            indx_sign_sacc(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        n_sacc = sum(indx_sign_sacc);
        
        t = units(1).pro.neural.sacc.ts_pst_win;
        for i = 1:length(units(indx_area))
            peak_pro(i) = units(indx_area(i)).pro.neural.sacc.peak_resp;
            peak_anti(i) = units(indx_area(i)).anti.neural.sacc.peak_resp;
            peak_time_pro(i) = units(indx_area(i)).pro.neural.sacc.peak_resp_time;
            peak_time_anti(i) = units(indx_area(i)).anti.neural.sacc.peak_resp_time;
        end
        % with baseline subtracted
        for i = 1:length(units(indx_area))
            peak_pro_diff(i) = units(indx_area(i)).pro.neural.sacc.peak_resp-units(i).pro.neural.base.rate_mu;
            peak_anti_diff(i) = units(indx_area(i)).anti.neural.sacc.peak_resp-units(i).anti.neural.base.rate_mu;
        end
        
        
        % plot (baseline subtracted)
        figure; hold on;
        plot(peak_pro_diff,peak_anti_diff, '.k', 'MarkerSize', 16);
        plot(peak_pro_diff(indx_sign_sacc), peak_anti_diff(indx_sign_sacc), '.c', 'MarkerSize', 16);
        %set(gca,'XScale','Log','YScale','Log' ,'FontSize', 18, 'TickDir', 'out');axis ([1e0 1e2 1e0 1e2]);
        set(gca,'FontSize', 18, 'TickDir', 'out');
        plot([-150 100],[-150 100]);
        title('Sacc peak resp (baseline subtracted)'); xlabel('Prosaccade'); ylabel('Antisaccade');
        box off; axis square
        [h,p] = ttest(abs(peak_pro_diff),abs(peak_anti_diff))
        % plot peak time resp
        h_pro = hist(peak_time_pro,t);
        h_anti = hist(peak_time_anti,t);
        figure; hold on;
        plot(t,h_pro, 'r', 'LineWidth', 2);
        plot(t,h_anti, 'g', 'LineWidth', 2);
        set(gca, 'TickDir', 'out'); box off;
        title('Peak resp time'); xlabel('Time (s)'); ylabel('Neuron nr')
        
    case 'firingVSamp'
        %gather
        amp_pro = [units(cellNum).pro.behav.trial.saccAmplitude];
        amp_anti = [units(cellNum).anti.behav.trial.saccAmplitude];
        
        spk_pro = [units(cellNum).pro.neural.sacc.nspk];
        spk_anti = [units(cellNum).anti.neural.sacc.nspk];
        
        % correlation
        [rho_pro,pval_pro] = corr(amp_pro',spk_pro)
        [rho_anti,pval_anti] = corr(amp_anti',spk_anti)
        
        %plot
        figure; hold on; box off
        plot(amp_pro,spk_pro, '.r', 'MarkerSize',16);
        plot(amp_anti,spk_anti, '.g','MarkerSize',16);
        xlabel('Amplitude (deg)'); ylabel('Firing rate (spk/s)');
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        lsline(gca)
        
        
        %Avg per amplitude blocks
        block1 = 4.9;
        block2 = 9.9;
        block3 = 15;
        
        %pro
        block1_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]<=block1);
        block2_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]>block1 & [units(cellNum).pro.behav.trial.saccAmplitude]<= block2);
        block3_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]>block2 & [units(cellNum).pro.behav.trial.saccAmplitude]<= block3);
        
        for i = 1:length(block1_pro_indx)
            block1_pro(i,1) = units(cellNum).pro.neural.sacc.nspk(block1_pro_indx(i));
        end
        for i = 1:length(block2_pro_indx)
            block2_pro(i,1) = units(cellNum).pro.neural.sacc.nspk(block2_pro_indx(i));
        end
        for i = 1:length(block3_pro_indx)
            block3_pro(i,1) = units(cellNum).pro.neural.sacc.nspk(block1_pro_indx(i));
        end
        
        
        %anti
        block1_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]<=block1);
        block2_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]>block1 & [units(cellNum).anti.behav.trial.saccAmplitude]<= block2);
        block3_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]>block2 & [units(cellNum).anti.behav.trial.saccAmplitude]<= block3);
        
        if ~isempty(block1_anti_indx)
            for i = 1:length(block1_anti_indx)
                block1_anti(i,1) = units(cellNum).anti.neural.sacc.nspk(block1_anti_indx(i));
            end
        else
            block1_anti(i,1) = NaN;
        end
        for i = 1:length(block2_anti_indx)
            block2_anti(i,1) = units(cellNum).anti.neural.sacc.nspk(block2_anti_indx(i));
        end
        for i = 1:length(block3_anti_indx)
            block3_anti(i,1) = units(cellNum).anti.neural.sacc.nspk(block3_anti_indx(i));
        end
        
        % mean pro
        %         block1_pro_mu = nanmean(block1_pro); block1_pro_sig = nanstd(block1_pro);
        %         block2_pro_mu = nanmean(block2_pro); block2_pro_sig = nanstd(block2_pro);
        %         block3_pro_mu = nanmean(block3_pro); block1_pro_sig = nanstd(block3_pro);
        %
        %         % mean anti
        %         block1_anti_mu = nanmean(block1_anti); block1_anti_sig = nanstd(block1_anti);
        %         block2_anti_mu = nanmean(block2_anti); block2_anti_sig = nanstd(block2_anti);
        %         block3_anti_mu = nanmean(block3_anti); block1_anti_sig = nanstd(block3_anti);
        %
        %         %gather
        %         pro_blocks = [block1_pro_mu block2_pro_mu block3_pro_mu ; block1_pro_sig block2_pro_sig block1_pro_sig];
        %         anti_blocks = [block1_anti_mu block2_anti_mu block3_anti_mu ; block1_anti_sig block2_anti_sig block1_anti_sig];
        %
        %         % plot means
        %         figure; hold on; box off
        %         errorbar(pro_blocks(1,:),pro_blocks(2,:), 'Color', 'r','LineWidth', 2);
        %         errorbar(anti_blocks(1,:),anti_blocks(2,:), 'Color', 'g','LineWidth', 2);
        %         set(gca, 'TickDir', 'out', 'xlim',[0.5 3.5]);
        
        % plot rho for all cells
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        for i = 1:length(units(indx_area))
            amp_pro_all{i,:} = [units(indx_area(i)).pro.behav.trial.saccAmplitude];
            amp_anti_all{i,:} = [units(indx_area(i)).anti.behav.trial.saccAmplitude];
            
            spk_pro_all{i,:} = [units(indx_area(i)).pro.neural.sacc.nspk];
            spk_anti_all{i,:} = [units(indx_area(i)).anti.neural.sacc.nspk];
            
            [rho_pro_all(i,:),pval_pro_all(i,:)] = corr(amp_pro_all{i,:}',spk_pro_all{i,:});
            [rho_anti_all(i,:),pval_anti_all(i,:)] = corr(amp_anti_all{i,:}',spk_anti_all{i,:});
        end
        
        z=1;
        
    case 'firingVSvel'
        pv_pro = [units(cellNum).pro.behav.trial.saccPeakVel];
        pv_anti = [units(cellNum).anti.behav.trial.saccPeakVel];
        
        spk_pro = [units(cellNum).pro.neural.sacc.nspk];
        spk_anti = [units(cellNum).anti.neural.sacc.nspk];
        
        %plot
        figure; hold on; box off
        plot(pv_pro,spk_pro, '.r', 'MarkerSize',16);
        plot(pv_anti,spk_anti, '.g','MarkerSize',16);
        xlabel('Peak Vel (deg/s)'); ylabel('Firing rate (spk/s)');
        set(gca, 'TickDir', 'out', 'FontSize', 18)
        lsline(gca)
        
        % plot rho for all cells
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        for i = 1:length(units(indx_area))
            vel_pro_all{i,:} = [units(indx_area(i)).pro.behav.trial.saccPeakVel];
            vel_anti_all{i,:} = [units(indx_area(i)).anti.behav.trial.saccPeakVel];
            
            spk_pro_all{i,:} = [units(indx_area(i)).pro.neural.sacc.nspk];
            spk_anti_all{i,:} = [units(indx_area(i)).anti.neural.sacc.nspk];
            
            [rho_pro_all(i,:),pval_pro_all(i,:)] = corr(vel_pro_all{i,:}',spk_pro_all{i,:});
            [rho_anti_all(i,:),pval_anti_all(i,:)] = corr(vel_anti_all{i,:}',spk_anti_all{i,:});
        end
        
        z=1;
        
    case 'firingVSdur'
        dur_pro = [units(cellNum).pro.behav.trial.saccDuration];
        dur_anti = [units(cellNum).anti.behav.trial.saccDuration];
        
        spk_pro = [units(cellNum).pro.neural.sacc.nspk];
        spk_anti = [units(cellNum).anti.neural.sacc.nspk];
        
        %plot
        figure; hold on; box off
        plot(dur_pro,spk_pro, '.r', 'MarkerSize',16);
        plot(dur_anti,spk_anti, '.g','MarkerSize',16);
        xlabel('Duration (ms)'); ylabel('Firing rate (spk/s)');
        set(gca, 'TickDir', 'out', 'FontSize', 18)
        lsline(gca)
        
        % plot rho for all cells
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        
        for i = 1:length(units(indx_area))
            dur_pro_all{i,:} = [units(indx_area(i)).pro.behav.trial.saccDuration];
            dur_anti_all{i,:} = [units(indx_area(i)).anti.behav.trial.saccDuration];
            
            spk_pro_all{i,:} = [units(indx_area(i)).pro.neural.sacc.nspk];
            spk_anti_all{i,:} = [units(indx_area(i)).anti.neural.sacc.nspk];
            
            [rho_pro_all(i,:),pval_pro_all(i,:)] = corr(dur_pro_all{i,:}',spk_pro_all{i,:});
            [rho_anti_all(i,:),pval_anti_all(i,:)] = corr(dur_anti_all{i,:}',spk_anti_all{i,:});
        end
        
        z=1;
        
        
    case 'binomial_pb_dist'
        % statistical test to compare if pro == anti - aligned to instr/sacc using spk count
        % Plot Z-statistic over time.
        % instr
        % gather
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        %% Uncomment to plot per neuron
        for i=1:length(indx_area)
            %
            stat_instr(i,:) = units(indx_area(i)).stats.instr.pval.pbDist_testStat;
            t_instr = units(indx_area(i)).pro.neural.instr.ts_pst;
            position_stat_sign_instr(i,:) = abs(stat_instr(i,:))>=1.96;
            %
            % plot instr
%             figure; subplot(2,1,1);
%             plot(t_instr, stat_instr(i,:), 'k','MarkerSize', 15);
%             hline(-1.96, 'k');hline(1.96, 'k');vline(0, 'c');
%             set(gca, 'xlim',[-0.150 0.151],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
%             title ('Binomial pb dist - Instruction')
%             xlabel('time (s)')
%             %
%             %             % sacc
%             %             %gather
            t_sacc = units(indx_area(i)).pro.neural.sacc.ts_pst;
            stat_sacc(i,:) = units(indx_area(i)).stats.sacc.pval.pbDist_testStat;
            position_stat_sign_sacc(i,:) = abs(stat_sacc(i,:))>=1.96;
%             %             % plot
%             subplot(2,1,2)
%             plot(t_sacc, stat_sacc(i,:), 'k','MarkerSize', 15);
%             hline(-1.96, 'k');hline(1.96, 'k');vline(0, 'c');
%             set(gca, 'xlim',[-0.150 0.151],'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
%             title (['Binomial pb dist - Saccade cell: ' num2str(indx_area(i))])
%             xlabel('time')
%             fname = 'Binomial_pb_dist';
%             print(fname,'-append', '-dpsc2')
%             waitforbuttonpress; close all;
            %
        end
        
        %% Disused %% Take the absolute value of the Z-statistic and average across neurons - average for each time point and plot it as a function of time
        % This will reveal how well an average neuron can discriminate between the two conditions.
        
        %         sacc_Z_stat = mean(abs(stat_sacc));
        %         instr_Z_stat = mean(abs(stat_instr));
        %
        %         % plot instr Z stat
        %         figure; hold on;
        %         plot(abs(stat_instr), '.k'); hold on; hline(1.96);  % data points 11 to 41
        %         %plot(nanmedian(abs(stat_instr)), '-r', 'LineWidth',2)
        %         set(gca, 'xlim',[10.5 41.5], 'TickDir', 'out', 'FontSize', 18);
        %
        %         % side histogram
        %         figure; hold on;
        %         histogram(abs(stat_instr), 25)
        %         set(gca,'Xdir','reverse','TickDir', 'out', 'FontSize', 18);
        %
        %         % Histogram of >1.96 only
        %         figure; hold on
        %         histogram(stat_instr(position_stat_sign_instr),25);
        %         plot(t_sacc(),abs(stat_instr(position_stat_sign_instr)), '.k')
        %
        %         figure; hold on
        %         histogram(abs(stat_instr(position_stat_sign_instr)),25);
        %          set(gca,'TickDir', 'out', 'FontSize', 18);
        %
        %          % plot sacc Z stat
        %          figure; hold on;
        %         plot(abs(stat_sacc), '.k'); hold on; hline(1.96);  % data points 11 to 41
        %         %plot(nanmedian(abs(stat_instr)), '-r', 'LineWidth',2)
        %         set(gca, 'xlim',[10.5 41.5], 'TickDir', 'out', 'FontSize', 18);
        %
        %         % side histogram
        %         figure; hold on;
        %         histogram(abs(stat_sacc), 25)
        %         set(gca,'Xdir','reverse','TickDir', 'out', 'FontSize', 18);
        
        %% Take the absolute value of the Z-statistic and average across significantly different neurons
        %                    instr -> get significantly different neurons
        % get significantly diff instr
        nunits = 1:length(indx_area);
        for i = 1:length(nunits)
            indx_sign_instr(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr);
        end
        n_instr = sum(indx_sign_instr);
        
        % plot
        z_sign_instr = stat_instr(indx_sign_instr,:);
        figure; hold on;
        plot(t_instr,smooth(nanmean(abs(z_sign_instr)),3),'LineWidth', 2,'Color','k'); % smoothed
        hline(1.96,'k')
        set(gca, 'xlim',[0 0.3],'ylim',[0 2], 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Instr => ' recArea ' n= ' num2str(n_instr)]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        % histogram of z stat instr
        figure; hold on;
        histogram(z_sign_instr,50);
        vline([-1.96 1.96])
        title ('Z stat instruction signif neurons')
        
        % get Z stat for those neurons and plot
        %plot
        figure; hold on;
        plot(t_instr,abs(z_sign_instr), '.k');
        set(gca, 'xlim',[-0.1 0.3], 'TickDir', 'out', 'FontSize', 18);
        hline(1.96);xlabel('time(s)'); ylabel('Z-stat')
        
        plot(t_instr,smooth(nanmean(abs(z_sign_instr)),3),'LineWidth', 2,'Color','k'); % smoothed
        hline(1.96,'k')
        set(gca, 'xlim',[0 0.3],'ylim',[0 2], 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Instr => ' recArea ' n= ' num2str(n_instr)]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        %   sacc -> get significantly different neurons
        for i = 1:length(nunits)
            indx_sign_sacc(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        n_sacc = sum(indx_sign_sacc);
        
        % get Z stat for those neurons and plot
        z_sign_sacc = stat_sacc(indx_sign_sacc,:);
        figure; hold on;
        plot(t_sacc,smooth(nanmean(abs(z_sign_sacc)),3),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[-0.150 0.151], 'TickDir', 'out', 'FontSize', 18);
        hline(1.96,'--k')
        title(['Abs Z stat Sacc => ' recArea ' n= ' num2str(n_sacc)]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        % histogram of z stat instr
        figure; hold on;
        histogram(z_sign_sacc,100);
        vline([-1.96 1.96])
        title ('Z stat sacc signif neurons')
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        
        %plot(t_sacc,nanmean(abs(z_sign_sacc)));
        figure; hold on;
        plot(t_sacc,abs(z_sign_sacc), '.k');
        set(gca, 'xlim',[-0.1 0.3], 'TickDir', 'out', 'FontSize', 18);
        hline(1.96);
        xlabel('time(s)'); ylabel('Z-stat')
        
        plot(t_sacc,smooth(nanmean(abs(z_sign_sacc)),3),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[-0.150 0.151], 'ylim', ([0.6 2.2]), 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Sacc => ' recArea ' n= ' num2str(n_sacc)]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        %% Same as above but take mean of Z-stat only on sacc window (0-0.2s)
        sacc_win = t_sacc>-0.100 & t_sacc<0.21;
        t_sacc_win = t_sacc(sacc_win);
        z_sign_sacc_win = z_sign_sacc(:,sacc_win);
        
        for i=1:size(z_sign_sacc_win,1)
            z_win_mu(i) = mean(abs(z_sign_sacc_win(i,:))); % mean for every cell in its window
        end
        
        figure; hold on;
        plot(z_win_mu,'.k', 'MarkerSize', 18)
        set(gca, 'TickDir', 'out', 'FontSize', 18); hline(1.96)
        xlabel('cell'); ylabel('Z-stat')
        title(['Z stat sacc window(0.1-0.2s) recarea => ' recArea])
        
        
    case 'spk_pb_stat'
        % instr
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        for i=1:length(indx_area)
            % gather
            stat_instr= units(indx_area(i)).stats.instr.pval.spk_count_bigWin;
            t_instr = units(indx_area(i)).pro.neural.instr.ts_spkCount_win(:,1)';
            signif_instr(i,:) = stat_instr<=0.05;
        end
        
        signif_instr = sum(signif_instr);
        %plot
        figure; hold on;
        plot(t_instr, smooth(signif_instr,5),'Linewidth',2);
        set(gca, 'xlim',[-0.05 0.5], 'TickDir', 'out', 'FontSize', 18);  vline(0);
        title ('spk pb stat - Instruction')
        xlabel('Time (s)'); ylabel('P <0.05 counts')
        
        figure; hold on;
        bar(t_instr,signif_instr);
        set(gca, 'xlim',[-0.05 0.5], 'TickDir', 'out', 'FontSize', 18);  vline(0);
        title ('spk pb stat - Instruction')
        xlabel('Time (s)'); ylabel('P <0.05 counts')
        
        % sacc
        %gather
        for i=1:length(indx_area)
            % gather
            stat_sacc= units(indx_area(i)).stats.sacc.pval.spk_count_bigWin;
            t_sacc = units(indx_area(i)).pro.neural.sacc.ts_spkCount_win(:,1)';
            signif_sacc(i,:) = stat_sacc<=0.05;
        end
        
        signif_sacc = sum(signif_sacc);
        %plot
        figure; hold on;
        plot(t_sacc, smooth(signif_sacc,5),'Linewidth',2);
        set(gca, 'xlim',[-0.150 0.250], 'TickDir', 'out', 'FontSize', 18);  vline(0);
        title ('spk pb stat - Sacc')
        xlabel('Time (s)'); ylabel('P <0.05 counts')
        
        figure; hold on;
        bar(t_sacc,signif_sacc);
        set(gca, 'xlim',[-0.05 0.250], 'TickDir', 'out', 'FontSize', 18);  vline(0);
        title ('spk pb stat - Sacc')
        xlabel('Time (s)'); ylabel('P <0.05 counts')
        
        % only for signif diff neurons sacc 
        for i = 1:length(indx_area)
            indx_sign_sacc(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        unit = units(indx_sign_sacc); 
        
        for i=1:length(unit)
            % gather
            stat_sacc_diff= unit(i).stats.sacc.pval.spk_count_bigWin;
            t_sacc = unit(i).pro.neural.sacc.ts_spkCount_win(:,1)';
            signif_sacc_diff(i,:) = stat_sacc_diff<=0.05;
        end
        signif_sacc_diff = sum(signif_sacc_diff);
        
        figure; hold on;
        plot(t_sacc, smooth(signif_sacc_diff,5),'Linewidth',2);
        set(gca, 'xlim',[-0.150 0.250], 'TickDir', 'out', 'FontSize', 18);  vline(0);
        title ('spk pb stat - Sacc diff neurons')
        xlabel('Time (s)'); ylabel('P <0.05 counts')
        
        
    case 'rosePlots'
        
        % plot mean firing rate + STD of different directions for pro and
        % anti during sac or instruction period
        
        % TODO - Correct for different amplitudes (maybe in spikes/degree)
        
        
        proLegend =[10 8 15 14 4 11 5 1; 1 2 3 4 5 6 7 8] ;
        antiLegend = [2 16 7 6 12 3 13 9; 1 2 3 4 5 6 7 8];
        
        sacc= 01;
        instr=0;
        
        if sacc
            pro_trials =units(cellNum).pro.neural.sacc.nspk ;
            anti_trials =units(cellNum).pro.neural.sacc.nspk ;
            title_flag = 'saccade'            ;
        elseif instr
            pro_trials =units(cellNum).pro.neural.instr.nspk ;
            anti_trials =units(cellNum).pro.neural.instr.nspk ;
            title_flag = 'instruction';
            
        end
        
        % loop over all directions to gather average and std's of firing
        % freqeuncies
        for i=1:8
            
            thisDirection = proLegend(1,i);
            trialIndex = find([units(cellNum).pro.behav.trial.conditionCode]==thisDirection);
            meanRates(1,i) = mean(pro_trials(trialIndex));
            upperbound(1,i) = meanRates(1,i) + std(pro_trials(trialIndex));
            lowerbound(1,i) = meanRates(1,i) - std(pro_trials(trialIndex));
            
            thisDirection = antiLegend(1,i);
            trialIndex = find([units(cellNum).anti.behav.trial.conditionCode]==thisDirection);
            meanRates(2,i) = mean(anti_trials(trialIndex));
            upperbound(2,i) = meanRates(2,i) + std(anti_trials(trialIndex));
            lowerbound(2,i) = meanRates(2,i) - std(anti_trials(trialIndex));
            
        end
        
        % do the actual plotting
        rosePlotProAnti(meanRates,lowerbound,upperbound  )
        
        set(gca,   'TickDir', 'out', 'FontSize', 18)
        title (['Average rates during ' title_flag ' window'])
        xlabel('Firing frequency(hz)')
        ylabel('Firing frequency(hz)')
        
    case 'bi_resp'
        
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        for i=1:length(indx_area)
            % gather
            indx_sign_sacc(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
            indx_sign_instr(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr);
        end
        
        signif_sacc = sum(indx_sign_sacc);
        signif_instr = sum(indx_sign_instr);
        
        signif_both = indx_sign_sacc & indx_sign_instr;
        num_both = sum(signif_both)
        
        
    case 'cv'
        
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        for i=1:length(indx_area)
            % gather
            cv_pro(i)= units(i).pro.neural.sacc.rate_std/units(i).pro.neural.sacc.rate_mu;
            cv_anti(i) = units(i).anti.neural.sacc.rate_std/units(i).anti.neural.sacc.rate_mu;
            indx_sign_sacc(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        
        %         histogram(cv_pro(indx_sign_sacc),20); hold on;
        %         histogram(cv_anti(indx_sign_sacc),20);
        
        figure; hold on;
        plot(cv_pro, cv_anti, '.k', 'MarkerSize', 18);
        plot(cv_pro(indx_sign_sacc), cv_anti(indx_sign_sacc), '.c', 'MarkerSize', 18);
        xlim([0 0.6]); ylim([0 0.6]); plot([0:0.1:1],[0:0.1:1])
        set (gca, 'TickDir', 'out','FontSize', 18); box off;
        xlabel('CV Pro');ylabel('CV Anti');
        title(recArea); axis square
        
        cv_pro_mu = mean(cv_pro(indx_sign_sacc)); cv_pro_sig = std(cv_pro(indx_sign_sacc))/sqrt(sum(indx_sign_sacc));
        cv_anti_mu = mean(cv_anti(indx_sign_sacc)); cv_anti_sig = std(cv_anti(indx_sign_sacc))/sqrt(sum(indx_sign_sacc));
        cv_all = [cv_pro_mu cv_anti_mu]; cv_std_all = [cv_pro_sig cv_anti_sig];
        
        axes('Position',[.7 .2 .2 .2]); hold on;
        bar(cv_all); errorbar(cv_all,cv_std_all);
        set (gca, 'TickDir', 'out','FontSize', 12, 'xlim',[0.5 2.5], 'xTick', [], 'xTickLabel', [], 'ylim', [0 0.4]); box off;
        plot(1,cv_pro(indx_sign_sacc),'.k', 'MarkerSize',14); plot(2,cv_anti(indx_sign_sacc),'.k', 'MarkerSize',14);
        title('CV signif cells')
        
        [h,p] = ttest(cv_pro(indx_sign_sacc),cv_anti(indx_sign_sacc))
        [h,p] = ttest(cv_pro,cv_anti)
        
    case 'cv2'
        
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        for i=1:length(indx_area)
            % gather
            cv2_pro(i)= mean([units(i).pro.neural.trial.cv2]);
            cv2_anti(i) = mean([units(i).anti.neural.trial.cv2]);
            indx_sign_sacc(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        
        figure; hold on;
        plot(cv2_pro, cv2_anti, '.k', 'MarkerSize', 18);
        plot(cv2_pro(indx_sign_sacc), cv2_anti(indx_sign_sacc), '.c', 'MarkerSize', 18);
        xlim([0.3 0.6]); ylim([0.3 0.6]); plot([0:0.1:1],[0:0.1:1])
        set (gca, 'TickDir', 'out','FontSize', 18, 'xTick', [0.3 0.4 0.5 0.6], 'yTick', [0.3 0.4 0.5 0.6]); box off;
        xlabel('CV2 Pro');ylabel('CV2 Anti');
        title(recArea); axis square
        
        cv2_pro_mu = mean(cv2_pro(indx_sign_sacc)); cv2_pro_sig = std(cv2_pro(indx_sign_sacc))/sqrt(sum(indx_sign_sacc));
        cv2_anti_mu = mean(cv2_anti(indx_sign_sacc)); cv2_anti_sig = std(cv2_anti(indx_sign_sacc))/sqrt(sum(indx_sign_sacc));
        cv2_all = [cv2_pro_mu cv2_anti_mu]; cv2_std_all = [cv2_pro_sig cv2_anti_sig];
        
        axes('Position',[.7 .2 .2 .2]); hold on;
        bar(cv2_all); errorbar(cv2_all,cv2_std_all);
        set (gca, 'TickDir', 'out','FontSize', 12, 'xlim',[0.5 2.5], 'xTick', [1, 2], 'xTickLabel', ['Pro' ,'Anti'], 'ylim', [0 0.7]); box off;
        plot(1,cv2_pro(indx_sign_sacc),'.k', 'MarkerSize',14); plot(2,cv2_anti(indx_sign_sacc),'.k', 'MarkerSize',14);
        title('CV2 signif cells')
        
        % [h,p] = ttest(cv2_pro(indx_sign_sacc),cv2_anti(indx_sign_sacc))
        [h,p] = ttest(cv2_pro,cv2_anti)
        
        
    case 'eye_move'
      % LOAD FILE SENT BY NICO  
        
        t = units(cellNum).pro.timepointsEye/1000; % time
        
        % extract eye trace around saccade onset
        % pro
        c = 1;
        for i = 1:length(units(cellNum).pro.behav.trial)
            sacc_onset = units(cellNum).pro.behav.trial(i).saccadeOnset;
            align_sacc = t - sacc_onset; t_sacc = align_sacc(align_sacc >= -0.150 & align_sacc <= 0.150);
            x = units(cellNum).pro.behav.trial(i).eyePositionX; y = units(cellNum).pro.behav.trial(i).eyePositionY;
            if ~isnan(x)
                heye(c,:) = x(align_sacc >= -0.150 & align_sacc <= 0.150);
                veye(c,:) = y(align_sacc >= -0.150 & align_sacc <= 0.150);
                c = c+1;
            end
        end
        
        figure; hold on;
        plot(t_sacc,heye','b', 'LineWidth', 2);
        plot(t_sacc,veye','c', 'LineWidth', 2);
        set(gca, 'yTick',[], 'TickDir', 'out', 'FontSize',18);
        title(['Eye trace pro cell ' num2str(cellNum)])
        
        % anti
        c = 1;
        for i = 1:length(units(cellNum).anti.behav.trial)
            sacc_onset = units(cellNum).anti.behav.trial(i).saccadeOnset;
            align_sacc = t - sacc_onset; t_sacc = align_sacc(align_sacc >= -0.150 & align_sacc <= 0.150);
            x = units(cellNum).anti.behav.trial(i).eyePositionX; y = units(cellNum).anti.behav.trial(i).eyePositionY;
            if ~isnan(x)
                heye(c,:) = x(align_sacc >= -0.150 & align_sacc <= 0.150);
                veye(c,:) = y(align_sacc >= -0.150 & align_sacc <= 0.150);
                c = c+1;
            end
        end
        
        figure; hold on;
        plot(t_sacc,heye','b', 'LineWidth', 2);
        plot(t_sacc,veye','c', 'LineWidth', 2);
        set(gca, 'yTick',[], 'TickDir', 'out', 'FontSize',18);
        title(['Eye trace anti cell ' num2str(cellNum)])
        
        % print('eye_trace','-depsc2', '-painters', '-cmyk')
        
    case 'peak_vel_distr'
        % get peak vel
        for i = 1:length(units(cellNum).pro.behav.trial), peak_pro(i,:) = units(cellNum).pro.behav.trial(i).saccPeakVel; end
        for i = 1:length(units(cellNum).anti.behav.trial), peak_anti(i,:) = units(cellNum).anti.behav.trial(i).saccPeakVel; end
        
        histogram(peak_pro,30); hold on
        histogram(peak_anti,30);
        
    case 'rate_sorted_vel'
        
        for i = 1:length(units(cellNum).pro.behav.trial), peak_pro(i,:) = units(cellNum).pro.behav.trial(i).saccPeakVel; end
        for i = 1:length(units(cellNum).anti.behav.trial), peak_anti(i,:) = units(cellNum).anti.behav.trial(i).saccPeakVel; end
        for j=1:length(units(cellNum).pro.behav.trial)
            r_pro(j,:) = units(indx(j)).pro.neural.sacc.rate_pst_win; % psth
            
        end
         r_anti(j,:) = units(indx(j)).anti.neural.sacc.rate_pst_win; % psth
         
    case 'mod_ratio'
        % gather indx vermis and lateral
        for cellNum = 1:length(units)
            indx_area_vermis(cellNum) = strcmp(units(cellNum).area, 'vermis');
            indx_area_lat(cellNum) = strcmp(units(cellNum).area, 'lateral');
        end
        indx_area_vermis = find(indx_area_vermis); unit_vermis = units(indx_area_vermis); 
        indx_area_lat = find(indx_area_lat); unit_lat = units(indx_area_lat); 
        
        % for all neurons
        for i = 1:length(unit_vermis)
        mod_r_vermis(i) = unit_vermis(i).stats.mod_ratio; 
        end
        
        for i = 1:length(unit_lat)
        mod_r_lat(i) = unit_lat(i).stats.mod_ratio;
        end
        mod_r_lat(12)=[]; %outlier. Check what's going on
        median_vermis = median(mod_r_vermis); median_lat = median(mod_r_lat);
        % plot
        figure; hold on; 
        cdfplot(mod_r_vermis);
        cdfplot(mod_r_lat);
        vline(1, '-k');
        set(gca, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', 'FontSize', 20); title('Modulation ratio OMV vs Lateral')
        
        
        % for significantly different
        
        for i = 1:length(indx_area_vermis)
            indx_sign_vermis(i) = logical(units(indx_area_vermis(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        unit_vermis_sign = unit_vermis(indx_sign_vermis);
        
        for i = 1:length(indx_area_lat)
            indx_sign_lat(i) = logical(units(indx_area_lat(i)).stats.sacc.flags.proVsAnti_sacc);
        end
        unit_lat_sign = unit_lat(indx_sign_lat);
        
        for i = 1:length(unit_vermis_sign)
        mod_r_vermis_sign(i) = unit_vermis_sign(i).stats.mod_ratio; 
        end
        
        for i = 1:length(unit_lat_sign)
        mod_r_lat_sign(i) = unit_lat_sign(i).stats.mod_ratio;
        end
        mod_r_lat_sign(7)=[];
        
        figure; hold on; 
        cdfplot(mod_r_vermis_sign);
        cdfplot(mod_r_lat_sign);
        vline(1);
        set(gca, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', 'FontSize', 20); title('Modulation ratio OMV vs Lateral signif')
        
        
        
end
        
        







        
end
