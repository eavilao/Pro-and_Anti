function plotUnit_ProAnti(units, plotType, cellNum, recArea)

% Needed: units.mat
% Input:    units  - output from extractWholeNeuronResults.m
%           cellNum - cell you want to plot. If all leave empty = []
%           plotType - plot you want (e.g. raster, psth, etc.)
%           recArea - OMV, lateral Cb

% if running population (all neurons), just input 0 in cellNum


%%
% shadedErrorBar(units(1).anti.neural.sacc_ts_pst,units(1).anti.neural.sacc_rate_pst,repmat(units(1).anti.neural.sacc_rate_sig,1,21))
% hold on
% shadedErrorBar(units(1).pro.neural.sacc_ts_pst,units(1).pro.neural.sacc_rate_pst,repmat(units(1).pro.neural.sacc_rate_sig,1,21))

%%

switch plotType
    
    case 'raster_sacc'
        % saccade aligned
        % pro
        [~,indx] = sort([units(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).pro.behav.trial(indx).reactionTime];
        r_pro= units(cellNum).pro.neural.trial;
        
        figure;subplot (2,1,1); hold on;box off
        
        if strcmp(units(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_SS_align_sacc)
                    plot(sorted_RT(j),j,'.r');
                    %plot(r_pro(indx(j)).tspk_SS_align_sacc(1:2:end),j,'.k'); %plot every n spikes
                    plot(r_pro(indx(j)).tspk_SS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.3 0.3]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Prosaccade (aligned to saccade onset) SS');xlabel('Time (s)');ylabel('Trial Num')
        else
            for j= 1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_CS_align_sacc)
                    plot(sorted_RT(j),j,'.r');
                    plot(r_pro(indx(j)).tspk_CS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.3 0.3]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Prosaccade (aligned to saccade onset) CS');xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
        % anti
        [~,indx] = sort([units(cellNum).anti.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).anti.behav.trial(indx).reactionTime];
        r_anti=units(cellNum).anti.neural.trial;
        subplot (2,1,2); hold on;box off
        
        if strcmp(units(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_SS_align_sacc)
                    plot(sorted_RT(j),j,'.r');
                    %plot(r_anti(indx(j)).tspk_SS_align_sacc(1:2:end),j,'.k'); %plot every n spikes
                    plot(r_anti(indx(j)).tspk_SS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.3 0.3]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Antisaccade (aligned to saccade onset) SS');xlabel('Time (s)');ylabel('Trial Num')
        else
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_CS_align_sacc)
                    plot(sorted_RT(j),j,'.r');
                    plot(r_anti(indx(j)).tspk_CS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.3 0.3]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Antisaccade (aligned to saccade onset) CS');xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
    case 'raster_instr'
        % instr aligned
        % pro
        [~,indx] = sort([units(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).pro.behav.trial(indx).reactionTime];
        r_pro= units(cellNum).pro.neural.trial;
        
        figure;subplot (2,1,1); hold on;box off
        
        if strcmp(units(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_SS)
                    plot(r_pro(indx(j)).tspk_SS(1:2:end),j,'.k'); %plot every n spikes
                    %plot(r_pro(indx(j)).tspk_SS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Prosaccade (aligned to instruction onset) SS');xlabel('Time (s)');ylabel('Trial Num')
        else
            for j=1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_CS)
                    plot(r_pro(indx(j)).tspk_CS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Prosaccade (aligned to instruction onset) CS');xlabel('Time (s)');ylabel('Trial Num')
            
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
            set (gca, 'xlim', ([-0.1 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Antisaccade (aligned to instruction onset) SS');xlabel('Time (s)');ylabel('Trial Num')
        else
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_CS)
                    plot(r_anti(indx(j)).tspk_CS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Antisaccade (aligned to saccade onset) CS');xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
        
    case 'psth'
        
        %gather
        t= units(cellNum).pro.neural.sacc.ts_pst;
        r_pro= units(cellNum).pro.neural.sacc.rate_pst;
        sem_pro = std(units(cellNum).pro.neural.sacc.rate_pst)/sqrt(length(units(cellNum).pro.neural.trial));
        sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
        r_anti = units(cellNum).anti.neural.sacc.rate_pst;
        sem_anti = std(units(cellNum).anti.neural.sacc.rate_pst)/sqrt(length(units(cellNum).anti.neural.trial));
        sem_anti = repmat(sem_anti,[1 size(r_anti,2)]);
        mean_base= units(cellNum).pro.neural.base.rate_mu;
        mean_base = repmat(mean_base,[1 size(r_pro,2)]);
        
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
        shadedErrorBar(t,r_pro,sem_pro,'lineprops','r');
        shadedErrorBar(t,r_anti,sem_anti,'lineprops','g');
        plot(t,mean_base,'--k','LineWidth', 0.3);
        set (gca, 'xlim',([-0.1 0.2]), 'TickDir', 'out', 'FontSize',18); % analysis window size
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
        vline(0, 'k-');
        box off
        title('Aligned to saccade onset')
        
        %% instruction period
        %gather
         
        t= units(cellNum).pro.neural.sacc.ts_pst;
        r_pro= units(cellNum).pro.neural.instr.rate_pst;
        sem_pro = std(units(cellNum).pro.neural.instr.rate_pst)/sqrt(length(units(cellNum).pro.neural.trial));
        sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
        r_anti = units(cellNum).anti.neural.instr.rate_pst;
        sem_anti = std(units(cellNum).anti.neural.instr.rate_pst)/sqrt(length(units(cellNum).anti.neural.trial));
        sem_anti = repmat(sem_anti,[1 size(r_anti,2)]);
        
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
        title('Aligned to instruction onset')
        
    case 'psth_sacc_all'
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
        
        % plot w/sem
        figure; hold on;
        shadedErrorBar(t,r_pro,sem_pro,'lineprops','r');
        shadedErrorBar(t,r_anti,sem_anti,'lineprops','g');
        plot(t,mean_base,'--k','LineWidth', 0.3);
        set (gca, 'xlim',([-0.1 0.2]), 'TickDir', 'out', 'FontSize',18); % analysis window size
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
        vline(0, 'k-');
        box off
        title(['Aligned to saccade onset ' recArea ' unit= ' num2str(indx_area(i))])
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
        
        % plot w/sem
        figure; hold on;
        shadedErrorBar(t,r_pro,sem_pro,'lineprops','r');
        shadedErrorBar(t,r_anti,sem_anti,'lineprops','g');
        plot(t,mean_base,'--k','LineWidth', 0.3);
        set (gca, 'xlim',([-0.1 0.2]), 'TickDir', 'out', 'FontSize',18); % analysis window size
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
        vline(0, 'k-');
        box off
        title(['Aligned to instruction onset ' recArea ' unit= ' num2str(indx_area(i))])
        waitforbuttonpress; close all;
        end
        
    case 'colormap_sacc'  % TODO pick sacc related
        
        units = units; nunits = 1:length(units);
        %pro
        t = units(1).pro.neural.sacc.align_ts_pst; clear r;
        for j=1:length(nunits)
            r(j,:) = units(j).pro.neural.sacc.norm_rate_pst(1,:);
        end
        
        [r,t] = smooth_colormap(r,t);
        % plot colormap
        B = goodcolormap('bwr');
        figure; set(gcf,'Position',[100 200 300 300]);
        hold on; colormap(B');
        %hold on; colormap(bluewhitered);
        imagesc(t,1:nunits,r,[0,1]);
        set(gca,'xlim',[-0.5 0.5],'ylim',[1 nunits(end)],...
            'YTickLabel',[],'TickDir','Out','Fontsize',16);
        xlabel('Time (s)'); ylabel('Neuron');
        title('Pro')
        vline(0,'k')
        
        %anti
        t = units(1).anti.neural.sacc.align_ts_pst; clear r;
        for j=1:length(nunits)
            r(j,:) = units(j).anti.neural.sacc.norm_rate_pst(1,:);
        end
        
        [r,t] = smooth_colormap(r,t);
        % plot colormap
        B = goodcolormap('bwr');
        figure; set(gcf,'Position',[100 200 300 300]);
        hold on; colormap(B');
        %hold on; colormap(bluewhitered);
        imagesc(t,1:nunits,r,[0,1]);
        set(gca,'xlim',[-0.5 0.5],'ylim',[1 nunits(end)],...
            'YTickLabel',[],'TickDir','Out','Fontsize',16);
        xlabel('Time (s)'); ylabel('Neuron');
        title('Anti')
        vline(0,'k')
        
    case 'waterfall_sacc'  % TODO sort by max
        units = units; nunits = 1:length(units);
        %pro
        t = units(1).pro.neural.sacc.align_ts_pst; clear r;
        for j=1:length(nunits)
            r(j,:) = units(j).pro.neural.sacc.norm_rate_pst(1,:);
            indx_neg=find(r(j,:)<0);
            r(j,indx_neg) = 0;
        end
        
        figure; hold on; %colormap(d);
        for j=1:length(nunits)
            waterfall(t,j,r(j,:));
        end
        set(gca,'xlim', [-0.2 0.2], 'zlim', [0 2],'CameraViewAngle', 9, 'zTick', [], 'ytick', []);
        view(gca,[0 31]);
        grid off
    case 'scatter_pro_anti'
        %pro
        %gather
        rate_anti = units(cellNum).anti.neural.instr.nspk;
        rate_pro = units(cellNum).pro.neural.instr.nspk(1:size(rate_anti)); 
        
        %plot
        figure; hold on;
        plot(rate_pro, rate_anti, '.k','MarkerSize', 18);
        xlim([20 110]); ylim([20 110]); 
        title('Instruction'); xlabel('Prosaccade'); ylabel('Antisaccade');
        
    case 'scatter_pro_anti_pop'
        % plot avg firing rate pro vs anti for all cells
        %% instr
        % gather
        for i = 1:length(units)
            rate_pro(i) = units(i).pro.neural.instr.rate_mu; 
            rate_anti(i) = units(i).anti.neural.instr.rate_mu;
            sig_pro(i) = units(i).pro.neural.instr.rate_sig;
            sig_anti(i) = units(i).anti.neural.instr.rate_sig;
        end
        % get significantly diff
        for i = 1:length(units)
            indx_sign(i) = logical(units(i).stats.instr.flags.proVsAnti_instr);
        end
        
        % plot
        figure; hold on;
        plot(rate_pro(indx_sign),rate_anti(indx_sign), '.k','MarkerSize', 18);
        errorbar(rate_pro,rate_anti,sig_pro,sig_pro,sig_anti,sig_anti, 'ok', 'MarkerSize', 18,'Fontsize', 18);
        plot([0 150],[0 150]);
        title('Instruction'); xlabel('Prosaccade'); ylabel('Antisaccade');
        
        %% sacc
        for i = 1:length(units)
            rate_pro(i) = units(i).pro.neural.sacc.rate_mu; 
            rate_anti(i) = units(i).anti.neural.sacc.rate_mu;
            sig_pro(i) = units(i).pro.neural.sacc.rate_sig;
            sig_anti(i) = units(i).anti.neural.sacc.rate_sig;
        end
        % get significantly diff
        for i = 1:length(units)
            indx_sign(i) = logical(units(i).stats.sacc.flags.proVsAnti_sacc);
        end
        
         % plot
        figure; hold on;
        plot(rate_pro(indx_sign),rate_anti(indx_sign), '.k','MarkerSize', 14);
        errorbar(rate_pro,rate_anti,sig_pro,sig_pro,sig_anti,sig_anti, 'ok', 'MarkerSize', 18, 'Fontsize', 18);
        plot([0 150],[0 150]);
        title('Saccade'); xlabel('Prosaccade'); ylabel('Antisaccade');
  
        
        case 'peak_resp_instr'
        % gather
        t = units(i).pro.neural.instr.ts_pst_win;
        for i = 1:length(units)
            peak_pro(i) = units(i).pro.neural.instr.peak_resp;
            peak_anti(i) = units(i).anti.neural.instr.peak_resp;
            peak_time_pro(i) = units(i).pro.neural.instr.peak_resp_time;
            peak_time_anti(i) = units(i).anti.neural.instr.peak_resp_time;
        end

        % plot 
        figure; hold on; 
        plot(peak_pro,peak_anti, '.k', 'MarkerSize', 16);
        set(gca, 'xlim',[0 160], 'ylim',[0 160],'FontSize', 18, 'TickDir', 'out');
        plot([0 160],[0 160]); 
        title('Instr peak resp'); xlabel('Prosaccade'); ylabel('Antisaccade');
        box off;
        % plot peak time resp
        h_pro = hist(peak_time_pro,t);
        h_anti = hist(peak_time_anti,t);
        figure; hold on;
        plot(t,h_pro, 'r', 'LineWidth', 2);
        plot(t,h_anti, 'g', 'LineWidth', 2);
        set(gca, 'xlim',[0 25], 'TickDir', 'out'); box off;
        title('Peak resp time instr'); xlabel('Time (s)'); ylabel('Neuron nr')
        
    case 'peak_resp_sacc'
        % gather
        t = units(i).pro.neural.sacc.ts_pst_win;
        for i = 1:length(units)
            peak_pro(i) = units(i).pro.neural.sacc.peak_resp;
            peak_anti(i) = units(i).anti.neural.sacc.peak_resp;
            peak_time_pro(i) = units(i).pro.neural.sacc.peak_resp_time;
            peak_time_anti(i) = units(i).anti.neural.sacc.peak_resp_time;
        end
        % with baseline subtracted
        for i = 1:length(units)
            peak_pro_diff(i) = units(i).pro.neural.sacc.peak_resp-units(i).pro.neural.base.rate_mu;
            peak_anti_diff(i) = units(i).anti.neural.sacc.peak_resp-units(i).anti.neural.base.rate_mu;
        end
        
        % plot 
        figure; hold on; 
        plot(peak_pro,peak_anti, '.k', 'MarkerSize', 16);
        set(gca, 'xlim',[0 160], 'ylim',[0 160],'FontSize', 18, 'TickDir', 'out');
        plot([0 160],[0 160]); 
        title('Sacc peak resp'); xlabel('Prosaccade'); ylabel('Antisaccade');
        box off;
        % plot peak time resp
        h_pro = hist(peak_time_pro,t);
        h_anti = hist(peak_time_anti,t);
        figure; hold on;
        plot(t,h_pro, 'r', 'LineWidth', 2);
        plot(t,h_anti, 'g', 'LineWidth', 2);
        set(gca, 'xlim',[0 25], 'TickDir', 'out'); box off;
        title('Peak resp time'); xlabel('Time (s)'); ylabel('Neuron nr')
        
    case 'firingVSamp'
        %gather
        amp_pro = [units(cellNum).pro.behav.trial.saccAmplitude];
        amp_anti = [units(cellNum).anti.behav.trial.saccAmplitude];
        
        spk_pro = [units(cellNum).pro.neural.trial.nspk_sacc];
        spk_anti = [units(cellNum).anti.neural.trial.nspk_sacc];
        
        %plot
        figure; hold on; box off
        plot(amp_pro,spk_pro, '.r');
        plot(amp_anti,spk_anti, '.g');
        xlabel('Amplitude (deg)'); ylabel('Firing rate (spk/s)');
        set(gca, 'TickDir', 'out', 'FontSize', 18)
        
        %Avg per amplitude blocks
        block1 = 4.9;
        block2 = 9.9;
        block3 = 20;
        
        %pro
        block1_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]<=block1);
        block2_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]>block1 & [units(cellNum).pro.behav.trial.saccAmplitude]<= block2);
        block3_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]>block2 & [units(cellNum).pro.behav.trial.saccAmplitude]<= block3);
        
        for i = 1:length(block1_pro_indx)
            block1_pro(i,1) = units(cellNum).pro.neural.trial(block1_pro_indx(i)).nspk_sacc;
        end
        for i = 1:length(block2_pro_indx)
            block2_pro(i,1) = units(cellNum).pro.neural.trial(block2_pro_indx(i)).nspk_sacc;
        end
        for i = 1:length(block3_pro_indx)
            block3_pro(i,1) = units(cellNum).pro.neural.trial(block1_pro_indx(i)).nspk_sacc;
        end
        
        
        %anti
        block1_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]<=block1);
        block2_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]>block1 & [units(cellNum).anti.behav.trial.saccAmplitude]<= block2);
        block3_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]>block2 & [units(cellNum).anti.behav.trial.saccAmplitude]<= block3);
        
        if ~isempty(block1_anti_indx)
            for i = 1:length(block1_anti_indx)
                block1_anti(i,1) = units(cellNum).anti.neural.trial(block1_anti_indx(i)).nspk_sacc;
            end
        else
            block1_anti(i,1) = NaN;
        end
        for i = 1:length(block2_anti_indx)
            block2_anti(i,1) = units(cellNum).anti.neural.trial(block2_anti_indx(i)).nspk_sacc;
        end
        for i = 1:length(block3_anti_indx)
            block3_anti(i,1) = units(cellNum).anti.neural.trial(block3_anti_indx(i)).nspk_sacc;
        end
        
        % mean pro
        block1_pro_mu = nanmean(block1_pro); block1_pro_sig = nanstd(block1_pro);
        block2_pro_mu = nanmean(block2_pro); block2_pro_sig = nanstd(block2_pro);
        block3_pro_mu = nanmean(block3_pro); block1_pro_sig = nanstd(block3_pro);
        
        % mean anti
        block1_anti_mu = nanmean(block1_anti); block1_anti_sig = nanstd(block1_anti);
        block2_anti_mu = nanmean(block2_anti); block2_anti_sig = nanstd(block2_anti);
        block3_anti_mu = nanmean(block3_anti); block1_anti_sig = nanstd(block3_anti);
        
        %gather
        pro_blocks = [block1_pro_mu block2_pro_mu block3_pro_mu ; block1_pro_sig block2_pro_sig block1_pro_sig];
        anti_blocks = [block1_anti_mu block2_anti_mu block3_anti_mu ; block1_anti_sig block2_anti_sig block1_anti_sig];
        
        % plot means
        figure; hold on; box off
        errorbar(pro_blocks(1,:),pro_blocks(2,:), 'Color', 'r');
        errorbar(anti_blocks(1,:),anti_blocks(2,:), 'Color', 'g');
        
    case 'binomial_pb_dist'
        % instr
        % gather
         
        stat_instr = units(cellNum).stats.instr.pval.pbDist_testStat;
        t = units(cellNum).pro.neural.instr.ts_pst;
        % plot
        figure; subplot(2,1,1);
        plot(t, stat_instr, '.k','MarkerSize', 15);
        hline(-1.96, 'k');hline(1.96, 'k');vline(0, 'c');
        set(gca, 'xlim',[-0.1 0.3],'ylim',[-4 4], 'TickDir', 'out', 'FontSize', 18)
        title ('Binomial pb dist - Instruction')
        xlabel('time (s)')
        % sacc
        %gather
        t = units(cellNum).pro.neural.sacc.ts_pst;
        stat_sacc = units(cellNum).stats.sacc.pval.pbDist_testStat;
        % plot
        subplot(2,1,2)
        plot(t, stat_sacc, '.k','MarkerSize', 15);
        hline(-1.96, 'k');hline(1.96, 'k');vline(0, 'c');
        set(gca, 'xlim',[-0.3 0.2],'ylim',[-4 4], 'TickDir', 'out', 'FontSize', 18)
        title ('Binomial pb dist - Saccade')
        xlabel('time')
        
    case 'window_pb_dist'
        % instr
        % gather
        stat_instr = units(cellNum).stats.instr.pval.spk_count_bigWin; 
        % plot
        figure; subplot(2,1,1);
        plot(stat_instr, '.k','MarkerSize', 15);
        hline(-1.96, 'k');hline(1.96, 'k');vline(0, 'c');
        set(gca, 'ylim',[0 1], 'TickDir', 'out', 'FontSize', 18)
        title ('Binomial pb dist - Instruction')
        xlabel('')
        % sacc
        %gather
        stat_sacc = units(cellNum).stats.sacc.pval.spk_count_bigWin; 
        % plot
        subplot(2,1,2)
        plot(stat_sacc, '.k','MarkerSize', 15);
        hline(-1.96, 'k');hline(1.96, 'k');vline(0, 'c');
        set(gca, 'ylim',[0 1], 'TickDir', 'out', 'FontSize', 18)
        title ('Binomial pb dist - Saccade')
        xlabel('time')
        
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
       
        
        
end
