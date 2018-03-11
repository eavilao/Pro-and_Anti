function plotUnit_ProAnti(trialData, plotType, cellNum)

% Needed: trialData.mat
% Input:    trialData  - output from extractWholeNeuronResults.m
%           cellNum - cell you want to plot
%           plotType - plot you want (e.g. raster, psth, etc.)


%%
% shadedErrorBar(trialData(1).anti.neural.sacc_ts_pst,trialData(1).anti.neural.sacc_rate_pst,repmat(trialData(1).anti.neural.sacc_rate_sig,1,21))
% hold on
% shadedErrorBar(trialData(1).pro.neural.sacc_ts_pst,trialData(1).pro.neural.sacc_rate_pst,repmat(trialData(1).pro.neural.sacc_rate_sig,1,21))

%%

switch plotType
    
    case 'raster_sacc'
        % saccade aligned
        % pro
        [~,indx] = sort([trialData(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[trialData(cellNum).pro.behav.trial(indx).reactionTime];
        r_pro= trialData(cellNum).pro.neural.trial;
        
        subplot (2,1,1); hold on;box off
        
        if strcmp(trialData(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_SS_align_sacc)
                    plot(sorted_RT(j),j,'.r');
                    %plot(r_pro(indx(j)).tspk_SS_align_sacc(1:2:end),j,'.k'); %plot every n spikes
                    plot(r_pro(indx(j)).tspk_SS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.5 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Prosaccade (aligned to saccade onset) SS');xlabel('Time (s)');ylabel('Trial Num')
        else
            for j= 1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_CS_align_sacc)
                    plot(sorted_RT(j),j,'.r');
                    plot(r_pro(indx(j)).tspk_CS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.5 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Prosaccade (aligned to saccade onset) CS');xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
        % anti
        [~,indx] = sort([trialData(cellNum).anti.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[trialData(cellNum).anti.behav.trial(indx).reactionTime];
        r_anti=trialData(cellNum).anti.neural.trial;
        subplot (2,1,2); hold on;box off
        
        if strcmp(trialData(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_SS_align_sacc)
                    plot(sorted_RT(j),j,'.r');
                    %plot(r_anti(indx(j)).tspk_SS_align_sacc(1:2:end),j,'.k'); %plot every n spikes
                    plot(r_anti(indx(j)).tspk_SS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.5 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Antisaccade (aligned to saccade onset) SS');xlabel('Time (s)');ylabel('Trial Num')
        else
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_CS_align_sacc)
                    plot(sorted_RT(j),j,'.r');
                    plot(r_anti(indx(j)).tspk_CS_align_sacc,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.5 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Antisaccade (aligned to saccade onset) CS');xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
    case 'raster_instr'
        % saccade aligned
        % pro
        [~,indx] = sort([trialData(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[trialData(cellNum).pro.behav.trial(indx).reactionTime];
        r_pro= trialData(cellNum).pro.neural.trial;
        
        subplot (2,1,1); hold on;box off
        
        if strcmp(trialData(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_pro(indx(j)).tspk_SS)
                    plot(sorted_RT(j),j,'.r');
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
                    plot(sorted_RT(j),j,'.r');
                    plot(r_pro(indx(j)).tspk_CS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Prosaccade (aligned to instruction onset) CS');xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
        % anti
        [~,indx] = sort([trialData(cellNum).anti.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[trialData(cellNum).anti.behav.trial(indx).reactionTime];
        r_anti=trialData(cellNum).anti.neural.trial;
        subplot (2,1,2); hold on;box off
        
        if strcmp(trialData(cellNum).id,'SS') % either SS or CS
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_SS)
                    plot(sorted_RT(j),j,'.r');
                    plot(r_anti(indx(j)).tspk_SS(1:2:end),j,'.k'); %plot every n spikes
                    %plot(trialData(cellNum).anti.neural.trial(indx).tspk_SS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Antisaccade (aligned to saccade onset) SS');xlabel('Time (s)');ylabel('Trial Num')
        else
            for j=1:length(indx)
                if ~isempty(r_anti(indx(j)).tspk_CS)
                    plot(sorted_RT(j),j,'.r');
                    plot(r_anti(indx(j)).tspk_CS,j,'.k');
                end
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.1 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
            title('Antisaccade (aligned to saccade onset) CS');xlabel('Time (s)');ylabel('Trial Num')
            
        end
        
        
    case 'psth_sacc'
        
        %gather
        t= trialData(cellNum).pro.neural.sacc.align_ts_pst;
        r_pro= trialData(cellNum).pro.neural.sacc.align_rate_pst;
        sem_pro = std(trialData(cellNum).pro.neural.sacc.align_rate_pst)/sqrt(length(trialData(cellNum).pro.neural.trial));
        sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
        r_anti = trialData(cellNum).anti.neural.sacc.align_rate_pst;
        sem_anti = std(trialData(cellNum).anti.neural.sacc.align_rate_pst)/sqrt(length(trialData(cellNum).anti.neural.trial));
        sem_anti = repmat(sem_anti,[1 size(r_anti,2)]);
        mean_base= trialData(cellNum).pro.neural.base.rate_mu;
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
        figure; hold on
        shadedErrorBar(t,r_pro,sem_pro,'lineprops','r');
        shadedErrorBar(t,r_anti,sem_anti,'lineprops','g');
        plot(t,mean_base,'--k','LineWidth', 0.3);
        set (gca, 'xlim',([-0.5 0.5]), 'TickDir', 'out', 'FontSize',18);
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s');
        vline(0, 'k-');
        box off
        title('Aligned to saccade onset')
        
        
    case 'psth_instr'
        %gather
        t= trialData(cellNum).pro.neural.sacc.align_ts_pst;
        r_pro= trialData(cellNum).pro.neural.rate_pst;
        sem_pro = std(trialData(cellNum).pro.neural.rate_pst)/sqrt(length(trialData(cellNum).pro.neural.trial));
        sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
        r_anti = trialData(cellNum).anti.neural.rate_pst;
        sem_anti = std(trialData(cellNum).anti.neural.rate_pst)/sqrt(length(trialData(cellNum).anti.neural.trial));
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
        figure; hold on
        shadedErrorBar(t, r_pro,sem_pro,'lineprops','r');
        shadedErrorBar(t, r_anti,sem_anti,'lineprops','g');
        set (gca, 'xlim',([-0.1 0.5]), 'TickDir', 'out', 'FontSize',18);
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s');
        vline(0, 'k--');
        box off
        title('Aligned to instruction onset')
        
    case 'colormap_sacc'
        units = trialData; nunits = 1:length(trialData);
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
        units = trialData; nunits = 1:length(trialData);
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
        
        
    case 'firingVSamp'
        %gather
        amp_pro = [trialData(cellNum).pro.behav.trial.saccAmplitude];
        amp_anti = [trialData(cellNum).anti.behav.trial.saccAmplitude];
        
        spk_pro = [trialData(cellNum).pro.neural.trial.nspk_sacc];
        spk_anti = [trialData(cellNum).anti.neural.trial.nspk_sacc];
        
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
        block1_pro_indx = find([trialData(cellNum).pro.behav.trial.saccAmplitude]<=block1);
        block2_pro_indx = find([trialData(cellNum).pro.behav.trial.saccAmplitude]>block1 & [trialData(cellNum).pro.behav.trial.saccAmplitude]<= block2);
        block3_pro_indx = find([trialData(cellNum).pro.behav.trial.saccAmplitude]>block2 & [trialData(cellNum).pro.behav.trial.saccAmplitude]<= block3);
        
        for i = 1:length(block1_pro_indx)
            block1_pro(i,1) = trialData(cellNum).pro.neural.trial(block1_pro_indx(i)).nspk_sacc;
        end
        for i = 1:length(block2_pro_indx)
            block2_pro(i,1) = trialData(cellNum).pro.neural.trial(block2_pro_indx(i)).nspk_sacc;
        end
        for i = 1:length(block3_pro_indx)
            block3_pro(i,1) = trialData(cellNum).pro.neural.trial(block1_pro_indx(i)).nspk_sacc;
        end
        
        
        %anti
        block1_anti_indx = find([trialData(cellNum).anti.behav.trial.saccAmplitude]<=block1);
        block2_anti_indx = find([trialData(cellNum).anti.behav.trial.saccAmplitude]>block1 & [trialData(cellNum).anti.behav.trial.saccAmplitude]<= block2);
        block3_anti_indx = find([trialData(cellNum).anti.behav.trial.saccAmplitude]>block2 & [trialData(cellNum).anti.behav.trial.saccAmplitude]<= block3);
        
        if ~isempty(block1_anti_indx)
            for i = 1:length(block1_anti_indx)
                block1_anti(i,1) = trialData(cellNum).anti.neural.trial(block1_anti_indx(i)).nspk_sacc;
            end
        else
            block1_anti(i,1) = NaN;
        end
        for i = 1:length(block2_anti_indx)
            block2_anti(i,1) = trialData(cellNum).anti.neural.trial(block2_anti_indx(i)).nspk_sacc;
        end
        for i = 1:length(block3_anti_indx)
            block3_anti(i,1) = trialData(cellNum).anti.neural.trial(block3_anti_indx(i)).nspk_sacc;
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
    
end
