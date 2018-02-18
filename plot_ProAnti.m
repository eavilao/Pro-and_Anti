function plotUnit_ProAnti(trialData, cellNum, plotType)

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
    case 'raster'   %TODO Sort goCueTime and realign trials to that sorting
        
        %% Raster
        
        % saccade aligned
        % pro
        [~,indx] = sort([trialData(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[trialData(cellNum).pro.behav.trial(indx).reactionTime]; 
        subplot (2,1,1); hold on;box off
        for j=1:length(indx)
            if ~isempty(trialData(cellNum).pro.neural.trial(1).tspk_SS_align_sacc)
                plot(sorted_RT(j),j,'.r');
                plot(trialData(cellNum).pro.neural.trial(indx(j)).tspk_SS_align_sacc(1:2:end),j,'.k'); %plot every n spikes
                %plot(trialData(cellNum).pro.neural.trial(indx).tspk_SS_align_sacc,j,'.k');
            end
        end
        vline(0, 'c');
        set (gca, 'xlim', ([-0.5 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
        title('Prosaccade (aligned to saccade onset)');xlabel('Time (s)');ylabel('Trial Num')
        
        % anti
       [~,indx] = sort([trialData(cellNum).anti.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[trialData(cellNum).anti.behav.trial(indx).reactionTime]; 
        subplot (2,1,2); hold on;box off
        for j=1:length(indx)
            if ~isempty(trialData(cellNum).anti.neural.trial(1).tspk_SS_align_sacc)
                plot(sorted_RT(j),j,'.r');
                plot(trialData(cellNum).anti.neural.trial(indx(j)).tspk_SS_align_sacc(1:2:end),j,'.k'); %plot every n spikes
                %plot(trialData(cellNum).anti.neural.trial(indx).tspk_SS_align_sacc,j,'.k');
            end
        end
        vline(0, 'c');
        set (gca, 'xlim', ([-0.5 0.5]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
        title('Antisaccade (aligned to saccade onset)');xlabel('Time (s)');ylabel('Trial Num')
        
        
    case 'psth_saccade'
        
        %gather
        t= trialData(cellNum).pro.neural.ts_sacc;
        r_pro= trialData(cellNum).pro.neural.nspk_sacc;
        r_anti = trialData(cellNum).anti.neural.nspk_sacc;
        
        % smooth
        r_pro_smooth = smooth(r_pro,50);
        r_anti_smooth = smooth(r_anti,50);
        % Plot single units aligned to saccade
        figure;
        plot(t, r_pro_smooth, 'r', 'LineWidth', 3); hold on
        plot(t, r_anti_smooth, 'g', 'LineWidth', 3);
        set (gca, 'xlim',([-0.5 0.5]), 'TickDir', 'out', 'FontSize',18);
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s');
        vline(0, 'k--')
        box off
        
    case 'psth_instruction'
        %gather
        t= trialData(cellNum).pro.neural.ts_rate_instr;
        r_pro= trialData(cellNum).pro.neural.rate_instr;
        r_anti = trialData(cellNum).anti.neural.rate_instr;

        % smooth
        r_pro_smooth = smooth(r_pro,50);
        r_anti_smooth = smooth(r_anti,50);
        
        
        % Plot single units aligned to saccade
        plot(t, r_pro_smooth, 'r', 'LineWidth', 3); hold on
        plot(t, r_anti_smooth, 'g', 'LineWidth', 3);
        set (gca, 'xlim',([0 0.2]), 'TickDir', 'out', 'FontSize',18);
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s');
        box off
        title('Instruction')
        
        
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
        xlabel ('Amplitude (block per 5 deg)'); ylabel('Firing rate (spk/s)');
        set(gca, 'TickDir', 'out', 'FontSize', 18, 'xlim',([0.5 3.5]), 'ylim',([0 200]), 'XTick',[0 1 2 3]); 
        
    case 'colormap'
end
