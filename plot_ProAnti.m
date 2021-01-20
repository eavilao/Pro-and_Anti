function plot_ProAnti(units, pop, plotType, cellNum, recArea)

% Needed: units.mat
% Input:    units  - output from extractWholeNeuronResults.m
%           cellNum - cell you want to plot. If all leave empty = []
%           plotType - plot you want (e.g. raster, psth, etc.)
%           recArea - OMV, lateral Cb
%            -- for eye move -- load file sent by Nico

% if running population (all neurons), just inpu t [] in cellNum.
% if running single neuron, on recArea input []. Example: plot_ProAnti(units, 'raster_sacc',1,  [])

% List of plots

% 'eyeKin' : plot eye kinematics for all cells and trial
% 'kin_regress_cell' : plot eye kin multiple linear regression for each cell
% 'kin_regress_vermis' : plot eye kin multiple linear regression for all together
% 'kin_regress_lateral' : plot eye kin multiple linear regression for all together
% 'raster_sacc': saccade aligned raster plot for the chosen cell
% 'raster_sacc_single': plot raster for a single trial
% 'raster_instr': instruction aligned raster plot for the chosen cell
% 'raster_instr_single'
% 'psth': psth for the chosen cell, aligned to saccade and instruction
% 'pplsth_sacc_all': psth aligned to saccade onset for all cells. Press any key to plot next cell
% 'psth_instr_all': psth aligned to instruction onset for all cells. Press any key to plot next cell
% 'psth_mean_instr' :
% 'psth_mean': mean psth of two separate populations, exc and sup.
% 'psth_mean_selected' : mean psth of two separate populations, exc and sup for selected neurons
% 'psth_mean_instr_selected' : mean psth of two separate populations, exc and sup for selected neurons that also were signif during instruction period
% 'psth_mean_instr_selected_sorted' : sorted by instruction length.
% 'psth_mean_comparison_baselines'
% 'psth_mean_comparison_baselines_omv_lat'
% 'delta_rate_mean': plot mean change in firing rate for all and only significant cells
% 'delta_rate_mean_norm_pro': plot mean change in firing rate for all and only significant cells normalized by prosaccade rate.
% 'delta_rate_mean_norm_instr'
% 'change_in_fr_anti_pro_selected': change in FR from baseline (or comparison windows).
% 'change_in_fr_anti_pro_instr_selected'
% 'max_abs_change_fr'
% 'max_abs_change_fr_instr'
% 'max_change_instr_sacc_selected'
% 'change_anti-pro': absolute change in firing rate
% 'change_anti-pro_selected'
% 'change_anti-pro_instr'
% 'max_delta_rate'
% 'latencies_sacc_selected'
% 'delta_rate': first instance, contains all plots with time-course of change in firing rate (firing rate-baseline) for all cells and then plots max change in a scatter plot for significant cells.
% 'latencies_sacc': plot max period of response in sacc window
% 'latencies_instr': plot max period of response in instr window
% 'DDI_sacc'
% 'DDI_instr'
% 'DDI': discrimination index as described in Takahashi et al. 2007
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
% 'binomial_pb_dist_selected'
% 'binomial_pb_dist_selected_instr'
% 'binomial_pb_dist_signif'
% 'spk_pb_stat': spike probability density for 50 ms
% 'rosePlots'
% 'bi_resp' : cells responsive to both sacc and instruction
% 'cv'
% 'cv2'
% 'eye_move' plot eye traces for all trials for one neuron % LOAD FILE SENT BY NICO
% 'eye_move_single' plot single eye trace for a specific trial in a recording % LOAD FILE SENT BY NICO
% 'eye_move_instr'
% 'peak_vel_distr'
% 'rate_sorted_vel'
% 'mod_ratio': plot modulation ratio for all neurons and for signif diff neurons
% 'mod_ratio_selected'
% 'mod_depth' : plot modulation depth
% 'sorted_colormap_sacc'
% 'sorted_colormap_sacc_selected'
% 'sorted_colormap_sacc_selected_both_lateral'
% 'sorted_colormap_sacc_selected_both'
% 'sorted_colormap_instr'
% 'sorted_colormap_sacc_rand' (cosine sim plots here)
% 'sorted_colormap_sacc_rand_selected_exc'
% 'sorted_colormap_sacc_rand_selected_sup'
% 'latencies'
% 'mag_sensitivity_per_area'
% 'mag_sensitivity'
% 'rec_location'


%%
%%

switch plotType
    case 'eye_kin'
        
        %% plot
        % gather
        eyeKin = eyeKinematics_ProAnti(units);
        % eyeKin = eyeKinematics_ProAnti_perMonk(units); % MODIFY THIS !!! Here pick monkey/area/both.
        
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
        set(h1(1),'FaceColor', [1 0 1], 'EdgeColor', [1 0 1]);
        set(h2(1),'FaceColor', [0 0 1],'EdgeColor', [0 0 1]);
        set(h1(2),'Color',[1 0 1]);
        set(h2(2),'Color',[0 0 1]);
        xlim([0 25]);
        vline(mean(proAmp),'m'); % draw line on mean
        vline(mean(antiAmp),'b'); % draw line on mean
        [p,h,stats] = ranksum(proAmp, antiAmp);
        
        %plot dur
        figure; hold on
        h1 = histfit(proDur,20,'kernel');
        h2 = histfit(antiDur,20,'kernel');
        xlabel('Saccade duration (ms)')
        ylabel('Number of trials')
        set (gca, 'TickDir', 'out','FontSize', 18);
        alpha(0.25)
        set(h1(1),'FaceColor', [1 0 1], 'EdgeColor', [1 0 1]);
        set(h2(1),'FaceColor', [0 0 1],'EdgeColor', [0 0 1]);
        set(h1(2),'Color',[1 0 1]);
        set(h2(2),'Color',[0 0 1]);
        xlim([0 150]);
        vline(mean(proDur),'m'); % draw line on mean
        vline(mean(antiDur),'b'); % draw line on mean
        [p,h] = ranksum(proDur, antiDur);
        
        % plot PV
        figure; hold on
        h1 = histfit(proPV,20,'kernel');
        h2 = histfit(antiPV,20,'kernel');
        xlabel('Saccade peak velocity (deg/s)')
        ylabel('Number of trials')
        set (gca, 'TickDir', 'out','FontSize', 18);
        alpha(0.25)
        set(h1(1),'FaceColor', [1 0 1], 'EdgeColor', [1 0 1]);
        set(h2(1),'FaceColor', [0 0 1],'EdgeColor', [0 0 1]);
        set(h1(2),'Color',[1 0 1]);
        set(h2(2),'Color',[0 0 1]);
        xlim([0 700]);
        vline(mean(proPV),'m'); % draw line on mean
        vline(mean(antiPV),'b'); % draw line on mean
        [p,h] = ranksum(proPV, antiPV);
        
        % plot RT
        figure; hold on
        h1 = histfit(proRT,20,'kernel');
        h2 = histfit(antiRT,20,'kernel');
        xlabel('Reaction time (ms)')
        ylabel('Number of trials')
        set (gca, 'TickDir', 'out','FontSize', 18);
        alpha(0.25)
        set(h1(1),'FaceColor', [1 0 1], 'EdgeColor', [1 0 1]);
        set(h2(1),'FaceColor', [0 0 1],'EdgeColor', [0 0 1]);
        set(h1(2),'Color',[1 0 1]);
        set(h2(2),'Color',[0 0 1]);
        xlim([0 500]);
        vline(mean(proRT),'m'); % draw line on mean
        vline(mean(antiRT),'b'); % draw line on mean
        [p,h] = ranksum(proRT, antiRT);
        
    case 'kin_regress_cell'
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        for u = 1:length(indx_area)
            r_pro(u) = mean(pop.kin(u).pro.r_all_pro);
            r_anti(u) = mean(pop.kin(u).anti.r_all_anti);
            coeff_pro(:,u) = units(indx_area(u)).stats.sacc.regress.coeff_pro; 
            coeff_anti(:,u) = units(indx_area(u)).stats.sacc.regress.coeff_anti;
            
            x1_pro(:,u) = units(indx_area(u)).stats.sacc.regress.coeff_pro(2);
            
            
            
            x1_anti(:,u) = units(indx_area(u)).stats.sacc.regress.coeff_anti();
            
        end
        
        
        
        
        LM_pro = stepwiselm([x1_pro x2_pro x3_pro x4_pro], y_pro)
        step_LM_pro = stepwise([x1_pro x2_pro x3_pro x4_pro], y_pro)
        
        LM_anti = stepwiselm([x1_anti x2_anti x3_anti x4_anti], y_anti)
        step_LM_anti = stepwise([x1_anti x2_anti x3_anti x4_anti], y_anti)
        
        
    case 'kin_regress_vermis'
        x1_pro = pop.pro.vermis.kin_all(:,1); % amp
        x2_pro = pop.pro.vermis.kin_all(:,2); % dur
        x3_pro = pop.pro.vermis.kin_all(:,3); % pv
        x4_pro = pop.pro.vermis.kin_all(:,4)*1000; % reaction time
        y_pro = pop.pro.vermis.r_all;
        coeff_pro = pop.stats.sacc.pro.regress.vermis.coeff_pro;
        
        [~,hfig] = plotmatrix([y_pro x1_pro x2_pro x3_pro x4_pro]);
        hfig(1,1).YLabel.String = 'firing'; hfig(5,1).XLabel.String = 'firing';
        hfig(2,1).YLabel.String = 'amp'; hfig(5,2).XLabel.String = 'amp';
        hfig(3,1).YLabel.String = 'dur'; hfig(5,3).XLabel.String = 'dur';
        hfig(4,1).YLabel.String = 'pv'; hfig(5,4).XLabel.String = 'pv';
        hfig(5,1).YLabel.String = 'rt'; hfig(5,5).XLabel.String = 'rt';
        
        LM_pro = stepwiselm([x1_pro x2_pro x3_pro x4_pro], y_pro)
        stepwise([x1_pro x2_pro x3_pro x4_pro], y_pro)
        
        scatter3(x1_pro,x2_pro,y_pro);
        hold on;
        x1fit_pro = linspace(min(x1_pro),max(x1_pro),25);
        x2fit_pro = linspace(min(x2_pro),max(x2_pro),25);
        [X1FIT_pro,X2FIT_pro] = meshgrid(x1fit_pro,x2fit_pro);
        YFIT_pro = coeff_pro(1) + coeff_pro(2)*X1FIT_pro + coeff_pro(3)*X2FIT_pro + coeff_pro(4)*X1FIT_pro.*X2FIT_pro;
        mesh(X1FIT_pro,X2FIT_pro,YFIT_pro);
        
        x1_anti = pop.anti.vermis.kin_all(:,1);
        x2_anti = pop.anti.vermis.kin_all(:,2);
        x3_anti = pop.anti.vermis.kin_all(:,3);
        x4_anti = pop.anti.vermis.kin_all(:,4)*1000;
        y_anti = pop.anti.vermis.r_all;
        coeff_pro = pop.stats.sacc.anti.regress.vermis.coeff_anti;
        
        LM_anti =stepwiselm([x1_anti x2_anti x3_anti x4_anti], y_anti)
        stepwise([x1_anti x2_anti x3_anti x4_anti], y_anti)
        
        [~,hfig] = plotmatrix([y_pro x1_pro x2_pro x3_pro x4_pro]);
        hfig(1,1).YLabel.String = 'firing'; hfig(5,1).XLabel.String = 'firing';
        hfig(2,1).YLabel.String = 'amp'; hfig(5,2).XLabel.String = 'amp';
        hfig(3,1).YLabel.String = 'dur'; hfig(5,3).XLabel.String = 'dur';
        hfig(4,1).YLabel.String = 'pv'; hfig(5,4).XLabel.String = 'pv';
        hfig(5,1).YLabel.String = 'rt'; hfig(5,5).XLabel.String = 'rt';
        
        scatter3(x1_anti,x2_anti,y_anti);
        hold on;
        x1fit_anti = min(x1_anti):1:max(x1_anti);
        x2fit_anti = min(x2_anti):1:max(x2_anti);
        [X1FIT_anti,X2FIT_anti] = meshgrid(x1fit_anti,x2fit_anti);
        YFIT_anti = coeff_anti(1) + coeff_anti(2)*X1FIT_anti + coeff_anti(3)*X2FIT_anti + coeff_anti(4)*X1FIT_anti.*X2FIT_anti;
        mesh(X1FIT_anti,X2FIT_anti,YFIT_anti);
        
    case 'kin_regress_lateral'
        x1_pro = pop.pro.lateral.kin_all(:,1);
        x2_pro = pop.pro.lateral.kin_all(:,2);
        x3_pro = pop.pro.lateral.kin_all(:,3);
        x4_pro = pop.pro.lateral.kin_all(:,4)*1000;
        y_pro = pop.pro.lateral.r_all;
        coeff_pro = pop.stats.sacc.pro.regress.lateral.coeff_pro;
        
        L = stepwise([x1_pro x2_pro x3_pro x4_pro], y_pro)
        LM_pro = stepwiselm([x1_pro x2_pro x3_pro x4_pro], y_pro)
        
        [~,hfig] = plotmatrix([y_pro x1_pro x2_pro x3_pro x4_pro]);
        hfig(1,1).YLabel.String = 'firing'; hfig(5,1).XLabel.String = 'firing';
        hfig(2,1).YLabel.String = 'amp'; hfig(5,2).XLabel.String = 'amp';
        hfig(3,1).YLabel.String = 'dur'; hfig(5,3).XLabel.String = 'dur';
        hfig(4,1).YLabel.String = 'pv'; hfig(5,4).XLabel.String = 'pv';
        hfig(5,1).YLabel.String = 'rt'; hfig(5,5).XLabel.String = 'rt';
        
        scatter3(x1_pro,x2_pro,y_pro);
        hold on;
        x1fit_pro = linspace(min(x1_pro),max(x1_pro),25);
        x2fit_pro = linspace(min(x2_pro),max(x2_pro),25);
        [X1FIT_pro,X2FIT_pro] = meshgrid(x1fit_pro,x2fit_pro);
        YFIT_pro = coeff_pro(1) + coeff_pro(2)*X1FIT_pro + coeff_pro(3)*X2FIT_pro + coeff_pro(4)*X1FIT_pro.*X2FIT_pro;
        mesh(X1FIT_pro,X2FIT_pro,YFIT_pro);
        
        x1_anti = pop.anti.lateral.kin_all(:,1);
        x2_anti = pop.anti.lateral.kin_all(:,2);
        x3_anti = pop.anti.lateral.kin_all(:,3);
        x4_anti = pop.anti.lateral.kin_all(:,4)*1000;
        y_anti = pop.anti.lateral.r_all;
        coeff_pro = pop.stats.sacc.anti.regress.lateral.coeff_anti;
        
        stepwise([x1_anti x2_anti x3_anti x4_anti], y_anti)
        % stepwiselm([x1_anti x2_anti x3_anti x4_anti], y_anti)
        
        [~,hfig] = plotmatrix([y_pro x1_pro x2_pro x3_pro x4_pro]);
        hfig(1,1).YLabel.String = 'firing'; hfig(5,1).XLabel.String = 'firing';
        hfig(2,1).YLabel.String = 'amp'; hfig(5,2).XLabel.String = 'amp';
        hfig(3,1).YLabel.String = 'dur'; hfig(5,3).XLabel.String = 'dur';
        hfig(4,1).YLabel.String = 'pv'; hfig(5,4).XLabel.String = 'pv';
        hfig(5,1).YLabel.String = 'rt'; hfig(5,5).XLabel.String = 'rt';
        
        scatter3(x1_anti,x2_anti,y_anti);
        hold on;
        x1fit_anti = min(x1_anti):1:max(x1_anti);
        x2fit_anti = min(x2_anti):1:max(x2_anti);
        [X1FIT_anti,X2FIT_anti] = meshgrid(x1fit_anti,x2fit_anti);
        YFIT_anti = coeff_anti(1) + coeff_anti(2)*X1FIT_anti + coeff_anti(3)*X2FIT_anti + coeff_anti(4)*X1FIT_anti.*X2FIT_anti;
        mesh(X1FIT_anti,X2FIT_anti,YFIT_anti);
        
        
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
        
        print('raster_sacc','-depsc2', '-painters', '-cmyk')
        
    case 'raster_sacc_single'
        [~,indx] = sort([units(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).pro.behav.trial(indx).reactionTime];
        r= units(cellNum).pro.neural.trial; % pro
        % r= units(cellNum).anti.neural.trial; % anti
        
        recArea = units(cellNum).area;
        
        for j=1:length(indx)
            figure('Position',[2454 782 350 119]);hold on;
            for ii = 1:length(r(indx(j)).tspk_SS_align_sacc)
                line([r(indx(j)).tspk_SS_align_sacc(ii) r(indx(j)).tspk_SS_align_sacc(ii)], [0 1],'Color', 'k', 'LineWidth',2)
            end
            vline(0, 'c');
            set (gca, 'xlim', ([-0.150 0.151]), 'TickDir', 'out', 'FontSize', 22); title(num2str(j));
            waitforbuttonpress; close all;
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
            set (gca, 'xlim', ([0 0.35]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
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
            set (gca, 'xlim', ([0 0.35]), 'ylim',([0 j]), 'TickDir', 'out', 'FontSize', 18);
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
        
        % print('raster','-depsc2', '-painters', '-cmyk')
        
    case 'raster_instr_single'
        [~,indx] = sort([units(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).pro.behav.trial(indx).reactionTime];
        r= units(cellNum).pro.neural.trial; % pro
        %         r= units(cellNum).anti.neural.trial; % anti
        
        recArea = units(cellNum).area;
        
        for j=1:length(indx)
            figure('Position',[2454 782 350 119]);hold on;
            for ii = 1:length(r(indx(j)).tspk_SS)
                line([r(indx(j)).tspk_SS(ii) r(indx(j)).tspk_SS(ii)], [0 1],'Color', 'k', 'LineWidth',2)
            end
            vline(0, 'c');
            set (gca, 'xlim', ([0 0.350]), 'TickDir', 'out', 'FontSize', 22); title(num2str(j));
            waitforbuttonpress; close all;
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
        t_instr = units(cellNum).pro.neural.instr.ts_pst;
        r_pro= units(cellNum).pro.neural.instr.rate_pst;
        sem_pro = std(units(cellNum).pro.neural.instr.rate_pst)/sqrt(length(units(cellNum).pro.neural.trial));
        std_pro = std(units(cellNum).pro.neural.instr.rate_pst); std_pro = repmat(std_pro,[1 size(r_pro,2)]);
        sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
        r_anti = units(cellNum).anti.neural.instr.rate_pst;
        sem_anti = std(units(cellNum).anti.neural.instr.rate_pst)/sqrt(length(units(cellNum).anti.neural.trial));
        std_anti = std(units(cellNum).anti.neural.instr.rate_pst); std_anti = repmat(std_anti,[1 size(r_anti,2)]);
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
        s_pro = shadedErrorBar(t_instr, r_pro,std_pro,'lineprops','r');
        s_anti = shadedErrorBar(t_instr, r_anti,std_anti,'lineprops','g');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set (gca, 'xlim',([0 0.35]), 'xTick', [0 .175 0.350] , 'TickDir', 'out', 'FontSize',22);
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
        %vline(0, 'k--');
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
            shadedErrorBar(t,r_pro,sem_pro,'lineprops','m');
            shadedErrorBar(t,r_anti,sem_anti,'lineprops','b');
            plot(t,mean_base,'--k','LineWidth', 0.3);
            set (gca, 'xlim',([-0.15 0.15]), 'TickDir', 'out', 'FontSize',18); % analysis window size
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
            std_pro = std(units(indx_area(i)).pro.neural.instr.rate_pst); std_pro = repmat(std_pro,[1 size(r_pro,2)]);
            sem_pro = std(units(indx_area(i)).pro.neural.instr.rate_pst)/sqrt(length(units(indx_area(i)).pro.neural.trial));
            sem_pro = repmat(sem_pro,[1 size(r_pro,2)]);
            r_anti = units(indx_area(i)).anti.neural.instr.rate_pst;
            std_anti = std(units(indx_area(i)).anti.neural.instr.rate_pst); std_anti = repmat(std_anti,[1 size(r_pro,2)]);
            sem_anti = std(units(indx_area(i)).anti.neural.instr.rate_pst)/sqrt(length(units(indx_area(i)).anti.neural.trial));
            sem_anti = repmat(sem_anti,[1 size(r_anti,2)]);
            mean_base= units(indx_area(i)).pro.neural.base.rate_mu;
            mean_base = repmat(mean_base,[1 size(r_pro,2)]);
            proVSanti_instr = units(indx_area(i)).stats.instr.flags.proVsAnti_instr;
            
            % plot w/sem
            figure; hold on;
            s_pro = shadedErrorBar(t,r_pro,sem_pro,'lineprops','m');
            s_anti = shadedErrorBar(t,r_anti,sem_anti,'lineprops','b');
            set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
            set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
            set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
            set (gca, 'xlim',([0 0.35]), 'TickDir', 'out', 'FontSize',18); % analysis window size
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
        
    case 'psth_mean_instr'
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc_instr==1
                indx_exc(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup_instr==1
                indx_sup(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        [~,indx] = sort([units(cellNum).pro.behav.trial.reactionTime],'descend'); % sort RT
        sorted_RT = -[units(cellNum).pro.behav.trial(indx).reactionTime];
        
        % get exc
        t = units(1).pro.neural.instr.ts_pst_win;
        for i = 1:length(indx_exc)
            r_exc_pro(i,:) = units(indx_exc(i)).pro.neural.instr.rate_pst_win;
            std_exc_pro(i) = std(units(indx_exc(i)).pro.neural.instr.rate_pst_win);
            sem_exc_pro(i,:)= std(units(indx_exc(i)).pro.neural.instr.rate_pst_win)/sqrt(length(indx_exc));
            r_exc_anti(i,:) = units(indx_exc(i)).anti.neural.instr.rate_pst_win;
            std_exc_anti(i,:) = std(units(indx_exc(i)).anti.neural.instr.rate_pst_win);
            sem_exc_anti(i,:)= std(units(indx_exc(i)).anti.neural.instr.rate_pst_win)/sqrt(length(indx_sup));
        end
        
        % get supp
        for i = 1:length(indx_sup)
            r_sup_pro(i,:) = units(indx_sup(i)).pro.neural.instr.rate_pst_win;
            std_sup_pro(i,:) = std(units(indx_sup(i)).pro.neural.instr.rate_pst_win);
            sem_sup_pro(i,:)= std(units(indx_sup(i)).pro.neural.instr.rate_pst_win)/sqrt(length(indx_sup));
            r_sup_anti(i,:) = units(indx_sup(i)).anti.neural.instr.rate_pst_win;
            std_sup_anti(i,:) = std(units(indx_sup(i)).anti.neural.instr.rate_pst_win);
            sem_sup_anti(i,:)= std(units(indx_sup(i)).anti.neural.instr.rate_pst_win)/sqrt(length(indx_sup));
        end
        
        % plot exc
        figure; hold on;
        plot(t,mean(r_exc_pro));
        plot(t,mean(r_exc_anti));
        s_pro = shadedErrorBar(t, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[0.06 0.350], 'ylim', [65 75], 'ytick', [45 90],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Firing rate (spk/s)'); xlabel('Time (s)')
        
        % plot sup
        figure; hold on;
        plot(t,mean(r_sup_pro));
        plot(t,mean(r_sup_anti));
        s_pro = shadedErrorBar(t, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[0.06 0.350], 'ylim', [45 65], 'ytick', [50 80],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Firing rate (spk/s)'); xlabel('Time (s)')
        
    case 'psth_mean'
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc==1
                indx_exc_pro(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup==1
                indx_sup_pro(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).anti.neural.exc==1
                indx_exc_anti(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).anti.neural.sup==1
                indx_sup_anti(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=-0.151 & t<=0.151);
        for i = 1:length(indx_exc_pro)
            r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151));
            sem_exc_pro(i,:)= std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); 
            r_exc_anti(i,:) = units(indx_exc_pro(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151));
            sem_exc_anti(i,:)= std(r_exc_anti(i,:))/sqrt(length(indx_exc_pro)); 
        end
        
        for i = 1:length(indx_exc_anti)
            r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_exc_anti(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151));
            sem_exc_anti(i,:)= std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); 
        end
        
        
        % get supp
        for i = 1:length(indx_sup_pro)
            r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151));
            sem_sup_pro(i,:)= std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); 
        end
        
        for i = 1:length(indx_sup_anti)
            r_sup_anti(i,:) = units(indx_sup_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_sup_anti(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151));
            sem_sup_anti(i,:)= std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti));
        end
        
        % plot exc
        figure; hold on;
        plot(t_win,mean(r_exc_pro));
        plot(t_win,mean(r_exc_anti));
        s_pro = shadedErrorBar(t_win, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [0.5 1.035], 'ytick', [0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Firing rate (spk/s)'); xlabel('Time (s)')
        
        % plot sup
        figure; hold on;
        plot(t_win,mean(r_sup_pro));
        plot(t_win,mean(r_sup_anti));
        s_pro = shadedErrorBar(t_win, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 1), set(s_anti.mainLine,'LineWidth', 1);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [0.5 1.035], 'ytick', [0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Firing rate (spk/s)'); xlabel('Time (s)')
        
    case 'psth_mean_selected'
        % gather
        indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=-0.151 & t<=0.151);
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)); ...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)); ...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_exc_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)); ...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti(i,:) = units(indx_sup_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_sup_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)); ...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        
        % plot exc
        figure; hold on;
        plot(t_win,mean(r_exc_pro));
        plot(t_win,mean(r_exc_anti));
        s_pro = shadedErrorBar(t_win, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [0.5 1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate'); xlabel('Time (s)')
        
        
        % plot sup
        figure; hold on;
        plot(t_win,mean(r_sup_pro));
        plot(t_win,mean(r_sup_anti));
        s_pro = shadedErrorBar(t_win, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [0.5 1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate'); xlabel('Time (s)')
        
        
        case 'psth_mean_selected_norm_on_pro'
        % gather
        indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=-0.151 & t<=0.151);
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)); ...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)); ...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_exc_anti(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)); ...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti(i,:) = units(indx_sup_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)/max(units(indx_sup_anti(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0.151)); ...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        
        % plot exc
        figure; hold on;
        plot(t_win,mean(r_exc_pro));
        plot(t_win,mean(r_exc_anti));
        s_pro = shadedErrorBar(t_win, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [0.5 1.1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate (to prosaccades)'); xlabel('Time (s)')
        
        
        % plot sup
        figure; hold on;
        plot(t_win,mean(r_sup_pro));
        plot(t_win,mean(r_sup_anti));
        s_pro = shadedErrorBar(t_win, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [0.5 1.03], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate (to prosaccades)'); xlabel('Time (s)')
        
        
    case 'psth_mean_instr_selected_norm_on_pro'
        indx_exc_pro = pop.indx_sel.(recArea).instr_back.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).instr_back.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).instr_back.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).instr_back.all.anti.sup;
        
        t_win = units(1).pro.neural.instr_back.ts_pst_win;
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win/max(units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win);,...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win/max(units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win); ...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.instr_back.rate_pst_win/max(units(indx_exc_anti(i)).pro.neural.instr_back.rate_pst_win); ...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti(i,:) = units(indx_sup_anti(i)).anti.neural.instr_back.rate_pst_win/max(units(indx_sup_anti(i)).pro.neural.instr_back.rate_pst_win); ...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        
        % plot exc
        figure; hold on;
        plot(t_win,mean(r_exc_pro));
        plot(t_win,mean(r_exc_anti));
        s_pro = shadedErrorBar(t_win, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.3 0.1], 'ylim', [0.5 1.1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate on prosaccades'); xlabel('Time from dir instruction (s)')
        
        % plot sup
        figure; hold on;
        plot(t_win,mean(r_sup_pro));
        plot(t_win,mean(r_sup_anti));
        s_pro = shadedErrorBar(t_win, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.3 0.1], 'ylim', [0.5 1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate on prosaccades'); xlabel('Time from dir instruction (s)')
        
        case 'psth_mean_instr_selected'
        indx_exc_pro = pop.indx_sel.(recArea).instr_back.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).instr_back.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).instr_back.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).instr_back.all.anti.sup;
        
        t_win = units(1).pro.neural.instr_back.ts_pst_win;
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win/max(units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win);,...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win/max(units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win); ...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.instr_back.rate_pst_win/max(units(indx_exc_anti(i)).anti.neural.instr_back.rate_pst_win); ...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti(i,:) = units(indx_sup_anti(i)).anti.neural.instr_back.rate_pst_win/max(units(indx_sup_anti(i)).anti.neural.instr_back.rate_pst_win); ...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        
        % plot exc
        figure; hold on;
        plot(t_win,mean(r_exc_pro));
        plot(t_win,mean(r_exc_anti));
        s_pro = shadedErrorBar(t_win, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.3 0], 'ylim', [0.5 1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate'); xlabel('Time from dir instruction (s)')
        
        % plot sup
        figure; hold on;
        plot(t_win,mean(r_sup_pro));
        plot(t_win,mean(r_sup_anti));
        s_pro = shadedErrorBar(t_win, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.3 0], 'ylim', [0.5 1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate'); xlabel('Time from dir instruction (s)')
        
        
    case 'psth_mean_instr_dir_selected'
        %% plot with instruction direction
        indx_exc_pro = pop.indx_sel.(recArea).instr_back.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).instr_back.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).instr_back.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).instr_back.all.anti.sup;
        
        win = [-0.3 0.1]; 
        t = units(indx_exc_pro(1)).pro.neural.instrDir.ts_pst; t_win = units(indx_exc_pro(1)).pro.neural.instrDir.ts_pst(t>=win(1) & t<=win(2));
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.instrDir.rate_pst(t>=win(1) & t<=win(2))/max(units(indx_exc_pro(i)).pro.neural.instrDir.rate_pst(t>=win(1) & t<=win(2)));,...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.instrDir.rate_pst(t>=win(1) & t<=win(2))/max(units(indx_sup_pro(i)).pro.neural.instrDir.rate_pst(t>=win(1) & t<=win(2)));,...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.instrDir.rate_pst(t>=win(1) & t<=win(2))/max(units(indx_exc_anti(i)).pro.neural.instrDir.rate_pst(t>=win(1) & t<=win(2)));,...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
         for i = 1:length(indx_sup_anti), r_sup_anti(i,:) = units(indx_sup_anti(i)).anti.neural.instrDir.rate_pst(t>=win(1) & t<=win(2))/max(units(indx_sup_anti(i)).pro.neural.instrDir.rate_pst(t>=win(1) & t<=win(2)));,...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        
        % plot exc
        figure; hold on;
        plot(t_win,mean(r_exc_pro));
        plot(t_win,mean(r_exc_anti));
        s_pro = shadedErrorBar(t_win, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.3 0.1], 'ylim', [0.5 1.5], 'ytick', [0 0.5 1 1.5],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate'); xlabel('Time from dir instruction (s)')
        
        % plot sup
        figure; hold on;
        plot(t_win,mean(r_sup_pro));
        plot(t_win,mean(r_sup_anti));
        s_pro = shadedErrorBar(t_win, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.3 0.1], 'ylim', [0.5 1.034], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Normalized firing rate'); xlabel('Time from dir instruction (s)')
        
        
    case 'psth_mean_comparison_baselines'
        indx_exc = pop.indx_sel.(recArea).sacc.all.exc;
        indx_sup = pop.indx_sel.(recArea).sacc.all.sup;
        
        for i = 1:length(indx_exc), r_exc_pro(i) = units(indx_exc(i)).pro.neural.sacc.rate_mu - units(indx_exc(i)).pro.neural.base.rate_mu; end
        for i = 1:length(indx_exc), r_exc_pro_base(i) = units(indx_exc(i)).pro.neural.base.rate_mu; end
        for i = 1:length(indx_exc), r_exc_anti(i) = units(indx_exc(i)).anti.neural.sacc.rate_mu - units(indx_exc(i)).anti.neural.base.rate_mu; end
        for i = 1:length(indx_exc), r_exc_anti_base(i) = units(indx_exc(i)).anti.neural.base.rate_mu; end
        
        for i = 1:length(indx_sup), r_sup_pro(i) = units(indx_sup(i)).pro.neural.sacc.rate_mu - units(indx_sup(i)).pro.neural.base.rate_mu; end
        for i = 1:length(indx_sup), r_sup_pro_base(i) = units(indx_sup(i)).pro.neural.base.rate_mu; end
        for i = 1:length(indx_sup), r_sup_anti(i) = units(indx_sup(i)).anti.neural.sacc.rate_mu - units(indx_sup(i)).anti.neural.base.rate_mu; end
        for i = 1:length(indx_sup), r_sup_anti_base(i) = units(indx_sup(i)).anti.neural.base.rate_mu; end
        
        % plot exc pro vs anti
        figure; hold on; plot(r_exc_pro_base,r_exc_anti_base, '.k','MarkerSize',30);plot([0 100],[0 100], 'k'); xlabel('Pro'); ylabel('Anti'); title('pro vs anti base for exc');  axis square; 
        figure; hold on; plot(r_exc_pro,r_exc_anti, '.k','MarkerSize',30); plot([0 80],[0 80], 'k'); title('pro vs anti diff from base for exc')
        xlabel('Pro'); ylabel('Anti'); set(gca, 'xlim', [0 80], 'ylim', [0 80]);  axis square; 
        
        % plot sup pro vs anti
        figure; hold on; plot(r_sup_pro_base,r_sup_anti_base, 'ok', 'MarkerSize',10); plot([0 180],[0 180], 'k'); xlabel('Pro'); ylabel('Anti'); title('pro vs anti base for sup');  axis square; 
        figure; hold on; plot(abs(r_sup_pro),abs(r_sup_anti), 'ok','MarkerSize',10); plot([0 40],[0 40], 'k'); title('pro vs anti diff from base for sup')
        xlabel('Pro'); ylabel('Anti'); set(gca, 'xlim', [0 40], 'ylim', [0 40]);  axis square; 
        
        [~,p_exc] = ttest2(r_exc_pro, r_exc_anti)
        [~,p_sup] = ttest2(abs(r_sup_pro), abs(r_sup_anti))
        
        % Plot all exc vs sup
        figure; hold on;  
        plot(r_exc_pro_base,r_exc_anti_base, '.k','MarkerSize',30);
        plot(r_sup_pro_base,r_sup_anti_base, 'ok', 'MarkerSize',10);plot([0 180],[0 180], 'k');
        xlabel('Pro'); ylabel('Anti');  title('pro vs anti base');
        axis square; 
        
        figure; hold on;
        plot(r_exc_pro,r_exc_anti, '.k','MarkerSize',30);
        plot(abs(r_sup_pro), abs(r_sup_anti), 'ok','MarkerSize',10);
        plot([0 80],[0 80], 'k'); axis square; 
        xlabel('Pro'); ylabel('Anti'); set(gca, 'xlim', [0 80], 'ylim', [0 80]);
        
        ranksum(r_exc_pro, abs(r_sup_pro))
        ranksum(r_exc_anti, abs(r_sup_anti))
        
        
    case 'psth_mean_comparison_baselines_omv_lat'
       indx_exc_omv = pop.indx_sel.vermis.sacc.all.exc;  indx_exc_lat = pop.indx_sel.lateral.sacc.all.exc;
       indx_sup_omv = pop.indx_sel.vermis.sacc.all.sup;  indx_sup_lat = pop.indx_sel.lateral.sacc.all.sup;
       
       for i = 1:length(indx_exc_omv), r_exc_pro_omv(i) = units(indx_exc_omv(i)).pro.neural.base.rate_mu; end ; for i = 1:length(indx_exc_omv), r_exc_anti_omv(i) = units(indx_exc_omv(i)).anti.neural.base.rate_mu; end
       for i = 1:length(indx_sup_omv), r_sup_pro_omv(i) = units(indx_sup_omv(i)).pro.neural.base.rate_mu; end ;  for i = 1:length(indx_sup_omv), r_sup_anti_omv(i) = units(indx_sup_omv(i)).anti.neural.base.rate_mu; end
       
       for i = 1:length(indx_exc_lat), r_exc_pro_lat(i) = units(indx_exc_lat(i)).pro.neural.base.rate_mu; end ; for i = 1:length(indx_exc_lat), r_exc_anti_lat(i) = units(indx_exc_lat(i)).anti.neural.base.rate_mu; end
       for i = 1:length(indx_sup_lat), r_sup_pro_lat(i) = units(indx_sup_lat(i)).pro.neural.base.rate_mu; end ;  for i = 1:length(indx_sup_lat), r_sup_anti_lat(i) = units(indx_sup_lat(i)).anti.neural.base.rate_mu; end
       
       
       
       % compare just baselines medial vs lateral
        exc_all_omv = [r_exc_pro_omv r_exc_anti_omv]; exc_all_omv_mu = mean(exc_all_omv); exc_all_omv_sem = (std(exc_all_omv))/(sqrt(length(exc_all_omv))); 
        sup_all_omv = [r_sup_pro_omv r_sup_anti_omv]; sup_all_omv_mu = mean(sup_all_omv); sup_all_omv_sem = (std(sup_all_omv))/(sqrt(length(sup_all_omv))); 
        
        exc_all_lat = [r_exc_pro_lat r_exc_anti_lat]; exc_all_lat_mu = mean(exc_all_lat); exc_all_lat_sem = (std(exc_all_lat))/(sqrt(length(exc_all_lat))); 
        sup_all_lat = [r_sup_pro_lat r_sup_anti_lat]; sup_all_lat_mu = mean(sup_all_lat); sup_all_lat_sem = (std(sup_all_lat))/(sqrt(length(sup_all_lat))); 
        
        % plot exc
        grp = [zeros(size(exc_all_omv)) ones(size(exc_all_lat))];
        figure; hold on;
        boxplot([exc_all_omv exc_all_lat], grp); 
        set(gca,'ylim',([0 120]), 'xTickLabel',  {'Medial','Lateral'});; axis square; title('Facilitation')
        
        
        % plot sup 
%         figure; hold on; boxplot(sup_all_omv); ylim([0 180])
%         figure; hold on; boxplot(sup_all_lat); ylim([0 180])
        grp = [zeros(size(sup_all_omv)) ones(size(sup_all_lat))];
        figure; hold on;
        boxplot([sup_all_omv sup_all_lat], grp); 
        set(gca,'ylim',([0 180]), 'xTickLabel',  {'Medial','Lateral'}); axis square; title('Suppression')
        
        ranksum(exc_all_omv, exc_all_lat)
        ranksum(sup_all_omv, sup_all_lat)
        ranksum([exc_all_omv' ; exc_all_lat'], [sup_all_omv' ; sup_all_lat'])
        ranksum([exc_all_omv' ; exc_all_lat'], [sup_all_omv' ; sup_all_lat'], 'Tail', 'left')
        % 
        ranksum(exc_all_omv, sup_all_omv, 'Tail', 'left')
        ranksum(exc_all_lat, sup_all_lat,'Tail', 'left')
        
        %% compare baseline firing rate vs max firing rate
        % extract max FR 
        for i = 1:length(indx_exc_omv), r_exc_pro_omv_peak(i,:) = units(indx_exc_omv(i)).pro.neural.sacc.peak_resp; end
        for i = 1:length(indx_exc_lat), r_exc_pro_lat_peak(i,:) = units(indx_exc_lat(i)).pro.neural.sacc.peak_resp; end
        for i = 1:length(indx_sup_omv), r_sup_pro_omv_peak(i,:) = units(indx_sup_omv(i)).pro.neural.sacc.peak_resp; end
        for i = 1:length(indx_sup_lat), r_sup_pro_lat_peak(i,:) = units(indx_sup_lat(i)).pro.neural.sacc.peak_resp; end
        
        % plot
        subplot(1,4,1); hold on; plot(r_exc_pro_omv,r_exc_pro_omv_peak, '.k', 'MarkerSize', 20); plot([0 160], [0 160]); set(gca, 'xlim', [0 160], 'ylim', [0 160], 'xTick',[0 160],'yTick',[0 160]); axis square; title('Fac Medial');
        subplot(1,4,2); hold on; plot(r_exc_pro_lat,r_exc_pro_lat_peak, '.k','MarkerSize', 20); plot([0 180], [0 180]); set(gca, 'xlim', [0 180], 'ylim', [0 180], 'xTick',[0 180],'yTick',[0 180]); axis square; title('Fac Lateral');
        subplot(1,4,3); hold on; plot(r_sup_pro_omv,r_sup_pro_omv_peak, '.k','MarkerSize', 20); plot([0 120], [0 120]); set(gca, 'xlim', [0 120], 'ylim', [0 120], 'xTick',[0 120],'yTick',[0 120]); axis square; title('Sup Medial');
        subplot(1,4,4); hold on; plot(r_sup_pro_lat,r_sup_pro_lat_peak, '.k','MarkerSize', 20); plot([0 180], [0 180]); set(gca, 'xlim', [0 180], 'ylim', [0 180], 'xTick',[0 180],'yTick',[0 180]); axis square; title('Sup Lateral');
        xlabel('Mean baseline firing (spks/s)'); ylabel('Peak firing rate (spsk/s)')
        
        
    case 'delta_rate_mean' % for exc and sup separately
        
        % get exc and sup
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc==1
                indx_exc(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup==1
                indx_sup(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        % get exc and sup signif
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc==1 && units(cellNum).stats.sacc.flags.proVsAnti_sacc_ks_nspk==1
                indx_exc_signif(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup==1 && units(cellNum).stats.sacc.flags.proVsAnti_sacc_ks_nspk==1
                indx_sup_signif(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        % get exc
        t = units(1).pro.neural.sacc.ts_pst_win;
        for i = 1:length(indx_exc)
            r_exc_pro(i,:) = units(indx_exc(i)).pro.neural.sacc.delta_rate_base;
            sem_exc_pro(i,:)= std(units(indx_exc(i)).pro.neural.sacc.delta_rate_base)/sqrt(length(indx_exc));
            r_exc_anti(i,:) = units(indx_exc(i)).anti.neural.sacc.delta_rate_base;
            sem_exc_anti(i,:)= std(units(indx_exc(i)).anti.neural.sacc.delta_rate_base)/sqrt(length(indx_sup));
        end
        
        % get exc signif
        for i = 1:length(indx_exc_signif)
            r_exc_pro_signif(i,:) = units(indx_exc_signif(i)).pro.neural.sacc.delta_rate_base;
            sem_exc_pro_signif(i,:)= std(units(indx_exc_signif(i)).pro.neural.sacc.delta_rate_base)/sqrt(length(indx_exc_signif));
            r_exc_anti_signif(i,:) = units(indx_exc_signif(i)).anti.neural.sacc.delta_rate_base;
            sem_exc_anti_signif(i,:)= std(units(indx_exc_signif(i)).anti.neural.sacc.delta_rate_base)/sqrt(length(indx_sup_signif));
        end
        
        % get sup
        for i = 1:length(indx_sup)
            r_sup_pro(i,:) = units(indx_sup(i)).pro.neural.sacc.delta_rate_base;
            sem_sup_pro(i,:)= std(units(indx_sup(i)).pro.neural.sacc.delta_rate_base)/sqrt(length(indx_sup));
            r_sup_anti(i,:) = units(indx_sup(i)).anti.neural.sacc.delta_rate_base;
            sem_sup_anti(i,:)= std(units(indx_sup(i)).anti.neural.sacc.delta_rate_base)/sqrt(length(indx_sup));
        end
        
        % get sup signif
        for i = 1:length(indx_sup_signif)
            r_sup_pro_signif(i,:) = units(indx_sup_signif(i)).pro.neural.sacc.delta_rate_base;
            sem_sup_pro_signif(i,:)= std(units(indx_sup_signif(i)).pro.neural.sacc.delta_rate_base)/sqrt(length(indx_sup_signif));
            r_sup_anti_signif(i,:) = units(indx_sup_signif(i)).anti.neural.sacc.delta_rate_base;
            sem_sup_anti_signif(i,:)= std(units(indx_sup_signif(i)).anti.neural.sacc.delta_rate_base)/sqrt(length(indx_sup_signif));
        end
        
        % plot exc
        figure; hold on;
        plot(t,mean(r_exc_pro));
        plot(t,mean(r_exc_anti));
        s_pro = shadedErrorBar(t, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','r');
        s_anti = shadedErrorBar(t, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','g');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        %set(gca, 'xlim',[-0.150 0.151], 'ylim', [-4 16], 'ytick', [-4 16],'TickDir', 'out', 'FontSize', 18);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [-10 15], 'ytick', [-10 0 15],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)')
        
        % plot sup
        %figure; hold on;
        plot(t,mean(r_sup_pro));
        plot(t,mean(r_sup_anti));
        s_pro = shadedErrorBar(t, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','r');
        s_anti = shadedErrorBar(t, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','g');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        %set(gca, 'xlim',[-0.150 0.151], 'ylim', [-10 2], 'ytick', [-10 2],'TickDir', 'out', 'FontSize', 18);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [-10 15], 'ytick', [-10 0 15],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)')
        
        % plot exc signif
        figure; hold on;
        plot(t,mean(r_exc_pro_signif));
        plot(t,mean(r_exc_anti_signif));
        s_pro = shadedErrorBar(t, mean(r_exc_pro_signif), repmat(mean(sem_exc_pro_signif),[size(mean(r_exc_pro_signif)) 1]), 'lineprops','r');
        s_anti = shadedErrorBar(t, mean(r_exc_anti_signif),repmat(mean(sem_exc_anti_signif),[size(mean(r_exc_anti_signif)) 1]), 'lineprops','g');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        % set(gca, 'xlim',[-0.150 0.151], 'ylim', [-10 30], 'ytick', [-10 30],'TickDir', 'out', 'FontSize', 18);
        % set(gca, 'xlim',[-0.150 0.151], 'ylim', [-5 20], 'ytick', [-5 20],'TickDir', 'out', 'FontSize', 18);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [-20 30], 'ytick', [-20 0 30],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)'); title('exc + signif ks ts')
        
        
        % plot sup signif
        %figure; hold on;
        plot(t,mean(r_sup_pro_signif));
        plot(t,mean(r_sup_anti_signif));
        s_pro = shadedErrorBar(t, mean(r_sup_pro_signif), repmat(mean(sem_sup_pro_signif),[size(mean(r_sup_pro_signif)) 1]), 'lineprops','r');
        s_anti = shadedErrorBar(t, mean(r_sup_anti_signif),repmat(mean(sem_sup_anti_signif),[size(mean(r_sup_anti_signif)) 1]), 'lineprops','g');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        %set(gca, 'xlim',[-0.150 0.151], 'ylim', [-20 5], 'ytick', [-20 5],'TickDir', 'out', 'FontSize', 18);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [-20 30], 'ytick', [-20 0 30],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)'); title('sup + signif ks')
        
    case 'delta_rate_mean_norm_pro'
        % get exc and sup
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc==1
                indx_exc(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup==1
                indx_sup(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        % get exc and sup signif
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc==1 && units(cellNum).stats.sacc.flags.proVsAnti_sacc_ks_nspk==1
                indx_exc_signif(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup==1 && units(cellNum).stats.sacc.flags.proVsAnti_sacc_ks_nspk==1
                indx_sup_signif(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        % get exc
        t = units(1).pro.neural.sacc.ts_pst_win;
        for i = 1:length(indx_exc)
            r_exc_pro(i,:) = units(indx_exc(i)).pro.neural.sacc.norm.delta_rate;
            std_exc_pro(i,:)= std(units(indx_exc(i)).pro.neural.sacc.norm.delta_rate);
            sem_exc_pro(i,:)= std(units(indx_exc(i)).pro.neural.sacc.norm.delta_rate)/sqrt(length(indx_exc));
            r_exc_anti(i,:) = units(indx_exc(i)).anti.neural.sacc.norm.delta_rate;
            std_exc_anti(i,:)= std(units(indx_exc(i)).anti.neural.sacc.norm.delta_rate);
            sem_exc_anti(i,:)= std(units(indx_exc(i)).anti.neural.sacc.norm.delta_rate)/sqrt(length(indx_exc));
        end
        
        % get exc signif
        for i = 1:length(indx_exc_signif)
            r_exc_pro_signif(i,:) = units(indx_exc_signif(i)).pro.neural.sacc.norm.delta_rate;
            std_exc_pro_signif(i,:)= std(units(indx_exc_signif(i)).pro.neural.sacc.norm.delta_rate);
            sem_exc_pro_signif(i,:)= std(units(indx_exc_signif(i)).pro.neural.sacc.norm.delta_rate)/sqrt(length(indx_exc_signif));
            r_exc_anti_signif(i,:) = units(indx_exc_signif(i)).anti.neural.sacc.norm.delta_rate;
            std_exc_anti_signif(i,:)= std(units(indx_exc_signif(i)).anti.neural.sacc.norm.delta_rate);
            sem_exc_anti_signif(i,:)= std(units(indx_exc_signif(i)).anti.neural.sacc.norm.delta_rate)/sqrt(length(indx_exc_signif));
        end
        
        % get sup
        for i = 1:length(indx_sup)
            r_sup_pro(i,:) = units(indx_sup(i)).pro.neural.sacc.norm.delta_rate;
            std_sup_pro(i,:)= std(units(indx_sup(i)).pro.neural.sacc.norm.delta_rate);
            sem_sup_pro(i,:)= std(units(indx_sup(i)).pro.neural.sacc.norm.delta_rate)/sqrt(length(indx_sup));
            r_sup_anti(i,:) = units(indx_sup(i)).anti.neural.sacc.norm.delta_rate;
            std_sup_anti(i,:)= std(units(indx_sup(i)).anti.neural.sacc.norm.delta_rate);
            sem_sup_anti(i,:)= std(units(indx_sup(i)).anti.neural.sacc.norm.delta_rate)/sqrt(length(indx_sup));
        end
        
        % get exc signif
        for i = 1:length(indx_sup_signif)
            r_sup_pro_signif(i,:) = units(indx_sup_signif(i)).pro.neural.sacc.norm.delta_rate;
            std_sup_pro_signif(i,:)= std(units(indx_sup_signif(i)).pro.neural.sacc.norm.delta_rate);
            sem_sup_pro_signif(i,:)= std(units(indx_sup_signif(i)).pro.neural.sacc.norm.delta_rate)/sqrt(length(indx_sup_signif));
            r_sup_anti_signif(i,:) = units(indx_sup_signif(i)).anti.neural.sacc.norm.delta_rate;
            std_sup_anti_signif(i,:)= std(units(indx_sup_signif(i)).anti.neural.sacc.norm.delta_rate);
            sem_sup_anti_signif(i,:)= std(units(indx_sup_signif(i)).anti.neural.sacc.norm.delta_rate)/sqrt(length(indx_sup_signif));
        end
        
        % plot exc
        figure; hold on;
        plot(t,mean(r_exc_pro));
        plot(t,mean(r_exc_anti));
        s_pro = shadedErrorBar(t, mean(r_exc_pro),repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        %set(gca, 'xlim',[-0.150 0.151], 'ylim', [-4 16], 'ytick', [-4 16],'TickDir', 'out', 'FontSize', 18);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [-2 2], 'ytick', [-2 2],'TickDir', 'out', 'FontSize', 26);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)')
        
        % plot sup
        %figure; hold on;
        plot(t,mean(r_sup_pro));
        plot(t,mean(r_sup_anti));
        s_pro = shadedErrorBar(t, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        %set(gca, 'xlim',[-0.150 0.151], 'ylim', [-10 2], 'ytick', [-10 2],'TickDir', 'out', 'FontSize', 18);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [-2.7 2.7], 'ytick', [-2.7 0 2.7],'TickDir', 'out', 'FontSize', 26);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)')
        
        % plot exc signif
        figure; hold on;
        plot(t,mean(r_exc_pro_signif));
        plot(t,mean(r_exc_anti_signif));
        s_pro = shadedErrorBar(t, mean(r_exc_pro_signif), repmat(mean(sem_exc_pro_signif),[size(mean(r_exc_pro_signif)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t, mean(r_exc_anti_signif),repmat(mean(sem_exc_anti_signif),[size(mean(r_exc_anti_signif)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [-2.5 2.5], 'ytick', [-2.5 0 2.5],'TickDir', 'out', 'FontSize', 26);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)');
        
        
        % plot sup signif
        %figure; hold on;
        plot(t,mean(r_sup_pro_signif));
        plot(t,mean(r_sup_anti_signif));
        s_pro = shadedErrorBar(t, mean(r_sup_pro_signif), repmat(mean(sem_sup_pro_signif),[size(mean(r_sup_pro_signif)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t, mean(r_sup_anti_signif),repmat(mean(sem_sup_anti_signif),[size(mean(r_sup_anti_signif)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.150 0.151], 'ylim', [-2.5 2.5], 'ytick', [-2.5 0 2.5],'TickDir', 'out', 'FontSize', 26);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)'); title('signif')
        
        
        %% insets with means for smaller windows pre saccade
        % pro small_win -0.1 to 0 ::: win -0.150 to 0
        %         r_small_win_exc_pro = mean(r_exc_pro_signif(:, t > -0.101 & t < 0.01));
        r_win_exc_pro = mean(r_exc_pro_signif(:, t > -0.151 & t < 0.01));
        %         r_small_win_sup_pro = mean(r_sup_pro_signif(:, t > -0.101 & t < 0.01));
        r_win_sup_pro = mean(r_sup_pro_signif(:, t > -0.151 & t < 0.01));
        
        % anti small_win -0.1 to 0 ::: win -0.150 to 0
        %         r_small_win_exc_anti = mean(r_exc_anti_signif(:, t > -0.101 & t < 0.01));
        r_win_exc_anti = mean(r_exc_anti_signif(:, t > -0.151 & t < 0.01));
        %         r_small_win_sup_anti = mean(r_sup_anti_signif(:, t > -0.101 & t < 0.01));
        r_win_sup_anti = mean(r_sup_anti_signif(:, t > -0.151 & t < 0.01));
        
        %         figure; hold on;
        %         errorbar(1,mean(r_small_win_exc_pro),std(r_small_win_exc_pro),'m','LineWidth',1, 'Marker', 'o');
        %         errorbar(2,mean(r_small_win_exc_anti),std(r_small_win_exc_anti),'b','LineWidth',1,'Marker', 'o');
        %         set(gca,'xlim', [0 3], 'xTick',[], 'ylim',[-0.2 1.2], 'yTick', [0 0.5 1], 'TickDir', 'out', 'FontSize',30);
        %         title('-0.1 to 0')
        
        figure; hold on;
        errorbar(1,mean(r_win_exc_pro),std(r_win_exc_pro),'g','LineWidth',1, 'Marker', '^','Capsize',0);
        errorbar(2,mean(r_win_exc_anti),std(r_win_exc_anti),'b','LineWidth',1, 'Marker', '^','Capsize',0);
        errorbar(1,mean(r_win_sup_pro),std(r_win_sup_pro),'m','LineWidth',1, 'Marker', 'v','Capsize',0);
        errorbar(2,mean(r_win_sup_anti),std(r_win_sup_anti),'b','LineWidth',1, 'Marker', 'v','Capsize',0);
        
        set(gca,'xlim', [0 3], 'xTick',[], 'ylim',[-1.2 1.2], 'yTick', [-1 -0.5 0 0.5 1], 'TickDir', 'out', 'FontSize',30);
        title('-0.150 to 0')
        
        %% insets with means for smaller windows post saccade
        % pro small_win 0 to 0.100 ::: win 1 to 0.150
        %         r_small_win_exc_pro = mean(r_exc_pro_signif(:, t > -0.01 & t < 0.1));
        r_win_exc_pro = mean(r_exc_pro_signif(:, t > -0.01 & t < 0.150));
        %         r_small_win_sup_pro = mean(r_sup_pro_signif(:, t > -0.01 & t < 0.1));
        r_win_sup_pro = mean(r_sup_pro_signif(:, t > -0.01 & t < 0.150));
        
        % anti small_win -0.1 to 0 ::: win -0.150 to 0
        %         r_small_win_exc_anti = mean(r_exc_anti_signif(:, t > -0.01 & t < 0.1));
        r_win_exc_anti = mean(r_exc_anti_signif(:, t > -0.01 & t < 0.150));
        %         r_small_win_sup_anti = mean(r_sup_anti_signif(:, t > -0.01 & t < 0.1));
        r_win_sup_anti = mean(r_sup_anti_signif(:, t > -0.01 & t < 0.150));
        
        %%% take grand mean
        
        %         figure; hold on;
        %         errorbar(1,mean(r_small_win_exc_pro),std(r_small_win_exc_pro),'m','LineWidth',1, 'Marker', 'o');
        %         errorbar(2,mean(r_small_win_exc_anti),std(r_small_win_exc_anti),'b','LineWidth',1,'Marker', 'o');
        %         set(gca,'xlim', [0 3], 'xTick',[], 'ylim',[0 2], 'yTick', [0 1 2], 'TickDir', 'out', 'FontSize',30);
        %         title('0 to 0.1')
        
        figure; hold on;
        errorbar(1,mean(r_win_exc_pro),std(r_win_exc_pro),'m','LineWidth',1, 'Marker', '^', 'Capsize',0);
        errorbar(2,mean(r_win_exc_anti),std(r_win_exc_anti),'b','LineWidth',1, 'Marker', '^','Capsize',0);
        errorbar(1,mean(r_win_sup_pro),std(r_win_sup_pro),'m','LineWidth',1, 'Marker', 'v','Capsize',0);
        errorbar(2,mean(r_win_sup_anti),std(r_win_sup_anti),'b','LineWidth',1, 'Marker', 'v','Capsize',0);
        
        set(gca,'xlim', [0 3], 'xTick',[], 'ylim',[-2 2], 'yTick', [-2 0 2], 'TickDir', 'out', 'FontSize',30);
        title('0 to 0.150')
        
        % stats : ranksum
        [h,p] = ranksum(r_win_exc_pro, r_win_exc_anti)
        
        
        
    case 'delta_rate_mean_norm_instr'
        % get exc and sup
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc_instr==1
                indx_exc(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup_instr==1
                indx_sup(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        % get exc and sup signif
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc_instr==1 && units(cellNum).stats.instr.flags.proVsAnti_instr_ks_nspk==1
                indx_exc_signif(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup_instr==1 && units(cellNum).stats.instr.flags.proVsAnti_instr_ks_nspk==1
                indx_sup_signif(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        % get exc
        t = units(1).pro.neural.instr.ts_pst_win;
        for i = 1:length(indx_exc)
            r_exc_pro(i,:) = units(indx_exc(i)).pro.neural.instr.norm.delta_rate;
            std_exc_pro(i,:)= std(units(indx_exc(i)).pro.neural.instr.norm.delta_rate);
            sem_exc_pro(i,:)= std(units(indx_exc(i)).pro.neural.instr.norm.delta_rate)/sqrt(length(indx_exc));
            r_exc_anti(i,:) = units(indx_exc(i)).anti.neural.instr.norm.delta_rate;
            std_exc_anti(i,:)= std(units(indx_exc(i)).anti.neural.instr.norm.delta_rate);
            sem_exc_anti(i,:)= std(units(indx_exc(i)).anti.neural.instr.norm.delta_rate)/sqrt(length(indx_exc));
        end
        
        % get exc signif
        for i = 1:length(indx_exc_signif)
            r_exc_pro_signif(i,:) = units(indx_exc_signif(i)).pro.neural.instr.norm.delta_rate;
            std_exc_pro_signif(i,:)= std(units(indx_exc_signif(i)).pro.neural.instr.norm.delta_rate);
            sem_exc_pro_signif(i,:)= std(units(indx_exc_signif(i)).pro.neural.instr.norm.delta_rate)/sqrt(length(indx_exc_signif));
            r_exc_anti_signif(i,:) = units(indx_exc_signif(i)).anti.neural.instr.norm.delta_rate;
            std_exc_anti_signif(i,:)= std(units(indx_exc_signif(i)).anti.neural.instr.norm.delta_rate);
            sem_exc_anti_signif(i,:)= std(units(indx_exc_signif(i)).anti.neural.instr.norm.delta_rate)/sqrt(length(indx_exc_signif));
        end
        
        % get sup
        for i = 1:length(indx_sup)
            r_sup_pro(i,:) = units(indx_sup(i)).pro.neural.instr.norm.delta_rate;
            std_sup_pro(i,:)= std(units(indx_sup(i)).pro.neural.instr.norm.delta_rate);
            sem_sup_pro(i,:)= std(units(indx_sup(i)).pro.neural.instr.norm.delta_rate)/sqrt(length(indx_sup));
            r_sup_anti(i,:) = units(indx_sup(i)).anti.neural.instr.norm.delta_rate;
            std_sup_anti(i,:)= std(units(indx_sup(i)).anti.neural.instr.norm.delta_rate);
            sem_sup_anti(i,:)= std(units(indx_sup(i)).anti.neural.instr.norm.delta_rate)/sqrt(length(indx_sup));
        end
        
        % get exc signif
        for i = 1:length(indx_sup_signif)
            r_sup_pro_signif(i,:) = units(indx_sup_signif(i)).pro.neural.instr.norm.delta_rate;
            std_sup_pro_signif(i,:)= std(units(indx_sup_signif(i)).pro.neural.instr.norm.delta_rate);
            sem_sup_pro_signif(i,:)= std(units(indx_sup_signif(i)).pro.neural.instr.norm.delta_rate)/sqrt(length(indx_sup_signif));
            r_sup_anti_signif(i,:) = units(indx_sup_signif(i)).anti.neural.instr.norm.delta_rate;
            std_sup_anti_signif(i,:)= std(units(indx_sup_signif(i)).anti.neural.instr.norm.delta_rate);
            sem_sup_anti_signif(i,:)= std(units(indx_sup_signif(i)).anti.neural.instr.norm.delta_rate)/sqrt(length(indx_sup_signif));
        end
        
        % plot exc
        figure; hold on;
        %         plot(t,mean(r_exc_pro));
        %         plot(t,mean(r_exc_anti));
        s_pro = shadedErrorBar(t, mean(r_exc_pro),repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[0.05 0.350], 'ylim', [0 1], 'ytick', [0 1],'TickDir', 'out', 'FontSize', 26 );
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)');hline(0, '--k');
        
        % plot sup
        %figure; hold on;
        %         plot(t,mean(r_sup_pro));
        %         plot(t,mean(r_sup_anti));
        s_pro = shadedErrorBar(t, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        figure; hold on;
        plot(t,mean(r_exc_pro_signif)); % plot(t,r_exc_pro_signif); %
        plot(t,mean(r_exc_anti_signif)); % plot(t,r_exc_anti_signif); %
        s_pro = shadedErrorBar(t, mean(r_exc_pro_signif), repmat(mean(sem_exc_pro_signif),[size(mean(r_exc_pro_signif)) 1]), 'lineprops','m'); % s_pro = shadedErrorBar(t, r_exc_pro_signif, repmat(mean(sem_exc_pro_signif),[size(r_exc_pro_signif) 1]), 'lineprops','r'); %
        s_anti = shadedErrorBar(t, mean(r_exc_anti_signif),repmat(mean(sem_exc_anti_signif),[size(mean(r_exc_anti_signif)) 1]), 'lineprops','b'); % s_anti = shadedErrorBar(t, r_exc_anti_signif,repmat(mean(sem_exc_anti_signif),[size(r_exc_anti_signif) 1]), 'lineprops','g');%
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[0.05 0.350], 'ylim', [-0.5 2.5], 'ytick', [-0.5 0 5],'TickDir', 'out', 'FontSize', 26);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)'); hline(0, '--k');
        
        
        % plot sup signif
        %figure; hold on;
        plot(t,mean(r_sup_pro_signif)); % plot(t,r_sup_pro_signif);%
        plot(t,mean(r_sup_anti_signif)); % plot(t,r_sup_anti_signif); %
        s_pro = shadedErrorBar(t, mean(r_sup_pro_signif), repmat(mean(sem_sup_pro_signif),[size(mean(r_sup_pro_signif)) 1]), 'lineprops','m'); % s_pro = shadedErrorBar(t, r_sup_pro_signif, repmat(mean(sem_sup_pro_signif),[size(r_sup_pro_signif) 1]), 'lineprops','r');%
        s_anti = shadedErrorBar(t, mean(r_sup_anti_signif),repmat(mean(sem_sup_anti_signif),[size(mean(r_sup_anti_signif)) 1]), 'lineprops','b'); % s_anti = shadedErrorBar(t, r_sup_anti_signif,repmat(mean(sem_sup_anti_signif),[size(r_sup_anti_signif) 1]), 'lineprops','g'); %
        set(s_pro.mainLine,'LineWidth', 4), set(s_anti.mainLine,'LineWidth', 4);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[0.05 0.350], 'ylim', [-2.5 2.5], 'ytick', [-2.5 0 2.5],'TickDir', 'out', 'FontSize', 26);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)'); title('signif');hline(0, '--k');
        
    case 'change_in_fr_anti_pro_selected'
        indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=0 & t<=0.150);
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.sacc.rate_pst_win - mean(units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0));...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro(i) = max(r_exc_pro(i,:)); max_sem_exc_pro(i) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.sacc.rate_pst_win - mean(units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro(i) = max(abs(r_sup_pro(i,:))); max_sem_sup_pro(i) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.sacc.rate_pst_win - mean(units(indx_exc_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti(i) = max(r_exc_anti(i,:)); max_sem_exc_anti(i) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti(i,:) =  units(indx_sup_anti(i)).anti.neural.sacc.rate_pst_win - mean(units(indx_sup_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti(i) = max(abs(r_sup_anti(i,:))); max_sem_sup_anti(i) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        % plot exc
        figure; hold on;
        plot(t_win,mean(r_exc_pro));
        plot(t_win,mean(r_exc_anti));
        s_pro = shadedErrorBar(t_win, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[0 0.151], 'ylim', [-20 35], 'ytick', [-20 0 35],'TickDir', 'out', 'FontSize', 18);
        
        % plot sup
        plot(t_win,mean(r_sup_pro));
        plot(t_win,mean(r_sup_anti));
        s_pro = shadedErrorBar(t_win, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[0 0.151], 'ylim', [-20 35], 'ytick', [-20 0 35],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)'); hline(0,'--k'); axis square
        
        % gather max modulation and average
        figure; hold on;
        errorbar(1,mean(max_exc_pro), mean(max_sem_exc_pro),'.m', 'MarkerSize', 18, 'CapSize', 0);
        errorbar(1,mean(max_sup_pro), mean(max_sem_exc_pro),'om', 'MarkerSize', 18,'CapSize', 0);
        errorbar(1,mean(max_exc_anti), mean(max_sem_exc_pro),'.b', 'MarkerSize', 18,'CapSize', 0);
        errorbar(1,mean(max_sup_anti), mean(max_sem_exc_pro),'ob', 'MarkerSize', 18,'CapSize', 0);
        set(gca,'xlim', [0 2], 'ylim', [0 40], 'ytick', [0 20 40],'TickDir', 'out', 'FontSize', 18);
        ylabel('Max abs change in firing rate (spk/s)'); axis square
        
         % anti - pro
        figure; hold on; 
        plot(t_win,(mean(r_exc_anti))-(mean(r_exc_pro)), 'k','LineWidth',2);
        plot(t_win,(mean(r_sup_anti))-(mean(r_sup_pro)), '--k','LineWidth',2);
        set(gca,'xlim', [0 0.150],'xTick',[0 0.05 0.1 0.15], 'ylim', [-8 6], 'ytick', [-8 0 6],'TickDir', 'out', 'FontSize', 18);
        ylabel('Change in firing rate (spk/s)'); xlabel('Time (s)')
        axis square; hline(0,'k'); title(recArea); 
        
    case 'change_in_fr_anti_pro_instr_selected'
        indx_exc_pro = pop.indx_sel.(recArea).instr_back.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).instr_back.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).instr_back.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).instr_back.all.anti.sup;
        
        t_win = units(1).pro.neural.instr_back.ts_pst_win;
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win - units(indx_exc_pro(i)).pro.neural.base.rate_mu; ...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro(i) = max(r_exc_pro(i,:)); max_sem_exc_pro(i) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win - units(indx_sup_pro(i)).pro.neural.base.rate_mu; ...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro(i) = max(abs(r_sup_pro(i,:))); max_sem_sup_pro(i) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.instr_back.rate_pst_win - units(indx_exc_anti(i)).anti.neural.base.rate_mu; ...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti(i) = max(r_exc_anti(i,:)); max_sem_exc_anti(i) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti(i,:) =  units(indx_sup_anti(i)).anti.neural.instr_back.rate_pst_win - units(indx_sup_anti(i)).anti.neural.base.rate_mu; ...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti(i) = max(abs(r_sup_anti(i,:))); max_sem_sup_anti(i) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        % plot exc
        figure; hold on; 
        plot(t_win,mean(r_exc_pro));
        plot(t_win,mean(r_exc_anti));
        s_pro = shadedErrorBar(t_win, mean(r_exc_pro), repmat(mean(sem_exc_pro),[size(mean(r_exc_pro)) 1]), 'lineprops','m');
        s_anti = shadedErrorBar(t_win, mean(r_exc_anti),repmat(mean(sem_exc_anti),[size(mean(r_exc_anti)) 1]), 'lineprops','b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.3 0.1], 'ylim', [-12 12], 'ytick', [-6 0 12],'TickDir', 'out', 'FontSize', 18);
        
        % plot sup
        plot(t_win,mean(r_sup_pro));
        plot(t_win,mean(r_sup_anti));
        s_pro = shadedErrorBar(t_win, mean(r_sup_pro), repmat(mean(sem_sup_pro),[size(mean(r_sup_pro)) 1]), 'lineprops','--m');
        s_anti = shadedErrorBar(t_win, mean(r_sup_anti),repmat(mean(sem_sup_anti),[size(mean(r_sup_anti)) 1]), 'lineprops','--b');
        set(s_pro.mainLine,'LineWidth', 2), set(s_anti.mainLine,'LineWidth', 2);
        set(s_pro.edge,'LineStyle', 'none'); set(s_anti.edge,'LineStyle', 'none');
        set(s_pro.patch, 'FaceAlpha', 0.1); set(s_anti.patch, 'FaceAlpha', 0.1);
        set(gca, 'xlim',[-0.3 0.1], 'ylim', [-12 12], 'ytick', [-12 0 12],'TickDir', 'out', 'FontSize', 18);
        ylabel ('Change in firing rate (spk/s)'); xlabel('Time (s)'); hline(0,'-k'); axis square
        
         % anti - pro
        figure; hold on; 
        plot(t_win,(mean(r_exc_anti))-(mean(r_exc_pro)), 'k','LineWidth',2);
        plot(t_win,(mean(r_sup_anti))-(mean(r_sup_pro)), '--k','LineWidth',2);
        set(gca,'xlim', [-0.3 0.1],'xTick',[-0.3 -0.2 -0.1 0], 'ylim', [-15 20], 'ytick', [-15 0 20],'TickDir', 'out', 'FontSize', 18);
        ylabel('Change in firing rate (spk/s)'); xlabel('Time (s)')
        axis square; hline(0,'k'); title(recArea); 
        
        % anti - pro plotted with abs values for both
        figure; hold on; 
        plot(t_win,(mean(r_exc_anti))-(mean(r_exc_pro)), 'k','LineWidth',2);
        plot(t_win,(mean(abs(r_sup_anti))-(mean(abs(r_sup_pro)))), '--k','LineWidth',2);
        set(gca,'xlim', [-0.3 0],'xTick',[-0.3 -0.2 -0.1 0], 'ylim', [-10 20], 'ytick', [-10 0 20],'TickDir', 'out', 'FontSize', 18);
        ylabel('Abs change in firing rate (spk/s)'); xlabel('Time (s)')
        axis square; hline(0,'k'); title(recArea);
        
        % test if change in fr is higher in facilitation than in suppression
        [h,p] = ranksum(mean(r_exc_anti)-mean(r_exc_pro),mean(abs(r_sup_anti))-(mean(abs(r_sup_pro))))
        
        early_instr_exc = mean(r_exc_anti(:,t_win>-0.31 & t_win<-0.19))-mean(r_exc_pro(:,t_win>-0.31 & t_win<-0.19));
        late_instr_exc = mean(r_exc_anti(:,t_win>-0.1 & t_win<=0))-mean(r_exc_pro(:,t_win>-0.1 & t_win<=0));
        
        early_instr_sup = mean(r_sup_anti(:,t_win>-0.31 & t_win<-0.19))-mean(r_sup_pro(:,t_win>-0.31 & t_win<-0.19));
        late_instr_sup = mean(abs(r_sup_anti(:,t_win>-0.1 & t_win<=0)))-mean(abs(r_sup_pro(:,t_win>-0.1 & t_win<=0)));
        
        %% test if fac is diff from supp at the end 
        exc_ch = (mean(r_exc_anti))-(mean(r_exc_pro)); 
        exc_ch_end = exc_ch(21:31); 
        
        sup_ch = (mean(abs(r_sup_anti))-(mean(abs(r_sup_pro)))); 
        sup_ch_end = sup_ch(21:31); 
        
        kstest2(exc_ch_end,sup_ch_end)
        
    case 'max_abs_change_fr'
        % vermis
        indx_exc_pro = pop.indx_sel.vermis.sacc.all.pro.exc;
        indx_sup_pro = pop.indx_sel.vermis.sacc.all.pro.sup;
        indx_exc_anti = pop.indx_sel.vermis.sacc.all.anti.exc;
        indx_sup_anti = pop.indx_sel.vermis.sacc.all.anti.sup;
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=0 & t<=0.150);
        for i = 1:length(indx_exc_pro), r_exc_pro_omv(i,:) = units(indx_exc_pro(i)).pro.neural.sacc.rate_pst_win - mean(units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0));...
                sem_exc_pro(i,:) = std(r_exc_pro_omv(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro_omv(i) = max(r_exc_pro_omv(i,:)); max_sem_exc_pro_omv(i) = std(r_exc_pro_omv(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro_omv(i,:) = units(indx_sup_pro(i)).pro.neural.sacc.rate_pst_win - mean(units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_sup_pro(i,:) = std(r_sup_pro_omv(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro_omv(i) = max(abs(r_sup_pro_omv(i,:))); max_sem_sup_pro_omv(i) = std(r_sup_pro_omv(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti_omv(i,:) = units(indx_exc_anti(i)).anti.neural.sacc.rate_pst_win - mean(units(indx_exc_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_exc_anti(i,:) = std(r_exc_anti_omv(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti_omv(i) = max(r_exc_anti_omv(i,:)); max_sem_exc_anti_omv(i) = std(r_exc_anti_omv(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti_omv(i,:) =  units(indx_sup_anti(i)).anti.neural.sacc.rate_pst_win - mean(units(indx_sup_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_sup_anti(i,:) = std(r_sup_anti_omv(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti_omv(i) = max(abs(r_sup_anti_omv(i,:))); max_sem_sup_anti_omv(i) = std(r_sup_anti_omv(i,:))/sqrt(length(indx_sup_anti)); end
        
        clear indx_exc_pro indx_sup_pro indx_exc_anti indx_sup_anti
        indx_exc_pro = pop.indx_sel.lateral.sacc.all.pro.exc;
        indx_sup_pro = pop.indx_sel.lateral.sacc.all.pro.sup;
        indx_exc_anti = pop.indx_sel.lateral.sacc.all.anti.exc;
        indx_sup_anti = pop.indx_sel.lateral.sacc.all.anti.sup;
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=0 & t<=0.150);
        for i = 1:length(indx_exc_pro), r_exc_pro_lat(i,:) = units(indx_exc_pro(i)).pro.neural.sacc.rate_pst_win - mean(units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0));...
                sem_exc_pro(i,:) = std(r_exc_pro_lat(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro_lat(i) = max(r_exc_pro_lat(i,:)); max_sem_exc_pro_lat(i) = std(r_exc_pro_lat(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro_lat(i,:) = units(indx_sup_pro(i)).pro.neural.sacc.rate_pst_win - mean(units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_sup_pro(i,:) = std(r_sup_pro_lat(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro_lat(i) = max(abs(r_sup_pro_lat(i,:))); max_sem_sup_pro_lat(i) = std(r_sup_pro_lat(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti_lat(i,:) = units(indx_exc_anti(i)).anti.neural.sacc.rate_pst_win - mean(units(indx_exc_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_exc_anti(i,:) = std(r_exc_anti_lat(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti_lat(i) = max(r_exc_anti_lat(i,:)); max_sem_exc_anti_lat(i) = std(r_exc_anti_lat(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti_lat(i,:) =  units(indx_sup_anti(i)).anti.neural.sacc.rate_pst_win - mean(units(indx_sup_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_sup_anti(i,:) = std(r_sup_anti_lat(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti_lat(i) = max(abs(r_sup_anti_lat(i,:))); max_sem_sup_anti_lat(i) = std(r_sup_anti_lat(i,:))/sqrt(length(indx_sup_anti)); end
        
        figure; hold on;
        errorbar(1,mean(max_exc_pro_omv), mean(max_sem_exc_pro_omv),'.m', 'MarkerSize', 18, 'CapSize', 0);
        errorbar(1,mean(max_sup_pro_omv), mean(max_sem_sup_pro_omv),'om', 'MarkerSize', 18,'CapSize', 0);
        errorbar(1,mean(max_exc_anti_omv), mean(max_sem_exc_anti_omv),'.b', 'MarkerSize', 18,'CapSize', 0);
        errorbar(1,mean(max_sup_anti_omv), mean(max_sem_sup_anti_omv),'ob', 'MarkerSize', 18,'CapSize', 0); 
        
        errorbar(2,mean(max_exc_pro_lat), mean(max_sem_exc_pro_lat),'.m', 'MarkerSize', 18, 'CapSize', 0);
        errorbar(2,mean(max_sup_pro_lat), mean(max_sem_sup_pro_lat),'om', 'MarkerSize', 18,'CapSize', 0);
        errorbar(2,mean(max_exc_anti_lat), mean(max_sem_exc_anti_lat),'.b', 'MarkerSize', 18,'CapSize', 0);
        errorbar(2,mean(max_sup_anti_lat), mean(max_sem_sup_anti_lat),'ob', 'MarkerSize', 18,'CapSize', 0);
        set(gca,'xlim', [0 2],'xTick',[], 'ylim', [0 40], 'ytick', [0 20 40],'TickDir', 'out', 'FontSize', 18);
        ylabel('Max abs change in firing rate (spk/s)'); axis square; title('Saccade')
        
        exc_pro_pval = ranksum(max_exc_pro_omv,max_exc_pro_lat)
        sup_pro_pval = ranksum(max_sup_pro_omv,max_sup_pro_lat)
        exc_anti_pval = ranksum(max_exc_anti_omv,max_exc_anti_lat)
        sup_anti_pval = ranksum(max_sup_anti_omv,max_sup_anti_lat)
        
        
    case 'max_abs_change_fr_instr'
        % vermis
        indx_exc_pro = pop.indx_sel.vermis.instr_back.all.pro.exc;
        indx_sup_pro = pop.indx_sel.vermis.instr_back.all.pro.sup;
        indx_exc_anti = pop.indx_sel.vermis.instr_back.all.anti.exc;
        indx_sup_anti = pop.indx_sel.vermis.instr_back.all.anti.sup;
        
        t_win = units(1).pro.neural.instr_back.ts_pst_win;
        for i = 1:length(indx_exc_pro), r_exc_pro_omv(i,:) = units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win - units(indx_exc_pro(i)).pro.neural.base.rate_mu; ...
                sem_exc_pro(i,:) = std(r_exc_pro_omv(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro_omv(i) = max(r_exc_pro_omv(i,:)); max_sem_exc_pro_omv(i) = std(r_exc_pro_omv(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro_omv(i,:) = units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win - units(indx_sup_pro(i)).pro.neural.base.rate_mu; ...
                sem_sup_pro(i,:) = std(r_sup_pro_omv(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro_omv(i) = max(abs(r_sup_pro_omv(i,:))); max_sem_sup_pro_omv(i) = std(r_sup_pro_omv(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti_omv(i,:) = units(indx_exc_anti(i)).anti.neural.instr_back.rate_pst_win - units(indx_exc_anti(i)).anti.neural.base.rate_mu; ...
                sem_exc_anti(i,:) = std(r_exc_anti_omv(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti_omv(i) = max(r_exc_anti_omv(i,:)); max_sem_exc_anti_omv(i) = std(r_exc_anti_omv(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti_omv(i,:) =  units(indx_sup_anti(i)).anti.neural.instr_back.rate_pst_win - units(indx_sup_anti(i)).anti.neural.base.rate_mu; ...
                sem_sup_anti(i,:) = std(r_sup_anti_omv(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti_omv(i) = max(abs(r_sup_anti_omv(i,:))); max_sem_sup_anti_omv(i) = std(r_sup_anti_omv(i,:))/sqrt(length(indx_sup_anti)); end
        
        % lateral
        clear indx_exc_pro indx_sup_pro indx_exc_anti indx_sup_anti
        indx_exc_pro = pop.indx_sel.lateral.instr_back.all.pro.exc;
        indx_sup_pro = pop.indx_sel.lateral.instr_back.all.pro.sup;
        indx_exc_anti = pop.indx_sel.lateral.instr_back.all.anti.exc;
        indx_sup_anti = pop.indx_sel.lateral.instr_back.all.anti.sup;
        
        t_win = units(1).pro.neural.instr_back.ts_pst_win;
        for i = 1:length(indx_exc_pro), r_exc_pro_lat(i,:) = units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win - units(indx_exc_pro(i)).pro.neural.base.rate_mu; ...
                sem_exc_pro(i,:) = std(r_exc_pro_lat(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro_lat(i) = max(r_exc_pro_lat(i,:)); max_sem_exc_pro_lat(i) = std(r_exc_pro_lat(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro_lat(i,:) = units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win - units(indx_sup_pro(i)).pro.neural.base.rate_mu; ...
                sem_sup_pro(i,:) = std(r_sup_pro_lat(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro_lat(i) = max(abs(r_sup_pro_lat(i,:))); max_sem_sup_pro_lat(i) = std(r_sup_pro_lat(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti_lat(i,:) = units(indx_exc_anti(i)).anti.neural.instr_back.rate_pst_win - units(indx_exc_anti(i)).anti.neural.base.rate_mu; ...
                sem_exc_anti(i,:) = std(r_exc_anti_lat(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti_lat(i) = max(r_exc_anti_lat(i,:)); max_sem_exc_anti_lat(i) = std(r_exc_anti_lat(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti_lat(i,:) =  units(indx_sup_anti(i)).anti.neural.instr_back.rate_pst_win - units(indx_sup_anti(i)).anti.neural.base.rate_mu; ...
                sem_sup_anti(i,:) = std(r_sup_anti_lat(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti_lat(i) = max(abs(r_sup_anti_lat(i,:))); max_sem_sup_anti_lat(i) = std(r_sup_anti_lat(i,:))/sqrt(length(indx_sup_anti)); end
        
        
        % gather max modulation and average
        figure; hold on;
        errorbar(1,mean(max_exc_pro_omv), mean(max_sem_exc_pro_omv),'.m', 'MarkerSize', 18, 'CapSize', 0);
        errorbar(1,mean(max_sup_pro_omv), mean(max_sem_sup_pro_omv),'om', 'MarkerSize', 18,'CapSize', 0);
        errorbar(1,mean(max_exc_anti_omv), mean(max_sem_exc_anti_omv),'.b', 'MarkerSize', 18,'CapSize', 0);
        errorbar(1,mean(max_sup_anti_omv), mean(max_sem_sup_anti_omv),'ob', 'MarkerSize', 18,'CapSize', 0); 
        
        errorbar(2,mean(max_exc_pro_lat), mean(max_sem_exc_pro_lat),'.m', 'MarkerSize', 18, 'CapSize', 0);
        errorbar(2,mean(max_sup_pro_lat), mean(max_sem_sup_pro_lat),'om', 'MarkerSize', 18,'CapSize', 0);
        errorbar(2,mean(max_exc_anti_lat), mean(max_sem_exc_anti_lat),'.b', 'MarkerSize', 18,'CapSize', 0);
        errorbar(2,mean(max_sup_anti_lat), mean(max_sem_sup_anti_lat),'ob', 'MarkerSize', 18,'CapSize', 0);
        set(gca,'xlim', [0 2],'xTick',[], 'ylim', [0 20.5], 'ytick', [0 10 20],'TickDir', 'out', 'FontSize', 18);
        ylabel('Max abs change in firing rate (spk/s)'); axis square; title('Instruction')
        
        exc_pro_pval = ranksum(max_exc_pro_omv,max_exc_pro_lat)
        sup_pro_pval = ranksum(max_sup_pro_omv,max_sup_pro_lat)
        exc_anti_pval = ranksum(max_exc_anti_omv,max_exc_anti_lat)
        sup_anti_pval = ranksum(max_sup_anti_omv,max_sup_anti_lat)
        
      case 'fr_comparison_instr_sacc_selected'
        indx_exc_pro = pop.indx_sel.(recArea).instr_back.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).instr_back.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).instr_back.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).instr_back.all.anti.sup;
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=0 & t<=0.150);
        % sacc .... psth during sacc window - mean of activity before sacc onset 
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.sacc.rate_pst_win; ...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro(i) = max(r_exc_pro(i,:)); max_sem_exc_pro(i) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.sacc.rate_pst_win; ...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro(i) = max(abs(r_sup_pro(i,:))); max_sem_sup_pro(i) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.sacc.rate_pst_win; ...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti(i) = max(r_exc_anti(i,:)); max_sem_exc_anti(i) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti(i,:) =  units(indx_sup_anti(i)).anti.neural.sacc.rate_pst_win; ...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti(i) = max(abs(r_sup_anti(i,:))); max_sem_sup_anti(i) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        
        % instr .. psth during instr back window - mean activity during baseline
         t_win = units(1).pro.neural.instr_back.ts_pst_win;
        for i = 1:length(indx_exc_pro), r_exc_pro_instr(i,:) = units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win; ...
                sem_exc_pro_instr(i,:) = std(r_exc_pro_instr(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro_instr(i) = max(r_exc_pro_instr(i,:)); max_sem_exc_pro_instr(i) = std(r_exc_pro_instr(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro_instr(i,:) = units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win; ...
                sem_sup_pro_instr(i,:) = std(r_sup_pro_instr(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro_instr(i) = max(abs(r_sup_pro_instr(i,:))); max_sem_sup_pro_instr(i) = std(r_sup_pro_instr(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti_instr(i,:) = units(indx_exc_anti(i)).anti.neural.instr_back.rate_pst_win; ...
                sem_exc_anti_instr(i,:) = std(r_exc_anti_instr(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti_instr(i) = max(r_exc_anti_instr(i,:)); max_sem_exc_anti_instr(i) = std(r_exc_anti_instr(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti_instr(i,:) =  units(indx_sup_anti(i)).anti.neural.instr_back.rate_pst_win; ...
                sem_sup_anti_instr(i,:) = std(r_sup_anti_instr(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti_instr(i) = max(abs(r_sup_anti_instr(i,:))); max_sem_sup_anti_instr(i) = std(r_sup_anti_instr(i,:))/sqrt(length(indx_sup_anti)); end
        
        if strcmp(recArea,'vermis')
            %plot vermis
            figure; hold on;
            plot(max_exc_pro, max_exc_pro_instr, '.m', 'MarkerSize', 30);
            plot(max_sup_pro, max_sup_pro_instr, '.m', 'MarkerSize', 30);
            plot(max_exc_anti, max_exc_anti_instr, '.b', 'MarkerSize', 30);
            plot(max_sup_anti, max_sup_anti_instr, '.b', 'MarkerSize', 30);
            plot([0 90], [0 90], 'k')
            set(gca,'xlim', [0 90],'xTick',[0 90], 'ylim', [0 90], 'ytick', [0 90],'TickDir', 'out', 'FontSize', 18);
            xlabel('Max change in firing -- Saccade'); ylabel('Max change in firing -- Instruction'); axis square;
            title('vermis')
            
            [~,pro_pval] = ttest2([max_exc_pro' ; max_sup_pro'],[max_exc_pro_instr' ; max_sup_pro_instr'])
            [~,anti_pval] = ttest2([max_exc_anti' ; max_sup_anti'],[max_exc_anti_instr' ; max_sup_anti_instr'])
            
        else
            %plot lateral
            figure; hold on;
            plot(max_exc_pro, max_exc_pro_instr, '.m', 'MarkerSize', 30);
            plot(max_sup_pro, max_sup_pro_instr, '.m', 'MarkerSize', 30);
            plot(max_exc_anti, max_exc_anti_instr, '.b', 'MarkerSize', 30);
            plot(max_sup_anti, max_sup_anti_instr, '.b', 'MarkerSize', 30);
            plot([0 83], [0 83], 'k')
            set(gca,'xlim', [0 83],'xTick',[0 80], 'ylim', [0 83], 'ytick', [0 80],'TickDir', 'out', 'FontSize', 18);
            xlabel('Max change in firing -- Saccade'); ylabel('Max change in firing -- Instruction'); axis square;
            title('lateral')
            
            [~,pro_pval] = ttest2([max_exc_pro' ; max_sup_pro'],[max_exc_pro_instr' ; max_sup_pro_instr'])
            [~,anti_pval] = ttest2([max_exc_anti' ; max_sup_anti'],[max_exc_anti_instr' ; max_sup_anti_instr'])
        end
        
    case 'max_change_instr_sacc_selected'
%         indx_exc_pro = pop.indx_sel.(recArea).instr_back.all.pro.exc;
%         indx_sup_pro = pop.indx_sel.(recArea).instr_back.all.pro.sup;
%         indx_exc_anti = pop.indx_sel.(recArea).instr_back.all.anti.exc;
%         indx_sup_anti = pop.indx_sel.(recArea).instr_back.all.anti.sup;
        
        indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
        
        % try it with fac y supp in sacc sel cells
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=0 & t<=0.150);
        % sacc .... psth during sacc window - mean of activity before sacc onset 
        for i = 1:length(indx_exc_pro), r_exc_pro(i,:) = units(indx_exc_pro(i)).pro.neural.sacc.rate_pst_win - mean(units(indx_exc_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0));...
                sem_exc_pro(i,:) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro(i) = max(r_exc_pro(i,:)); max_sem_exc_pro(i) = std(r_exc_pro(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro(i,:) = units(indx_sup_pro(i)).pro.neural.sacc.rate_pst_win - mean(units(indx_sup_pro(i)).pro.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_sup_pro(i,:) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro(i) = max(abs(r_sup_pro(i,:))); max_sem_sup_pro(i) = std(r_sup_pro(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti(i,:) = units(indx_exc_anti(i)).anti.neural.sacc.rate_pst_win - mean(units(indx_exc_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_exc_anti(i,:) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti(i) = max(r_exc_anti(i,:)); max_sem_exc_anti(i) = std(r_exc_anti(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti(i,:) =  units(indx_sup_anti(i)).anti.neural.sacc.rate_pst_win - mean(units(indx_sup_anti(i)).anti.neural.sacc.rate_pst(1,t>=-0.151 & t<=0)); ...
                sem_sup_anti(i,:) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti(i) = max(abs(r_sup_anti(i,:))); max_sem_sup_anti(i) = std(r_sup_anti(i,:))/sqrt(length(indx_sup_anti)); end
        
        % instr .. psth during instr back window - mean activity during baseline
         t_win = units(1).pro.neural.instr_back.ts_pst_win;
        for i = 1:length(indx_exc_pro), r_exc_pro_instr(i,:) = units(indx_exc_pro(i)).pro.neural.instr_back.rate_pst_win - units(indx_exc_pro(i)).pro.neural.base.rate_mu; ...
                sem_exc_pro_instr(i,:) = std(r_exc_pro_instr(i,:))/sqrt(length(indx_exc_pro)); max_exc_pro_instr(i) = max(r_exc_pro_instr(i,:)); max_sem_exc_pro_instr(i) = std(r_exc_pro_instr(i,:))/sqrt(length(indx_exc_pro)); end
        for i = 1:length(indx_sup_pro), r_sup_pro_instr(i,:) = units(indx_sup_pro(i)).pro.neural.instr_back.rate_pst_win - units(indx_sup_pro(i)).pro.neural.base.rate_mu; ...
                sem_sup_pro_instr(i,:) = std(r_sup_pro_instr(i,:))/sqrt(length(indx_sup_pro)); max_sup_pro_instr(i) = max(abs(r_sup_pro_instr(i,:))); max_sem_sup_pro_instr(i) = std(r_sup_pro_instr(i,:))/sqrt(length(indx_sup_pro)); end
        for i = 1:length(indx_exc_anti), r_exc_anti_instr(i,:) = units(indx_exc_anti(i)).anti.neural.instr_back.rate_pst_win - units(indx_exc_anti(i)).anti.neural.base.rate_mu; ...
                sem_exc_anti_instr(i,:) = std(r_exc_anti_instr(i,:))/sqrt(length(indx_exc_anti)); max_exc_anti_instr(i) = max(r_exc_anti_instr(i,:)); max_sem_exc_anti_instr(i) = std(r_exc_anti_instr(i,:))/sqrt(length(indx_exc_anti)); end
        for i = 1:length(indx_sup_anti), r_sup_anti_instr(i,:) =  units(indx_sup_anti(i)).anti.neural.instr_back.rate_pst_win - units(indx_sup_anti(i)).anti.neural.base.rate_mu; ...
                sem_sup_anti_instr(i,:) = std(r_sup_anti_instr(i,:))/sqrt(length(indx_sup_anti)); max_sup_anti_instr(i) = max(abs(r_sup_anti_instr(i,:))); max_sem_sup_anti_instr(i) = std(r_sup_anti_instr(i,:))/sqrt(length(indx_sup_anti)); end
        
        if strcmp(recArea,'vermis')
            %plot vermis
            figure; hold on;
            plot(max_exc_pro, max_exc_pro_instr, '.m', 'MarkerSize', 30);
            plot(max_sup_pro, max_sup_pro_instr, 'om', 'MarkerSize', 10);
            plot(max_exc_anti, max_exc_anti_instr, '.b', 'MarkerSize', 30);
            plot(max_sup_anti, max_sup_anti_instr, 'ob', 'MarkerSize', 10);
            plot([0 90], [0 90], 'k')
            set(gca,'xlim', [0 90],'xTick',[0 90], 'ylim', [0 90], 'ytick', [0 90],'TickDir', 'out', 'FontSize', 18);
            xlabel('Max change in firing -- Saccade'); ylabel('Max change in firing -- Instruction'); axis square;
            title('vermis')
            
            [~,pro_pval] = ttest2([max_exc_pro' ; max_sup_pro'],[max_exc_pro_instr' ; max_sup_pro_instr'])
            [~,anti_pval] = ttest2([max_exc_anti' ; max_sup_anti'],[max_exc_anti_instr' ; max_sup_anti_instr'])
            fac_pro
            
            
            % see if change in firing is more in pro vs anti from fac or supp
            fac_pro = max_exc_pro'-max_exc_pro_instr';
            fac_anti = max_exc_anti'-max_exc_anti_instr'; fac_anti(6)=NaN;
            sup_pro = max_sup_pro'-max_sup_pro_instr';
            sup_anti = max_sup_anti'-max_sup_anti_instr';
            
            figure; boxplot(fac_pro); ylim([-30 46]); 
            set(gca,'TickDir', 'out', 'FontSize', 18);
            ylabel('Saccade-Instr'); box off; axis square;
            title('facilitation pro'); hline(0, '--k')
            
            figure; boxplot(fac_anti); ylim([-30 46]); 
            set(gca,'TickDir', 'out', 'FontSize', 18);
            ylabel('Saccade-Instr'); box off; axis square;
            title('facilitation anti'); hline(0, '--k')
                        
            figure; boxplot(sup_pro); ylim([-30 46]); 
            set(gca,'TickDir', 'out', 'FontSize', 18);
            ylabel('Saccade-Instr'); box off; axis square;
            title('suppression pro'); hline(0, '--k')
            
            figure; boxplot(sup_anti); ylim([-30 46]); 
            set(gca,'TickDir', 'out', 'FontSize', 18);
            ylabel('Saccade-Instr'); box off; axis square;
            title('suppression anti'); hline(0, '--k')

            ranksum(fac_pro,fac_anti)
            ranksum(sup_pro,sup_anti)
            ranksum(fac_pro, sup_pro)
            ranksum(fac_anti, sup_anti)

        else
            %plot lateral
            figure; hold on;
            plot(max_exc_pro, max_exc_pro_instr, '.m', 'MarkerSize', 30);
            plot(max_sup_pro, max_sup_pro_instr, 'om', 'MarkerSize', 10);
            plot(max_exc_anti, max_exc_anti_instr, '.b', 'MarkerSize', 30);
            plot(max_sup_anti, max_sup_anti_instr, 'ob', 'MarkerSize', 10);
            plot([0 83], [0 83], 'k')
            set(gca,'xlim', [0 83],'xTick',[0 80], 'ylim', [0 83], 'ytick', [0 80],'TickDir', 'out', 'FontSize', 18);
            xlabel('Max change in firing -- Saccade'); ylabel('Max change in firing -- Instruction'); axis square;
            title('lateral')
            
            [~,pro_pval] = ttest2([max_exc_pro' ; max_sup_pro'],[max_exc_pro_instr' ; max_sup_pro_instr'])
            [~,anti_pval] = ttest2([max_exc_anti' ; max_sup_anti'],[max_exc_anti_instr' ; max_sup_anti_instr'])
            
            fac_pro = max_exc_pro'-max_exc_pro_instr';
            fac_anti = max_exc_anti'-max_exc_anti_instr';
            sup_pro = max_sup_pro'-max_sup_pro_instr';
            sup_anti = max_sup_anti'-max_sup_anti_instr';
            
            figure; boxplot(fac_pro); ylim([-30 46]); 
            set(gca,'TickDir', 'out', 'FontSize', 18);
            ylabel('Saccade-Instr'); box off; axis square;
            title('facilitation pro'); hline(0, '--k')
            
            figure; boxplot(fac_anti); ylim([-30 46]); 
            set(gca,'TickDir', 'out', 'FontSize', 18);
            ylabel('Saccade-Instr'); box off; axis square;
            title('facilitation anti'); hline(0, '--k')
                        
            figure; boxplot(sup_pro); ylim([-30 46]); 
            set(gca,'TickDir', 'out', 'FontSize', 18);
            ylabel('Saccade-Instr'); box off; axis square;
            title('suppression pro'); hline(0, '--k')
            
            figure; boxplot(sup_anti); ylim([-30 46]);
            set(gca,'TickDir', 'out', 'FontSize', 18);
            ylabel('Saccade-Instr'); box off; axis square;
            title('suppression anti'); hline(0, '--k')
            
            ranksum(fac_pro,fac_anti)
            ranksum(sup_pro,sup_anti)
            ranksum(fac_pro, sup_pro)
            ranksum(fac_anti, sup_anti)
        end
        
        
        
    case 'change_anti-pro'
        % gather indx
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        % gather indx separately
        
        for cellNum = 1:length(units)
            indx_vermis(cellNum) = strcmp(units(cellNum).area, 'vermis');
            indx_lat(cellNum) = strcmp(units(cellNum).area, 'lateral');
        end
        indx_vermis = find(indx_vermis); indx_lat = find(indx_lat);
        
        
        t = units(1).pro.neural.sacc.ts_pst_win;
        for i = 1:length(indx_area)
            r_pro(i,:) = units(indx_area(i)).pro.neural.sacc.delta_rate_base;
            r_anti(i,:) = units(indx_area(i)).anti.neural.sacc.delta_rate_base;
            indx_sign(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
            indx_sacc_ks_nspk(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk);
            %indx_sacc_ks_tspk(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc_ks_t_spk);
        end
        
        % per area
        for i = 1:length(indx_vermis)
            r_pro_vermis(i,:) = units(indx_vermis(i)).pro.neural.sacc.delta_rate_base;
            r_anti_vermis(i,:) = units(indx_vermis(i)).anti.neural.sacc.delta_rate_base;
            indx_sign_vermis(i) = logical(units(indx_vermis(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk);
        end
        for i = 1:length(indx_lat)
            r_pro_lat(i,:) = units(indx_lat(i)).pro.neural.sacc.delta_rate_base;
            r_anti_lat(i,:) = units(indx_lat(i)).anti.neural.sacc.delta_rate_base;
            indx_sign_lat(i) = logical(units(indx_lat(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk);
        end
        
        % plot all cells
        %         figure; hold on;
        %         plot(t,mean(abs(r_anti))-mean(abs(r_pro)),'Color','m', 'LineWidth', 2);
        %         set(gca,'xlim', [-0.150 0.151],'TickDir','out','ylim',[-1 3], 'ytick',[-1 0 3], 'FontSize', 26)
        %         xlabel('Time (s)'); ylabel('Abs change in FR anti-pro'); title('all cells')
        %
        %         % plot all signif cells
        %         figure; hold on;
        %         plot(t,mean(abs(r_anti(indx_sign,:)))-mean(abs(r_pro(indx_sign,:))),'Color','k', 'LineWidth', 2);
        %         set(gca,'xlim', [-0.150 0.151],'TickDir','out','ylim',[-1 4], 'ytick',[-1 0 4], 'FontSize', 26)
        %         xlabel('Time (s)'); ylabel('Abs change in FR anti-pro')
        
        % plot together omv and lat
        figure; hold on;
        plot(t,mean(abs(r_anti_vermis(indx_sign_vermis,:)))-mean(abs(r_pro_vermis(indx_sign_vermis,:))),'-k', 'LineWidth', 2);
        plot(t,mean(abs(r_anti_lat(indx_sign_lat,:)))-mean(abs(r_pro_lat(indx_sign_lat,:))),'--k', 'LineWidth', 2);
        set(gca,'xlim', [-0.150 0.151],'TickDir','out','ylim',[0 7], 'ytick',[0 7], 'FontSize', 26)
        xlabel('Time (s)'); ylabel('Abs change in FR anti-pro')
        
    case 'change_anti-pro_selected'
        
        
        
    case 'change_anti-pro_instr'
        % gather indx
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        % gather indx separately
        for cellNum = 1:length(units)
            indx_vermis(cellNum) = strcmp(units(cellNum).area, 'vermis');
            indx_lat(cellNum) = strcmp(units(cellNum).area, 'lateral');
        end
        indx_vermis = find(indx_vermis); indx_lat = find(indx_lat);
        
        
        t = units(1).pro.neural.instr.ts_pst_win;
        for i = 1:length(indx_area)
            r_pro(i,:) = units(indx_area(i)).pro.neural.instr.delta_rate_base;
            r_anti(i,:) = units(indx_area(i)).anti.neural.instr.delta_rate_base;
            indx_sign(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr);
            indx_instr_ks_nspk(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr_ks_nspk);
        end
        
        % per area
        for i = 1:length(indx_vermis)
            r_pro_vermis(i,:) = units(indx_vermis(i)).pro.neural.instr.delta_rate_base;
            r_anti_vermis(i,:) = units(indx_vermis(i)).anti.neural.instr.delta_rate_base;
            indx_sign_vermis(i) = logical(units(indx_vermis(i)).stats.instr.flags.proVsAnti_instr_ks_nspk);
        end
        for i = 1:length(indx_lat)
            r_pro_lat(i,:) = units(indx_lat(i)).pro.neural.instr.delta_rate_base;
            r_anti_lat(i,:) = units(indx_lat(i)).anti.neural.instr.delta_rate_base;
            indx_sign_lat(i) = logical(units(indx_lat(i)).stats.instr.flags.proVsAnti_instr_ks_nspk);
        end
        
        % plot all cells
        figure; hold on;
        plot(t,mean(abs(r_anti))-mean(abs(r_pro)),'Color','m', 'LineWidth', 2);
        set(gca,'xlim', [0.05 0.350],'TickDir','out','ylim',[-1 3], 'ytick',[-1 0 3], 'FontSize', 26)
        xlabel('Time (s)'); ylabel('Abs change in FR anti-pro'); title('all cells')
        
        % plot all signif cells
        figure; hold on;
        plot(t,mean(abs(r_anti(indx_sign,:)))-mean(abs(r_pro(indx_sign,:))),'Color','k', 'LineWidth', 2);
        set(gca,'xlim', [0.05 0.350],'TickDir','out','ylim',[-1 4], 'ytick',[-1 0 4], 'FontSize', 26)
        xlabel('Time (s)'); ylabel('Abs change in FR anti-pro')
        
        % plot together omv and lat
        figure; hold on;
        plot(t,mean(abs(r_anti_vermis(indx_sign_vermis,:)))-mean(abs(r_pro_vermis(indx_sign_vermis,:))),'-k', 'LineWidth', 2);
        plot(t,mean(abs(r_anti_lat(indx_sign_lat,:)))-mean(abs(r_pro_lat(indx_sign_lat,:))),'--k', 'LineWidth', 2);
        set(gca,'xlim', [0.05 0.350],'TickDir','out','ylim',[-3 3], 'ytick',[-3 0 3], 'FontSize', 26);
        xlabel('Time (s)'); ylabel('Abs change in FR anti-pro')
        
    case 'max_delta_rate'
        fprintf(['        >>> loading ' recArea ' cells <<< \n']);
        for cellNum = 1:length(units)
            indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx_area = find(indx_area);
        
        t = units(1).pro.neural.sacc.ts_pst_win;
        for i = 1:length(indx_area)
            max_pro(i) = max(abs(units(indx_area(i)).pro.neural.sacc.delta_rate_base));
            max_anti(i) = max(abs(units(indx_area(i)).anti.neural.sacc.delta_rate_base));
        end
        
        for i = 1:length(indx_area)
            indx_sign(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc);
            indx_sacc_ks_nspk(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk);
            indx_sacc_ks_tspk(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc_ks_t_spk);
        end
        
        figure; hold on;
        plot(max_pro,max_anti, '.k','MarkerSize', 18);
        %plot(max_pro(indx_sacc_ks_tspk),max_anti(indx_sacc_ks_tspk), '.c','MarkerSize', 18);
        set(gca,'XScale','Log','YScale','Log' ,'FontSize', 18, 'TickDir', 'out');axis ([1e0 1e2 1e0 1e2]);
        plot([1e0 1e2],[1e0 1e2]);
        xlabel('Max change pro'); ylabel('Max change anti');
        % title(['Max change in firing rate from base >> ' recArea])
        axis square; %title('ts_ks')
        [h,p] = ttest(max_pro,max_anti)
        
        
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
        
    case 'latencies_sacc_selected'
         % gather relevant neurons
        indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc;
        indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup;
        indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
        indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
        
         t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=-0.151 & t<=0.151);
          for i = 1:length(indx_exc_pro), r_exc_pro(i) = units(cellNum).pro.neural.sacc.peak_resp_time; end
          for i = 1:length(indx_sup_pro), r_sup_pro(i) = units(cellNum).pro.neural.sacc.peak_resp_time; end
          for i = 1:length(indx_exc_anti), r_exc_anti(i) = units(cellNum).anti.neural.sacc.peak_resp_time; end
          for i = 1:length(indx_sup_anti), r_sup_anti(i) = units(cellNum).anti.neural.sacc.peak_resp_time; end
          
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
            r(j,:) = units(indx_area(j)).pro.neural.sacc.norm.rate_pst(1,:);
        end
        
        [r,t] = smooth_colormap(r,t);
        
        % plot colormap
        B = goodcolormap('mk');  %B = goodcolormap('bwr');
        figure; set(gcf,'Position',[100 200 300 300]);
        hold on; colormap(B');
        % hold on; colormap(bluewhitered);
        % hold on; %colormap(bone(4));
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
        [rho_pro_cell,pval_pro_cell] = corr(amp_pro',spk_pro);
        [rho_anti_cell,pval_anti_cell] = corr(amp_anti',spk_anti);
        
        %plot
        figure; hold on; box off
        plot(amp_pro,spk_pro, '.r', 'MarkerSize',16);
        plot(amp_anti,spk_anti, '.g','MarkerSize',16);
        xlabel('Amplitude (deg)'); ylabel('Firing rate (spk/s)');
        set(gca, 'TickDir', 'out', 'FontSize', 18);
        lsline(gca);
        
        % gather and plot for all cells
        amp_pro = []; amp_anti=[]; spk_pro=[]; spk_anti=[];
        for cell = 1:length(units)
            indx_area(cell) = strcmp(units(cell).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        for i=1:length(nunits)
            temp_amps_pro = [units(indx_area(i)).pro.behav.trial.saccAmplitude]';
            temp_amps_anti = [units(indx_area(i)).anti.behav.trial.saccAmplitude]';
            temp_spks_pro = [units(indx_area(i)).pro.neural.sacc.nspk];
            temp_spks_anti = [units(indx_area(i)).anti.neural.sacc.nspk];
            
            amp_pro = [amp_pro ; temp_amps_pro];
            amp_anti = [amp_anti ; temp_amps_anti];
            spk_pro = [spk_pro ; temp_spks_pro];
            spk_anti = [spk_anti ; temp_spks_anti];
        end
        %plot
        figure; hold on;
        plot(spk_anti,amp_anti, '.b', 'MarkerSize',16);
        plot(spk_pro,amp_pro, '.m', 'MarkerSize',16);
        set(gca, 'xlim',[0 250], 'TickDir', 'out', 'FontSize', 18);
        xlabel('Firing rate (spks/s)'); ylabel('Amplitude (deg)');
        [rho_pro,pval_pro] = corr(amp_pro,spk_pro)
        [rho_anti,pval_anti] = corr(amp_anti,spk_anti)
        
        %Avg per amplitude blocks
        block1 = 4; %4.9;
        block2 = 6; %9.9;
        block3 = 8; %15;
        
        %pro
        block1_pro_spk = spk_pro(amp_pro<block1);
        block2_pro_spk = spk_pro(amp_pro>=block1 & amp_pro<block2);
        block3_pro_spk = spk_pro(amp_pro>=block2 & amp_pro<block3);
        block4_pro_spk = spk_pro(amp_pro>=block3);
        
        block1_pro_amp = amp_pro(amp_pro<block1);
        block2_pro_amp = amp_pro(amp_pro>=block1 & amp_pro<block2);
        block3_pro_amp = amp_pro(amp_pro>=block2 & amp_pro<block3);
        block4_pro_amp = amp_pro(amp_pro>=block3);
        
        % anti
        block1_anti_spk = spk_pro(amp_anti<block1);
        block2_anti_spk = spk_pro(amp_anti>=block1 & amp_anti<block2);
        block3_anti_spk = spk_pro(amp_anti>=block2 & amp_anti<block3);
        block4_anti_spk = spk_pro(amp_anti>=block3);
        
        block1_anti_amp = amp_anti(amp_anti<block1);
        block2_anti_amp = amp_anti(amp_anti>=block1 & amp_anti<block2);
        block3_anti_amp = amp_anti(amp_anti>=block2 & amp_anti<block3);
        block4_anti_amp = amp_anti(amp_anti>=block3);
        
        % plot pro
        %figure; hold on;
        errorbar(mean(block1_pro_spk),mean(block1_pro_amp), std(block1_pro_amp), std(block1_pro_amp),std(block1_pro_spk),std(block1_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block2_pro_spk),mean(block2_pro_amp), std(block2_pro_amp), std(block2_pro_amp),std(block2_pro_spk),std(block2_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block3_pro_spk),mean(block3_pro_amp), std(block3_pro_amp), std(block3_pro_amp),std(block3_pro_spk),std(block3_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block4_pro_spk),mean(block4_pro_amp), std(block4_pro_amp), std(block4_pro_amp),std(block4_pro_spk),std(block4_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        
        % plot anti
        errorbar(mean(block1_anti_spk),mean(block1_anti_amp), std(block1_anti_amp), std(block1_anti_amp),std(block1_anti_spk),std(block1_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block2_anti_spk),mean(block2_anti_amp), std(block2_anti_amp), std(block2_anti_amp),std(block2_anti_spk),std(block2_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block3_anti_spk),mean(block3_anti_amp), std(block3_anti_amp), std(block3_anti_amp),std(block3_anti_spk),std(block3_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block4_anti_spk),mean(block4_anti_amp), std(block4_anti_amp), std(block4_anti_amp),std(block4_anti_spk),std(block4_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        
        set(gca, 'xlim',[0 250], 'ylim', [0 30]);
        
        
        
        %         %pro
        %         block1_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]<=block1);
        %         block2_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]>block1 & [units(cellNum).pro.behav.trial.saccAmplitude]<= block2);
        %         block3_pro_indx = find([units(cellNum).pro.behav.trial.saccAmplitude]>block2 & [units(cellNum).pro.behav.trial.saccAmplitude]<= block3);
        %
        %         for i = 1:length(block1_pro_indx)
        %             block1_pro(i,1) = units(cellNum).pro.neural.sacc.nspk(block1_pro_indx(i));
        %         end
        %         for i = 1:length(block2_pro_indx)
        %             block2_pro(i,1) = units(cellNum).pro.neural.sacc.nspk(block2_pro_indx(i));
        %         end
        %         for i = 1:length(block3_pro_indx)
        %             block3_pro(i,1) = units(cellNum).pro.neural.sacc.nspk(block1_pro_indx(i));
        %         end
        %
        %
        %         %anti
        %         block1_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]<=block1);
        %         block2_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]>block1 & [units(cellNum).anti.behav.trial.saccAmplitude]<= block2);
        %         block3_anti_indx = find([units(cellNum).anti.behav.trial.saccAmplitude]>block2 & [units(cellNum).anti.behav.trial.saccAmplitude]<= block3);
        %
        %         if ~isempty(block1_anti_indx)
        %             for i = 1:length(block1_anti_indx)
        %                 block1_anti(i,1) = units(cellNum).anti.neural.sacc.nspk(block1_anti_indx(i));
        %             end
        %         else
        %             block1_anti(i,1) = NaN;
        %         end
        %         for i = 1:length(block2_anti_indx)
        %             block2_anti(i,1) = units(cellNum).anti.neural.sacc.nspk(block2_anti_indx(i));
        %         end
        %         for i = 1:length(block3_anti_indx)
        %             block3_anti(i,1) = units(cellNum).anti.neural.sacc.nspk(block3_anti_indx(i));
        %         end
        
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
%         
%         % plot rho for all cells
%         fprintf(['        >>> loading ' recArea ' cells <<< \n']);
%         for cellNum = 1:length(units)
%             indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
%         end
%         indx_area = find(indx_area);
%         nunits = 1:length(indx_area);
%         
%         for i = 1:length(units(indx_area))
%             amp_pro_all{i,:} = [units(indx_area(i)).pro.behav.trial.saccAmplitude];
%             amp_anti_all{i,:} = [units(indx_area(i)).anti.behav.trial.saccAmplitude];
%             
%             spk_pro_all{i,:} = [units(indx_area(i)).pro.neural.sacc.nspk];
%             spk_anti_all{i,:} = [units(indx_area(i)).anti.neural.sacc.nspk];
%             
%             [rho_pro_all(i,:),pval_pro_all(i,:)] = corr(amp_pro_all{i,:}',spk_pro_all{i,:});
%             [rho_anti_all(i,:),pval_anti_all(i,:)] = corr(amp_anti_all{i,:}',spk_anti_all{i,:});
%         end
        
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
        
        % gather and plot for all cells
        x_pro = []; x_anti=[]; spk_pro=[]; spk_anti=[];
        for cell = 1:length(units)
            indx_area(cell) = strcmp(units(cell).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        for i=1:length(nunits)
            temp_x_pro = [units(indx_area(i)).pro.behav.trial.saccPeakVel]';
            temp_x_anti = [units(indx_area(i)).anti.behav.trial.saccPeakVel]';
            temp_spks_pro = [units(indx_area(i)).pro.neural.sacc.nspk];
            temp_spks_anti = [units(indx_area(i)).anti.neural.sacc.nspk];
            
            x_pro = [x_pro ; temp_x_pro];
            x_anti = [x_anti ; temp_x_anti];
            spk_pro = [spk_pro ; temp_spks_pro];
            spk_anti = [spk_anti ; temp_spks_anti];
        end
        %plot
        figure; hold on;
        plot(spk_pro,x_pro, '.m', 'MarkerSize',16);
        plot(spk_anti,x_anti, '.b', 'MarkerSize',16);
        set(gca, 'xlim',[0 250], 'TickDir', 'out', 'FontSize', 18);
        xlabel('Firing rate (spks/s)'); ylabel('Peak Vel (deg/s)');
        [rho_pro,pval_pro] = corr(x_pro,spk_pro)
        [rho_anti,pval_anti] = corr(x_anti,spk_anti)
        
        block1 = 200;
        block2 = 300;
        block3 = 400;
        
        %pro
        block1_pro_spk = spk_pro(x_pro<block1);
        block2_pro_spk = spk_pro(x_pro>=block1 & x_pro<block2);
        block3_pro_spk = spk_pro(x_pro>=block2 & x_pro<block3);
        block4_pro_spk = spk_pro(x_pro>=block3);
        
        block1_pro_amp = x_pro(x_pro<block1);
        block2_pro_amp = x_pro(x_pro>=block1 & x_pro<block2);
        block3_pro_amp = x_pro(x_pro>=block2 & x_pro<block3);
        block4_pro_amp = x_pro(x_pro>=block3);
        
        % anti
        block1_anti_spk = spk_pro(x_anti<block1);
        block2_anti_spk = spk_pro(x_anti>=block1 & x_anti<block2);
        block3_anti_spk = spk_pro(x_anti>=block2 & x_anti<block3);
        block4_anti_spk = spk_pro(x_anti>=block3);
        
        block1_anti_amp = x_anti(x_anti<block1);
        block2_anti_amp = x_anti(x_anti>=block1 & x_anti<block2);
        block3_anti_amp = x_anti(x_anti>=block2 & x_anti<block3);
        block4_anti_amp = x_anti(x_anti>=block3);
        
        % plot pro
%         figure; hold on;
        errorbar(mean(block1_pro_spk),mean(block1_pro_amp), std(block1_pro_amp), std(block1_pro_amp),std(block1_pro_spk),std(block1_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block2_pro_spk),mean(block2_pro_amp), std(block2_pro_amp), std(block2_pro_amp),std(block2_pro_spk),std(block2_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block3_pro_spk),mean(block3_pro_amp), std(block3_pro_amp), std(block3_pro_amp),std(block3_pro_spk),std(block3_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block4_pro_spk),mean(block4_pro_amp), std(block4_pro_amp), std(block4_pro_amp),std(block4_pro_spk),std(block4_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        
        % plot anti
        errorbar(mean(block1_anti_spk),mean(block1_anti_amp), std(block1_anti_amp), std(block1_anti_amp),std(block1_anti_spk),std(block1_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block2_anti_spk),mean(block2_anti_amp), std(block2_anti_amp), std(block2_anti_amp),std(block2_anti_spk),std(block2_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block3_anti_spk),mean(block3_anti_amp), std(block3_anti_amp), std(block3_anti_amp),std(block3_anti_spk),std(block3_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block4_anti_spk),mean(block4_anti_amp), std(block4_anti_amp), std(block4_anti_amp),std(block4_anti_spk),std(block4_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        
        set(gca, 'xlim',[0 250], 'ylim', [0 1000]);
        
        
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
        
        % gather and plot for all cells
        x_pro = []; x_anti=[]; spk_pro=[]; spk_anti=[];
        for cell = 1:length(units)
            indx_area(cell) = strcmp(units(cell).area, recArea);
        end
        indx_area = find(indx_area);
        nunits = 1:length(indx_area);
        for i=1:length(nunits)
            temp_x_pro = [units(indx_area(i)).pro.behav.trial.saccDuration]';
            temp_x_anti = [units(indx_area(i)).anti.behav.trial.saccDuration]';
            temp_spks_pro = [units(indx_area(i)).pro.neural.sacc.nspk];
            temp_spks_anti = [units(indx_area(i)).anti.neural.sacc.nspk];
            
            x_pro = [x_pro ; temp_x_pro];
            x_anti = [x_anti ; temp_x_anti];
            spk_pro = [spk_pro ; temp_spks_pro];
            spk_anti = [spk_anti ; temp_spks_anti];
        end
        %plot
        figure; hold on;
        plot(spk_anti,x_anti, '.b', 'MarkerSize',16);
        plot(spk_pro,x_pro, '.m', 'MarkerSize',16);
        set(gca, 'xlim',[0 250], 'TickDir', 'out', 'FontSize', 18);
        xlabel('Firing rate (spks/s)'); ylabel('Duration (ms)');
        [rho_pro,pval_pro] = corr(x_pro,spk_pro)
        [rho_anti,pval_anti] = corr(x_anti,spk_anti)
        
        block1 = 20;
        block2 = 40;
        block3 = 60;
        
        %pro
        block1_pro_spk = spk_pro(x_pro<block1);
        block2_pro_spk = spk_pro(x_pro>=block1 & x_pro<block2);
        block3_pro_spk = spk_pro(x_pro>=block2 & x_pro<block3);
        block4_pro_spk = spk_pro(x_pro>=block3);
        
        block1_pro_amp = x_pro(x_pro<block1);
        block2_pro_amp = x_pro(x_pro>=block1 & x_pro<block2);
        block3_pro_amp = x_pro(x_pro>=block2 & x_pro<block3);
        block4_pro_amp = x_pro(x_pro>=block3);
        
        % anti
        block1_anti_spk = spk_pro(x_anti<block1);
        block2_anti_spk = spk_pro(x_anti>=block1 & x_anti<block2);
        block3_anti_spk = spk_pro(x_anti>=block2 & x_anti<block3);
        block4_anti_spk = spk_pro(x_anti>=block3);
        
        block1_anti_amp = x_anti(x_anti<block1);
        block2_anti_amp = x_anti(x_anti>=block1 & x_anti<block2);
        block3_anti_amp = x_anti(x_anti>=block2 & x_anti<block3);
        block4_anti_amp = x_anti(x_anti>=block3);
        
        % plot pro
%         figure; hold on;
        errorbar(mean(block1_pro_spk),mean(block1_pro_amp), std(block1_pro_amp), std(block1_pro_amp),std(block1_pro_spk),std(block1_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block2_pro_spk),mean(block2_pro_amp), std(block2_pro_amp), std(block2_pro_amp),std(block2_pro_spk),std(block2_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block3_pro_spk),mean(block3_pro_amp), std(block3_pro_amp), std(block3_pro_amp),std(block3_pro_spk),std(block3_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block4_pro_spk),mean(block4_pro_amp), std(block4_pro_amp), std(block4_pro_amp),std(block4_pro_spk),std(block4_pro_spk),'om', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        
        % plot anti
        errorbar(mean(block1_anti_spk),mean(block1_anti_amp), std(block1_anti_amp), std(block1_anti_amp),std(block1_anti_spk),std(block1_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block2_anti_spk),mean(block2_anti_amp), std(block2_anti_amp), std(block2_anti_amp),std(block2_anti_spk),std(block2_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block3_anti_spk),mean(block3_anti_amp), std(block3_anti_amp), std(block3_anti_amp),std(block3_anti_spk),std(block3_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        errorbar(mean(block4_anti_spk),mean(block4_anti_amp), std(block4_anti_amp), std(block4_anti_amp),std(block4_anti_spk),std(block4_anti_spk),'ob', 'MarkerSize', 18, 'Capsize', 0, 'LineWidth',1);
        
        set(gca, 'xlim',[0 250], 'ylim', [0 200]);
        
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
%             %             %
            stat_instr(i,:) = units(indx_area(i)).stats.instr.pval.pbDist_testStat;
            t_instr = units(indx_area(i)).pro.neural.instr.ts_pst;
%             position_stat_sign_instr(i,:) = abs(stat_instr(i,:))>=1.96;
%             %             %
%             %plot instr
%             subplot(2,1,1); hold on
%             plot(t_instr, stat_instr(i,:), 'k','MarkerSize', 15);
%             hline(-1.96, 'k');hline(1.96, 'k');
%             set(gca, 'xlim',[0 0.350],'xTick',[0 0.175 0.350] ,'ylim',[-6 6], 'TickDir', 'out', 'FontSize', 18)
%             title ('Binomial pb dist - Instruction')
%             xlabel('time (s)'); ylabel('z-stat');
%             
            % Find if there are two consecutive timepoints
            % above or below 1.96
            clear a0 ii
            a0 = abs(stat_instr(i,t_instr>0.049 & t_instr<0.351)); % input vector
            a0 = a0>1.96;
            ii = strfind(a0,[1 1]);
            if ~isempty(ii)
                good_cell_instr(i) = 1;
                good_cell_instr_timepoints{i} = ii;
            else
                good_cell_instr(i) = 0;
                good_cell_instr_timepoints{i} = [];
            end
%             
%             %             %             %
%             %             %             %             % sacc
%             %             %             %             %gather
            t_sacc = units(indx_area(i)).pro.neural.sacc.ts_pst;
            stat_sacc(i,:) = units(indx_area(i)).stats.sacc.pval.pbDist_testStat;
            position_stat_sign_sacc(i,:) = abs(stat_sacc(i,:))>=1.96;
%             %             %             %             % plot
%             subplot(2,1,2); hold on
%             plot(t_sacc, stat_sacc(i,:), 'k','MarkerSize', 15);
%             set(gca, 'xlim',[-0.150 0.151],'xTick',[-0.1 0 0.1],'ylim',[-8 8], 'TickDir', 'out', 'FontSize', 18)
%             title (['Binomial pb dist - Saccade cell: ' num2str(indx_area(i))])
%             xlabel('time (s)'); ylabel('z-stat'); hline(-1.96, 'k'); hline(1.96, 'k');vline(0, 'c');
%             fname = 'Binomial_pb_dist';
%             print(fname,'-append', '-dpsc2')
%             %waitforbuttonpress; %close all;
%             
            % Find if there are two consecutive timepoints
            % above or below 1.96
            clear a0 ii
            a0 = abs(stat_sacc(i,t_sacc>-0.151 & t_sacc<0.151)); % input vector
            a0 = a0>1.96;
            ii = strfind(a0,[1 1]);
            if ~isempty(ii)
                good_cell_sacc(i) = 1;
                good_cell_sacc_timepoints{i} = ii;
            else
                good_cell_sacc(i) = 0;
                good_cell_sacc_timepoints{i} = [];
            end
%             %             %
        end
        
        
        %% Take the absolute value of the Z-statistic and average across significantly different neurons
        %                    instr -> get significantly different neurons
        % get significantly diff instr
        nunits = 1:length(indx_area);
        for i = 1:length(nunits)
            indx_sign_instr(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr_ks_nspk);
            indx_sign_instr_exc(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr_ks_nspk) & units(cellNum).pro.neural.exc==1;
            indx_sign_instr_sup(i) = logical(units(indx_area(i)).stats.instr.flags.proVsAnti_instr_ks_nspk) & units(cellNum).pro.neural.sup==1;
        end
        n_instr = sum(indx_sign_instr);
        
        %% Disused %% Take the absolute value of the Z-statistic and average across neurons - average for each time point and plot it as a function of time
        % This will reveal how well an average neuron can discriminate between the two conditions.
        
        sacc_Z_stat = mean(abs(stat_sacc));
        instr_Z_stat = mean(abs(stat_instr));
        
        % plot instr Z stat
        figure; hold on;
        plot(t_instr,abs(stat_instr), '.', 'MarkerSize', 20, 'Color', [0.5 0.5 0.5]); hline(1.96);  % data points 11 to 41
        plot(t_instr, abs(stat_instr(logical(good_cell_instr),:)), '-', 'LineWidth', 1);
        %plot(t_instr,abs(stat_instr(indx_sign_instr,:)), '.', 'MarkerSize', 20); hold on; hline(1.96);  % data points 11 to 41
        %plot(nanmedian(abs(stat_instr)), '-r', 'LineWidth',2)
        set(gca, 'xlim',[0.05 0.35], 'xTick', [0.05 0.350],'yTick',[0 2 4], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'); xlabel('time (s)'); title('Instruction');
        
        % side histogram
        figure; hold on;
        histogram(abs(stat_instr), 25)
        %histogram(abs(stat_instr(indx_sign_instr,:)), 25)
        set(gca,'Xdir','reverse','TickDir', 'out', 'FontSize', 18);
        
        % Histogram of >1.96 only
        figure; hold on
        histogram(t_instr,stat_instr(position_stat_sign_instr),25);
        plot(t_sacc(),abs(stat_instr(position_stat_sign_instr)), '.k')
        set(gca, 'xlim',[0.05 0.35], 'TickDir', 'out', 'FontSize', 18);
        
        figure; hold on
        histogram(abs(stat_instr(position_stat_sign_instr)),25);
        set(gca,'TickDir', 'out', 'FontSize', 18);
        
        % plot sacc Z stat
        figure; hold on;
        plot(t_sacc,abs(stat_sacc), '.', 'MarkerSize', 20, 'Color', [0.5 0.5 0.5]); hline(1.96);  % data points 11 to 41
        plot(t_sacc, abs(stat_sacc(logical(good_cell_sacc),:)), '-', 'LineWidth', 1);
        %plot(nanmedian(abs(stat_instr)), '-r', 'LineWidth',2)
        set(gca, 'xlim',[-0.150 0.150],'xTick',[-0.1 0 0.1],'ytick',[0 2 4], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade');
        
        % side histogram
        figure; hold on;
        histogram(abs(stat_sacc), 25)
        set(gca,'Xdir','reverse','TickDir', 'out', 'FontSize', 18);
        
        %% exc and sup z stat
        cnt_exc=1; cnt_sup=1;
        for cellNum = 1:length(units)
            if strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.exc==1
                indx_exc(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
            elseif strcmp(units(cellNum).area, recArea) && units(cellNum).pro.neural.sup==1
                indx_sup(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
            end
        end
        
        % plot
        z_sign_instr = stat_instr(indx_sign_instr,:);
        
        figure; hold on;
        plot(t_instr,smooth(nanmean(abs(z_sign_instr)),3),'LineWidth', 2,'Color','k'); % smoothed
        hline(1.96,'k')
        set(gca, 'xlim',[0 0.3],'ylim',[0 3], 'TickDir', 'out', 'FontSize', 18);
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
        set(gca, 'xlim',[0 0.3],'ylim',[0 3], 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Instr => ' recArea ' n= ' num2str(n_instr)]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        % plot exc instr
        figure; hold on;
        z_sign_instr_exc = stat_instr(indx_sign_instr(),:);
        plot(t_instr,smooth(nanmean(abs(z_stat_exc_signif_instr)),3),'LineWidth', 2,'Color','m'); % smoothed
        set(gca, 'xlim',[-0.150 0.151], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Sacc exc => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot sup instr
        plot(t_instr,smooth(nanmean(abs(z_stat_sup_signif_instr)),3),'LineWidth', 2,'Color','b'); % smoothed
        set(gca, 'xlim',[-0.150 0.151], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Exc and sup => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        
        
        %   sacc -> get significantly different neurons
        for i = 1:length(nunits)
            indx_sign_sacc(i) = logical(units(indx_area(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk);
        end
        n_sacc = sum(indx_sign_sacc);
        
        % get Z stat for those neurons and plot
        z_sign_sacc = stat_sacc(indx_sign_sacc,:);
        figure; hold on;
        plot(t_sacc,smooth(nanmean(abs(z_sign_sacc)),3),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[-0.150 0.151], 'TickDir', 'out', 'FontSize', 18);
        hline(1.96,'--k'); ylim([0.6 2.4])
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
        set(gca, 'xlim',[-0.150 0.151], 'ylim', ([0.6 5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Abs Z stat Sacc => ' recArea ' n= ' num2str(n_sacc)]);
        xlabel('time(s)'); ylabel('Z-stat')
        
        %% Same as above but take mean of Z-stat only on sacc window (0-0.2s)
        sacc_win = t_sacc>-0.100 & t_sacc<0.21;
        t_sacc_win = t_sacc(sacc_win);
        z_sign_sacc_win = z_sign_sacc(:,sacc_win);
        
        for i=1:size(z_sign_sacc_win,1)
            z_win_mu(i) = mean(abs(z_sign_sacc_win(i,:))); % mean for every cell in its window
        end
        
        %         figure; hold on;
        %         plot(z_win_mu,'.k', 'MarkerSize', 18)
        %         set(gca, 'TickDir', 'out', 'FontSize', 18); hline(1.96)
        %         xlabel('cell'); ylabel('Z-stat')
        %         title(['Z stat sacc window(0.1-0.2s) recarea => ' recArea])
        
        
        % plot exc
        figure; hold on;
        plot(t_sacc,smooth(nanmean(abs(z_stat_exc_signif)),3),'LineWidth', 2,'Color','m'); % smoothed
        set(gca, 'xlim',[-0.150 0.151], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Sacc exc => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot sup
        plot(t_sacc,smooth(nanmean(abs(z_stat_sup_signif)),3),'LineWidth', 2,'Color','b'); % smoothed
        set(gca, 'xlim',[-0.150 0.151], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Exc and sup => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
    case 'binomial_pb_dist_selected'
%         indx_exc = pop.indx_sel.(recArea).sacc.all.exc;
%         indx_sup = pop.indx_sel.(recArea).sacc.all.sup;
        indx_exc = pop.indx_sel.(recArea).sacc.all.exc_signif;
        indx_sup = pop.indx_sel.(recArea).sacc.all.sup_signif;
        win = [-0.151 0.151]; 
        
        t_sacc = units(1).pro.neural.sacc.ts_pst;
        
        for i=1:length(indx_exc), stat_exc(i,:) = units(indx_exc(i)).stats.sacc.pval.pbDist_testStat; end
        for i=1:length(indx_sup), stat_sup(i,:) = units(indx_sup(i)).stats.sacc.pval.pbDist_testStat; end
        
         % plot exc
        figure; hold on;
        plot(t_sacc,(nanmean(abs(stat_exc))));
        plot(t_sacc,smooth(nanmean(abs(stat_exc)),'sgolay'),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Sacc exc => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot sup
        %plot(t_sacc,(nanmean(abs(stat_sup))))
        plot(t_sacc,smooth(nanmean(abs(stat_sup)),'sgolay'),'--k','LineWidth', 2); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Exc and sup => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot all
        stat_all = ([abs(stat_exc) ; abs(stat_sup)]); 
        figure; hold on;
        %plot(t_sacc,(nanmean(stat_all)));
        plot(t_sacc,smooth(nanmean(stat_all),3),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Sacc all => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        %% plot data points per cell
        % find cells with two consecutive timepoints above or below 1.96
        % exc
        a0 = abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))); % input vector
        a0 = a0>1.96;
        for j = 1:length(a0(:,1))
            ii = strfind(a0(j,:),[1 1]);
            if ~isempty(ii)
                good_cell_sacc(j) = 1;
                good_cell_sacc_timepoints{j} = ii;
            else
                good_cell_sacc(j) = 0;
                good_cell_sacc_timepoints{j} = [];
            end
        end
        % sup
        a0_sup = abs(stat_sup(:,t_sacc>win(1) & t_sacc<win(2))); % input vector
        a0_sup = a0_sup>1.96;
        for j = 1:length(a0_sup(:,1))
            ii = strfind(a0_sup(j,:),[1 1]);
            if ~isempty(ii)
                good_cell_sacc_sup(j) = 1;
                good_cell_sacc_timepoints_sup{j} = ii;
            else
                good_cell_sacc_sup(j) = 0;
                good_cell_sacc_timepoints_sup{j} = [];
            end
        end
        
        
        stat_exc(isnan(stat_exc))=0; stat_sup(isnan(stat_sup))=0;
        % plot exc
        figure; hold on;
        plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))), '.', 'MarkerSize', 20, 'Color', [0.5 0.5 0.5]); hline(1.96, '--k');
       % plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))), 'MarkerSize', 20); hline(1.96);
        plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2))), '-', 'LineWidth', 2);
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 5],'ytick',[0 5], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; 
        % sum and median per bin
        figure; plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),nanmedian(abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 3],'ytick',[0 3], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; title('median exc'); hline(1.96,'--k')
        
        figure; plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),nansum(abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 30],'ytick',[0 30], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; title('sum exc');
        
        % all selected cells
        pcolor(t_sacc(t_sacc>win(1) & t_sacc<win(2)),1:size(abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2)))),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))));
        shading interp;
        set(gca, 'clim', [0 3]); % 3.92 max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar;
        
        % cells with two consecutive timepoints above or below 1.96
        pcolor(t_sacc(t_sacc>win(1) & t_sacc<win(2)),1:size(abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2)))),abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2))));
        shading interp;
        set(gca, 'clim', [0 3.92], 'yTick', [1 size(abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2))),1)]); % max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar; title('exc')
        
        % plot sup
        figure; hold on;
        plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_sup(:,t_sacc>win(1) & t_sacc<win(2))), '.', 'MarkerSize', 20, 'Color', [0.5 0.5 0.5]); hline(1.96, '--k');
       % plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))), 'MarkerSize', 20); hline(1.96);
        plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2))), '-', 'LineWidth', 2);
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 5],'ytick',[0 5], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; hline(1.96,'--k')
         % sum and median per bin
        figure; plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),median(abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 3],'ytick',[0 3], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; title('median sup');hline(1.96,'--k')
        
        figure; plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),sum(abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 30],'ytick',[0 30], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; title('sum sup');
        
         % cells with two consecutive timepoints above or below 1.96
        pcolor(t_sacc(t_sacc>win(1) & t_sacc<win(2)),1:size(abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2)))),abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2))));
        shading interp;
        set(gca, 'clim', [0 3.92], 'yTick', [1 size(abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2))),1)]); % max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar; title('sup')
        
        
        
        
    case 'binomial_pb_dist_selected_instr' 
        %         indx_exc = pop.indx_sel.(recArea).sacc.all.exc;
        %         indx_sup = pop.indx_sel.(recArea).sacc.all.sup;
        
        %         indx_exc = pop.indx_sel.(recArea).instr_back.all.exc_signif; % also significantly diff between pro and anti
        %         indx_sup = pop.indx_sel.(recArea).instr_back.all.sup_signif; % also significantly diff between pro and anti
        
        indx_exc = pop.indx_sel.(recArea).instr_back.all.exc;
        indx_sup = pop.indx_sel.(recArea).instr_back.all.sup;
        
        
        t_instr = units(1).pro.neural.instrDir.ts_spkCount(1,:);
        t_sacc = units(1).pro.neural.sacc.ts_pst;
        
        for i=1:length(indx_exc), stat_exc(i,:) = units(indx_exc(i)).stats.instr_back.pval.pbDist_testStat; end
        for i=1:length(indx_sup), stat_sup(i,:) = units(indx_sup(i)).stats.instr_back.pval.pbDist_testStat; end
        
        for i=1:length(indx_exc), stat_exc_sacc(i,:) = units(indx_exc(i)).stats.sacc.pval.pbDist_testStat; end
        for i=1:length(indx_sup), stat_sup_sacc(i,:) = units(indx_sup(i)).stats.sacc.pval.pbDist_testStat; end
        
        % plot exc
        figure; hold on;
        plot(t_instr,(nanmean(abs(stat_exc))));
        plot(t_instr,smooth(nanmean(abs(stat_exc)),'sgolay'),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[-0.3 0], 'ylim', ([0 2.5]), 'TickDir', 'out', 'FontSize', 18)
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        title(['instr selected exc ' recArea])
        
        % plot sup
        plot(t_instr,(nanmean(abs(stat_sup))))
        plot(t_instr,smooth(nanmean(abs(stat_sup)),'sgolay'),'--k','LineWidth', 2); % smoothed
        set(gca, 'xlim',[-0.3 0], 'ylim', ([0 2.5]), 'TickDir', 'out', 'FontSize', 18);
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot all
        stat_all = ([abs(stat_exc) ; abs(stat_sup)]); 
        figure; hold on;
        %plot(t_sacc,(nanmean(stat_all)));
        plot(t_instr,smooth(nanmean(stat_all),3),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[-0.150 0.151], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Sacc all => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
           
        %% sacc
        % plot exc
        figure; hold on;
        plot(t_sacc,(nanmean(abs(stat_exc_sacc))));
        plot(t_sacc,smooth(nanmean(abs(stat_exc_sacc)),'sgolay'),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[-0.2 0.2], 'ylim', ([0 6]), 'TickDir', 'out', 'FontSize', 18);
        title(['Sacc exc instr selected => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot sup
        plot(t_sacc,(nanmean(abs(stat_sup_sacc))))
        plot(t_sacc,smooth(nanmean(abs(stat_sup_sacc)),'sgolay'),'--k','LineWidth', 2); % smoothed
        set(gca, 'xlim',[-0.2 0.2], 'ylim', ([0 6]), 'TickDir', 'out', 'FontSize', 18);
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
    case 'binomial_pb_dist_selected_signif'
        indx_exc = pop.indx_sel.(recArea).sacc.all.exc_signif;
        indx_sup = pop.indx_sel.(recArea).sacc.all.sup_signif;

        win = [-0.151 0.151]; 
        
        t_sacc = units(1).pro.neural.sacc.ts_pst;
        
        for i=1:length(indx_exc), stat_exc(i,:) = units(indx_exc(i)).stats.sacc.pval.pbDist_testStat; end
        for i=1:length(indx_sup), stat_sup(i,:) = units(indx_sup(i)).stats.sacc.pval.pbDist_testStat; end
        
         % plot exc
        figure; hold on;
        plot(t_sacc,(nanmean(abs(stat_exc))),'-k','LineWidth', 2);
        %plot(t_sacc,smooth(nanmean(abs(stat_exc)),'sgolay'),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0 4]), 'TickDir', 'out', 'FontSize', 18);
        title(['Sacc exc => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot sup
        plot(t_sacc,(nanmean(abs(stat_sup))),'--k','LineWidth', 2);
        %plot(t_sacc,smooth(nanmean(abs(stat_sup)),'sgolay'),'--k','LineWidth', 2); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0 4]), 'TickDir', 'out', 'FontSize', 18);
        title(['Exc and sup signif => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot all
        stat_all = ([abs(stat_exc) ; abs(stat_sup)]); 
        figure; hold on;
        %plot(t_sacc,(nanmean(stat_all)));
        plot(t_sacc,smooth(nanmean(stat_all),3),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Sacc all signif => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        %% plot data points per cell
        % find cells with two consecutive timepoints above or below 1.96
        % exc
        a0 = abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))); % input vector
        a0 = a0>1.96;
        for j = 1:length(a0(:,1))
            ii = strfind(a0(j,:),[1 1]);
            if ~isempty(ii)
                good_cell_sacc(j) = 1;
                good_cell_sacc_timepoints{j} = ii;
            else
                good_cell_sacc(j) = 0;
                good_cell_sacc_timepoints{j} = [];
            end
        end
        % sup
        a0_sup = abs(stat_sup(:,t_sacc>win(1) & t_sacc<win(2))); % input vector
        a0_sup = a0_sup>1.96;
        for j = 1:length(a0_sup(:,1))
            ii = strfind(a0_sup(j,:),[1 1]);
            if ~isempty(ii)
                good_cell_sacc_sup(j) = 1;
                good_cell_sacc_timepoints_sup{j} = ii;
            else
                good_cell_sacc_sup(j) = 0;
                good_cell_sacc_timepoints_sup{j} = [];
            end
        end
        
        
        stat_exc(isnan(stat_exc))=0; stat_sup(isnan(stat_sup))=0;
        % plot exc
        figure; hold on;
        plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))), '.', 'MarkerSize', 20, 'Color', [0.5 0.5 0.5]); hline(1.96, '--k');
       % plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))), 'MarkerSize', 20); hline(1.96);
        plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2))), '-', 'LineWidth', 2);
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 5],'ytick',[0 5], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; 
        % sum and median per bin
        figure; plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),nanmedian(abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 3],'ytick',[0 3], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; title('median exc'); hline(1.96,'--k')
        
        figure; plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),nansum(abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 30],'ytick',[0 30], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; title('sum exc');
        
        % all selected cells
        pcolor(t_sacc(t_sacc>win(1) & t_sacc<win(2)),1:size(abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2)))),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))));
        shading interp;
        set(gca, 'clim', [0 3]); % 3.92 max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar;
        
        % cells with two consecutive timepoints above or below 1.96
        pcolor(t_sacc(t_sacc>win(1) & t_sacc<win(2)),1:size(abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2)))),abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2))));
        shading interp;
        set(gca, 'clim', [0 3.92], 'yTick', [1 size(abs(stat_exc(logical(good_cell_sacc),t_sacc>win(1) & t_sacc<win(2))),1)]); % max is 1.96 x 2 = 3.92
        colormap([0 0 0; parula]);
        colorbar; title('exc')
        
        % plot sup
        figure; hold on;
        plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_sup(:,t_sacc>win(1) & t_sacc<win(2))), '.', 'MarkerSize', 20, 'Color', [0.5 0.5 0.5]); hline(1.96, '--k');
       % plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))), 'MarkerSize', 20); hline(1.96);
        plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2))), '-', 'LineWidth', 2);
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 5],'ytick',[0 5], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; hline(1.96,'--k')
         % sum and median per bin
        figure; plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),median(abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 3],'ytick',[0 3], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; title('median sup');hline(1.96,'--k')
        
        figure; plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),sum(abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 30],'ytick',[0 30], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; title('sum sup');
        
         % all selected cells
        pcolor(t_sacc(t_sacc>win(1) & t_sacc<win(2)),1:size(abs(stat_sup(:,t_sacc>win(1) & t_sacc<win(2)))),abs(stat_sup(:,t_sacc>win(1) & t_sacc<win(2))));
        shading interp;
        set(gca, 'clim', [0 3]); % 3.92 max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar;
        
         % cells with two consecutive timepoints above or below 1.96
        pcolor(t_sacc(t_sacc>win(1) & t_sacc<win(2)),1:size(abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2)))),abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2))));
        shading interp;
        set(gca, 'clim', [0 3.92], 'yTick', [1 size(abs(stat_sup(logical(good_cell_sacc_sup),t_sacc>win(1) & t_sacc<win(2))),1)]); % max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar; title('sup')
        
    case 'binomial_pb_dist_instr_selected_signif'
        indx_exc = pop.indx_sel.(recArea).instr_back.all.exc_signif;
        indx_sup = pop.indx_sel.(recArea).instr_back.all.sup_signif;

        win = [-0.3 0]; 
        
        t_instr = units(1).pro.neural.instrDir.ts_pst;
        
        for i=1:length(indx_exc), stat_exc(i,:) = units(indx_exc(i)).stats.instr_back.pval.pbDist_testStat; end
        for i=1:length(indx_sup), stat_sup(i,:) = units(indx_sup(i)).stats.instr_back.pval.pbDist_testStat; end
        
         % plot exc
        figure; hold on;
        plot(t_instr,(nanmean(abs(stat_exc))),'-k','LineWidth', 2);
        %plot(t_sacc,smooth(nanmean(abs(stat_exc)),'sgolay'),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0 4]), 'TickDir', 'out', 'FontSize', 18);
        title(['Instr exc => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot sup
        plot(t_instr,(nanmean(abs(stat_sup))),'--k','LineWidth', 2);
        %plot(t_sacc,smooth(nanmean(abs(stat_sup)),'sgolay'),'--k','LineWidth', 2); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0 4]), 'TickDir', 'out', 'FontSize', 18);
        title(['Exc and sup instr signif => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        % plot all
        stat_all = ([abs(stat_exc) ; abs(stat_sup)]); 
        figure; hold on;
        %plot(t_sacc,(nanmean(stat_all)));
        plot(t_instr,smooth(nanmean(stat_all),3),'LineWidth', 2,'Color','k'); % smoothed
        set(gca, 'xlim',[win(1) win(2)], 'ylim', ([0.6 2.5]), 'TickDir', 'out', 'FontSize', 18);
        title(['Instr all signif => ' recArea]); box off;
        xlabel('time(s)'); ylabel('Z-stat'); hline(1.96, '--k');
        
        %% plot data points per cell
        % find cells with two consecutive timepoints above or below 1.96
        % exc
        a0 = abs(stat_exc(:,t_instr>win(1) & t_instr<win(2))); % input vector
        a0 = a0>1.96;
        for j = 1:length(a0(:,1))
            ii = strfind(a0(j,:),[1 1]);
            if ~isempty(ii)
                good_cell_sacc(j) = 1;
                good_cell_sacc_timepoints{j} = ii;
            else
                good_cell_sacc(j) = 0;
                good_cell_sacc_timepoints{j} = [];
            end
        end
        % sup
        a0_sup = abs(stat_sup(:,t_instr>win(1) & t_instr<win(2))); % input vector
        a0_sup = a0_sup>1.96;
        for j = 1:length(a0_sup(:,1))
            ii = strfind(a0_sup(j,:),[1 1]);
            if ~isempty(ii)
                good_cell_sacc_sup(j) = 1;
                good_cell_sacc_timepoints_sup{j} = ii;
            else
                good_cell_sacc_sup(j) = 0;
                good_cell_sacc_timepoints_sup{j} = [];
            end
        end
        
        
        stat_exc(isnan(stat_exc))=0; stat_sup(isnan(stat_sup))=0;
        % plot exc
        figure; hold on;
        plot(t_instr(t_instr>win(1) & t_instr<win(2)),abs(stat_exc(:,t_instr>win(1) & t_instr<win(2))), '.', 'MarkerSize', 20, 'Color', [0.5 0.5 0.5]); hline(1.96, '--k');
       % plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))), 'MarkerSize', 20); hline(1.96);
        plot(t_instr(t_instr>win(1) & t_instr<win(2)),abs(stat_exc(logical(good_cell_sacc),t_instr>win(1) & t_instr<win(2))), '-', 'LineWidth', 2);
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 5],'ytick',[0 5], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; 
        % sum and median per bin
        figure; plot(t_instr(t_instr>win(1) & t_instr<win(2)),nanmedian(abs(stat_exc(logical(good_cell_sacc),t_instr>win(1) & t_instr<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 3],'ytick',[0 3], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; title('median exc'); hline(1.96,'--k')
        
        figure; plot(t_instr(t_instr>win(1) & t_instr<win(2)),nansum(abs(stat_exc(logical(good_cell_sacc),t_instr>win(1) & t_instr<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 30],'ytick',[0 30], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade exc'); box off; title('sum exc');
        
        % all selected cells
        pcolor(t_instr(t_instr>win(1) & t_instr<win(2)),1:size(abs(stat_exc(:,t_instr>win(1) & t_instr<win(2)))),abs(stat_exc(:,t_instr>win(1) & t_instr<win(2))));
        shading interp;
        set(gca, 'clim', [0 1.96]); % 3.92 max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar; title('exc')
        
        % cells with two consecutive timepoints above or below 1.96
        pcolor(t_instr(t_instr>win(1) & t_instr<win(2)),1:size(abs(stat_exc(logical(good_cell_sacc),t_instr>win(1) & t_instr<win(2)))),abs(stat_exc(logical(good_cell_sacc),t_instr>win(1) & t_instr<win(2))));
        shading interp;
        set(gca, 'clim', [0 3.92], 'yTick', [1 size(abs(stat_exc(logical(good_cell_sacc),t_instr>win(1) & t_instr<win(2))),1)]); % max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar; title('exc')
        
        % plot sup
        figure; hold on;
        plot(t_instr(t_instr>win(1) & t_instr<win(2)),abs(stat_sup(:,t_instr>win(1) & t_instr<win(2))), '.', 'MarkerSize', 20, 'Color', [0.5 0.5 0.5]); hline(1.96, '--k');
       % plot(t_sacc(t_sacc>win(1) & t_sacc<win(2)),abs(stat_exc(:,t_sacc>win(1) & t_sacc<win(2))), 'MarkerSize', 20); hline(1.96);
        plot(t_instr(t_instr>win(1) & t_instr<win(2)),abs(stat_sup(logical(good_cell_sacc_sup),t_instr>win(1) & t_instr<win(2))), '-', 'LineWidth', 2);
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 5],'ytick',[0 5], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; hline(1.96,'--k')
         % sum and median per bin
        figure; plot(t_instr(t_instr>win(1) & t_instr<win(2)),median(abs(stat_sup(logical(good_cell_sacc_sup),t_instr>win(1) & t_instr<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 3],'ytick',[0 3], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; title('median sup');hline(1.96,'--k')
        
        figure; plot(t_instr(t_instr>win(1) & t_instr<win(2)),sum(abs(stat_sup(logical(good_cell_sacc_sup),t_instr>win(1) & t_instr<win(2)))),'k')
        set(gca, 'xlim',[win(1) win(2)],'xTick',[-0.1 0 0.1], 'yLim', [0 30],'ytick',[0 30], 'TickDir', 'out', 'FontSize', 18);
        ylabel('z-stat'), xlabel('time (s)'); title('Saccade sup'); box off; title('sum sup');
        
         % all selected cells
        pcolor(t_instr(t_instr>win(1) & t_instr<win(2)),1:size(abs(stat_sup(:,t_instr>win(1) & t_instr<win(2)))),abs(stat_sup(:,t_instr>win(1) & t_instr<win(2))));
        shading interp;
        set(gca, 'clim', [0 1.96]); % 3.92 max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar; title('sup')
        
         % cells with two consecutive timepoints above or below 1.96
        pcolor(t_instr(t_instr>win(1) & t_instr<win(2)),1:size(abs(stat_sup(logical(good_cell_sacc_sup),t_instr>win(1) & t_instr<win(2)))),abs(stat_sup(logical(good_cell_sacc_sup),t_instr>win(1) & t_instr<win(2))));
        shading interp;
        set(gca, 'clim', [0 3.92], 'yTick', [1 size(abs(stat_sup(logical(good_cell_sacc_sup),t_instr>win(1) & t_instr<win(2))),1)]); % max is 1.96 x 2 
        colormap([0 0 0; parula]);
        colorbar; title('sup')
        
           
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
            align_sacc = t - sacc_onset-0.012; t_sacc = align_sacc(align_sacc >= -0.150 & align_sacc <= 0.150);
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
            align_sacc = t - sacc_onset-0.012; t_sacc = align_sacc(align_sacc >= -0.150 & align_sacc <= 0.150);
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
        
    case 'eye_move_single'
        % LOAD FILE SENT BY NICO
        
        t = units(cellNum).pro.timepointsEye/1000; % time
        
        % extract eye trace around saccade onset
        % pro
        c = 1;
        for i = 1:length(units(cellNum).pro.behav.trial)
            sacc_onset = units(cellNum).pro.behav.trial(i).saccadeOnset;
            align_sacc = t - sacc_onset-0.012; t_sacc = align_sacc(align_sacc >= -0.150 & align_sacc <= 0.150);
            x = units(cellNum).pro.behav.trial(i).eyePositionX; y = units(cellNum).pro.behav.trial(i).eyePositionY;
            if ~isnan(x)
                heye(c,:) = x(align_sacc >= -0.150 & align_sacc <= 0.150);
                veye(c,:) = y(align_sacc >= -0.150 & align_sacc <= 0.150);
                c = c+1;
            end
        end
        
        %                 for ii = 1:length(heye)
        %                 figure; hold on;
        %                 plot(t_sacc,smooth(heye(ii,:)',10),'b', 'LineWidth', 2);
        %                 plot(t_sacc,smooth(veye(ii,:)',10),'c', 'LineWidth', 2);
        %                 set(gca, 'yTick',[], 'TickDir', 'out', 'FontSize',22);
        %                 title(num2str(ii));
        %                 waitforbuttonpress; close all;
        %                 end
        
        % anti
        c = 1;
        for i = 1:length(units(cellNum).anti.behav.trial)
            sacc_onset = units(cellNum).anti.behav.trial(i).saccadeOnset;
            align_sacc = t - sacc_onset-0.012; t_sacc = align_sacc(align_sacc >= -0.150 & align_sacc <= 0.150);
            x = units(cellNum).anti.behav.trial(i).eyePositionX; y = units(cellNum).anti.behav.trial(i).eyePositionY;
            if ~isnan(x)
                heye(c,:) = x(align_sacc >= -0.150 & align_sacc <= 0.150);
                veye(c,:) = y(align_sacc >= -0.150 & align_sacc <= 0.150);
                c = c+1;
            end
        end
        
        for ii = 1:length(heye)
            figure; hold on;
            plot(t_sacc,smooth(heye(ii,:)',10),'b', 'LineWidth', 2);
            plot(t_sacc,smooth(veye(ii,:)',10),'c', 'LineWidth', 2);
            set(gca, 'yTick',[], 'TickDir', 'out', 'FontSize',22);
            title([num2str(ii) ' anti']);
            waitforbuttonpress; close all;
        end
        %
        % print('eye_trace','-depsc2', '-painters', '-cmyk')
        
    case 'eye_move_instr'
        % LOAD FILE SENT BY NICO
        
        t = units(cellNum).pro.timepointsEye-0.012/1000; % time
        t_instr = t(t>-0.150 & t<0.300);
        % extract eye trace
        % pro
        for i = 1:length(units(cellNum).pro.behav.trial)
            x_pro(i,:) = units(cellNum).pro.behav.trial(i).eyePositionX; y_pro(i,:) = units(cellNum).pro.behav.trial(i).eyePositionY;
            heye_pro(i,:) = x_pro(t>-0.150 & t<0.300); veye_pro(i,:) = y_pro(t>-0.150 & t<0.300);
        end
        
        % plot all
        figure; hold on;
        plot(t, x_pro, 'm'); plot(t, y_pro, 'k');
        set(gca, 'xlim', [0 0.350],'ylim',[-5 5], 'TickDir', 'out', 'FontSize',18);
        title(['Eye trace anti cell ' num2str(cellNum)])
        
        % print('eye_trace','-depsc2', '-painters', '-cmyk')
        
        % anti
        for i = 1:length(units(cellNum).pro.behav.trial)
            x_anti(i,:) = units(cellNum).pro.behav.trial(i).eyePositionX; y_anti(i,:) = units(cellNum).pro.behav.trial(i).eyePositionY;
            heye(i,:) = x_anti(t>-0.150 & t<0.300); veye(i,:) = y_anti(t>-0.150 & t<0.300);
        end
        
        % plot all
        figure; hold on;
        plot(t, x_anti, 'b'); plot(t, y_anti, 'c');
        set(gca, 'xlim', [0 0.350],'ylim',[-5 5], 'TickDir', 'out', 'FontSize',18);
        title(['Eye trace anti cell ' num2str(cellNum)])
        
        % plot single
        for ii = 1:length(x)
            figure; hold on;
            plot(t,smooth(x(ii,:)',10),'b', 'LineWidth', 2);
            plot(t,smooth(y(ii,:)',10),'c', 'LineWidth', 2);
            set(gca, 'yTick',[], 'TickDir', 'out', 'FontSize',22);
            title([num2str(ii) ' anti']);
            %waitforbuttonpress; close all;
        end
        
        % plot all
        plot(t, x);
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
        set(gca, 'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', 'FontSize', 30); title('Modulation ratio OMV vs Lateral');
        axis square;
        
        
        % for significantly different
        
        for i = 1:length(indx_area_vermis)
            indx_sign_vermis(i) = logical(units(indx_area_vermis(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk);
        end
        unit_vermis_sign = unit_vermis(indx_sign_vermis);
        unit_vermis_ns = unit_vermis(~indx_sign_vermis);
        
        for i = 1:length(indx_area_lat)
            indx_sign_lat(i) = logical(units(indx_area_lat(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk);
        end
        unit_lat_sign = unit_lat(indx_sign_lat);
        unit_lat_ns = unit_lat(~indx_sign_lat);
        
        % extract
        for i = 1:length(unit_vermis_sign)
            mod_r_vermis_sign(i) = unit_vermis_sign(i).stats.mod_ratio;
        end
        for i = 1:length(unit_vermis_ns)
            mod_r_vermis_ns(i) = unit_vermis_ns(i).stats.mod_ratio;
        end
        
        for i = 1:length(unit_lat_sign)
            mod_r_lat_sign(i) = unit_lat_sign(i).stats.mod_ratio;
        end
        mod_r_lat_sign(6)=[];
        for i = 1:length(unit_lat_ns)
            mod_r_lat_ns(i) = unit_lat_ns(i).stats.mod_ratio;
        end
        
        
        figure; hold on;
        cdfplot(mod_r_vermis_sign);
        cdfplot(mod_r_lat_sign);
        vline(1);
        set(gca, 'yTick',[0 0.5 1],'xlim',[0.5 1.5] ,'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out', 'FontSize', 30); title('Modulation ratio OMV vs Lateral signif');
        axis square;
        
        [h_mod,p_mod] = kstest2(mod_r_vermis_sign, mod_r_lat_sign);
        
        figure; hold on;
        [f_v,x_v] = ksdensity(mod_r_vermis_sign);
        y1_v = cumsum(f_v); y1_v = y1_v./max(y1_v);
        plot(x_v,y1_v,'k','linewidth',2);
        xlabel('Modulation ratio'); axis square; box off;
        
        [f_l,x_l] = ksdensity(mod_r_lat_sign);
        y1_l = cumsum(f_l); y1_l = y1_l./max(y1_l);
        plot(x_l,y1_l,'--k','linewidth',2);
        set(gca, 'yTick',[0 50 100], 'xlim',[0.25 1.6], 'xTick', [0.5 1 1.5], 'TickDir', 'out', 'FontSize', 30);
        xlabel('Modulation ratio'); axis square; box off; vline(1, 'k');
        
        % not significant
        [f_ns_v,x_ns_v] = ksdensity(mod_r_vermis_ns);
        y1_ns_v = cumsum(f_ns_v); y1_ns_v = y1_ns_v./max(y1_ns_v);
        plot(x_ns_v,y1_ns_v,'r','linewidth',2);
        xlabel('Modulation ratio'); axis square; box off;
        
        [f_ns_l,x_ns_l] = ksdensity(mod_r_lat_ns);
        y1_ns_l = cumsum(f_ns_l); y1_ns_l = y1_ns_l./max(y1_ns_l);
        plot(x_ns_l,y1_ns_l,'g','linewidth',2);
        xlabel('Modulation ratio'); axis square; box off;
        
        % stat
        [h,p] = kstest2(x_v, x_l)
        
        % get slope
        [coefficients_omv,s_omv] = polyfit(x_v, y1_v, 1);
        [coefficients_lat,s_lat] = polyfit(x_l, y1_l, 1);
        % Now get the slope, which is the first coefficient in the array:
        slope_omv = coefficients_omv(1);
        slope_lat = coefficients_lat(1);
        [~,SE_omv] = polyval(coefficients_omv,x_v,s_omv);
        [~,SE_lat] = polyval(coefficients_lat,x_l,s_lat);
        
        % caluclate normal z as we cannot assume homogeneity of the error
        % variances in the two samples. Difference between the two slopes
        % divided by the standard error of the difference between the
        % slopes.
        
        z_slope = (coefficients_lat(1)-coefficients_omv(1))/(sqrt(SE_lat(2)-SE_omv(2)))
        
    case 'mod_ratio_selected'
        indx_exc_all_omv = [pop.indx_sel.vermis.sacc.all.pro.exc pop.indx_sel.vermis.sacc.all.anti.exc];
        indx_sup_all_omv = [pop.indx_sel.vermis.sacc.all.pro.sup pop.indx_sel.vermis.sacc.all.anti.sup];
        
        indx_exc_all_lat = [pop.indx_sel.lateral.sacc.all.pro.exc pop.indx_sel.lateral.sacc.all.anti.exc];
        indx_sup_all_lat = [pop.indx_sel.lateral.sacc.all.pro.sup pop.indx_sel.lateral.sacc.all.anti.sup];
        indx_sup_all_lat(3)=[];
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=0 & t<=0.150);
        for i=1:length(indx_exc_all_omv),mod_exc_omv(i,:) = units(indx_exc_all_omv(i)).stats.mod_ratio; mod_exc_omv_mu(i) = mean(mod_exc_omv(i,:)); end
        for i=1:length(indx_sup_all_omv),mod_sup_omv(i,:) = units(indx_sup_all_omv(i)).stats.mod_ratio; mod_sup_omv_mu(i) = mean(mod_sup_omv(i,:)); end
        for i=1:length(indx_exc_all_lat),mod_exc_lat(i,:) = units(indx_exc_all_lat(i)).stats.mod_ratio; mod_exc_lat_mu(i) = mean(mod_exc_lat(i,:)); end
        for i=1:length(indx_sup_all_lat),mod_sup_lat(i,:) = units(indx_sup_all_lat(i)).stats.mod_ratio; mod_sup_lat_mu(i) = mean(mod_sup_lat(i,:)); end
       
        
        % plot over time
        figure; hold on;
        plot(t_win, mean(mod_exc_omv), '-k');
        plot(t_win, mean(mod_exc_lat), '--k');
        set(gca, 'xlim',[0 0.151], 'ylim', [0.9 1.15], 'ytick', [0.9 1.15],'TickDir', 'out', 'FontSize', 18);
        hline(1,'k'); title('Mod ratio exc omv vs lat')
        axis square;
        
        figure; hold on 
        plot(t_win, mean(mod_sup_omv), '-k');
        plot(t_win, mean(mod_sup_lat), '--k');
        set(gca, 'xlim',[0 0.151], 'ylim', [0.9 1.15], 'ytick', [0.9 1.15],'TickDir', 'out', 'FontSize', 18);
        hline(1,'k'); title('Mod ratio sup omv vs lat')
        axis square;
        
        % cdf exc
        figure; hold on;
        [f_v,x_v] = ksdensity(mod_exc_omv_mu);
        y1_v = cumsum(f_v); y1_v = y1_v./max(y1_v);
        plot(x_v,y1_v,'-k','linewidth',2);
        xlabel('Modulation ratio'); axis square; box off;
        
        [f_l,x_l] = ksdensity(mod_exc_lat_mu);
        y1_v = cumsum(f_l); y1_v = y1_v./max(y1_v);
        plot(x_l,y1_v,'--k','linewidth',2);
        set(gca, 'xlim',[0.2 1.52],'xTick', [0.25 0.5 1], 'ylim', [0 1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        title('mod ratio pro/anti exc omv vs lat'); vline(1,'k')
        
        [~,p_exc] = kstest2(x_v, x_l)
%         histogram(mod_exc_omv_mu,25); hold on
%         histogram(mod_exc_lat_mu,25);
        
         % cdf sup
        figure; hold on;
        [f_v,x_v] = ksdensity(mod_sup_omv_mu);
        y1_v = cumsum(f_v); y1_v = y1_v./max(y1_v);
        plot(x_v,y1_v,'-k','linewidth',2);
        xlabel('Modulation ratio'); axis square; box off;
        
        [f_l,x_l] = ksdensity(mod_sup_lat_mu);
        y1_l = cumsum(f_l); y1_l = y1_l./max(y1_l);
        plot(x_l,y1_l,'--k','linewidth',2);
        set(gca, 'xlim',[0.5 1.62],'xTick', [0.5 1 1.62], 'ylim', [0 1], 'ytick', [0 0.5 1],'TickDir', 'out', 'FontSize', 18);
        title('mod ratio pro/anti sup omv vs lat'); vline(1,'k')
        
        [~,p_sup] = kstest2(x_v, x_l)
         histogram(mod_sup_omv_mu,25); hold on
         histogram(mod_sup_lat_mu,25);
        
    case 'sorted_colormap_sacc'
        % get area
        for cellNum = 1:length(units)
            indx(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx = find(indx);
        nunits_area = 1:length(indx);
        
        win = [-0.251 0.251];
        
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=win(1) & t<=win(2));
        
        for j=1:length(indx)
            r_pro(j,:) = units(indx(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); % psth
            r_anti(j,:) = units(indx(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2));
        end
        % pro
        [maxRates,pos_max] = max(r_pro, [], 2);
        [~,indx_max] = sort(pos_max);
        r_pro_norm = r_pro./repmat(maxRates,[1 size(r_pro,2)]);
        r_pro_sorted = r_pro_norm(indx_max,:); [~, max_pos_pro] =  max(r_pro_sorted,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); % %colormap(flipud(B'));
        imagesc(t_win,1:size(r_pro_sorted,1),r_pro_sorted, [0 1]); hold on;
        scatter(t_win(max_pos_pro),1:size(r_pro_sorted,1),20,'k','filled');
        %plot(t(max_pos),1:size(r_pro_sorted,1),'w', 'LineWidth',2);
        %scatter(indx_max,1:size(r_pro_sorted,1),5,'k','filled');
        set(gca,'xlim',[win(1) win(2)], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted Pro ' recArea])
        %figure; histogram(t(max_pos_pro),t); set(gca, 'yTick', [0 17],'xTick', [], 'TickDir', 'out', 'FontSize', 18); box off;
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t(max_pos_pro),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro');
        
        % colormap unsorted
        %         figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); %colormap(flipud(B'));
        %         imagesc(t,1:size(r_pro_norm,1),r_pro_norm, [0 1]);
        %         set(gca,'xlim',[-0.15 0.151],'YTickLabel', [],'TickDir', 'out', 'FontSize', 18); box off;
        %         vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Pro ' recArea])
        
        % anti
        [maxRates_anti,pos_max_anti ] = max(r_anti, [], 2);
        [~,indx_max_anti] = sort(pos_max_anti);
        r_anti_norm = r_anti./repmat(maxRates_anti,[1 size(r_anti,2)]);
        r_anti_sorted = r_anti_norm(indx_max_anti,:); [~, max_pos_anti] =  max(r_anti_sorted,[],2);
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); %colormap(B'); %colormap(flipud(B'));
        imagesc(t_win,1:size(r_anti_sorted,1),r_anti_sorted, [0 1]); hold on;
        scatter(t_win(max_pos_anti),1:size(r_anti_sorted,1),20,'k','filled');
        %plot(t(max_pos),1:size(r_anti_sorted,1),'w', 'LineWidth',2);
        %scatter(indx_max,1:size(r_pro_sorted,1),5,'k','filled');
        set(gca,'xlim',[win(1) win(2)], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted anti ' recArea])
        %figure; histogram(t(max_pos_anti),t); set(gca, 'yTick', [0 14], 'TickDir', 'out', 'FontSize', 18); box off;
        
        %anti
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_anti),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti');
        
        %% plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(max_pos_pro), 'm', 'LineWidth',2);
        plot(t_win(max_pos_anti), 'b', 'LineWidth',2);
        set(gca,'ylim',[win(1) win(2)],'xlim', [0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square;
        
        % compute difference between both - zero no change
        figure; hold on;
        plot(t(max_pos_pro)-t(max_pos_anti), 'k','LineWidth',2); %% mean(t(max_pos_pro)-t(max_pos_anti)) % for difference between both
        set(gca,'ylim',[win(1) win(2)], 'xlim', [0 90], 'TickDir', 'out', 'FontSize', 18); box off;
        axis square; hline(0); vline(72); 
        
        
        %% plot a few neurons for sanity check
        cell = 46;
        t_check = units(cell).pro.neural.sacc.ts_pst;
        r_check_pro = units(cell).pro.neural.sacc.rate_pst;
        r_check_anti = units(cell).anti.neural.sacc.rate_pst;
        
        figure;
        plot(t_check,r_check_pro,'color','m'); hold on;
        plot(t_check,r_check_anti,'color','b');
        set(gca, 'xlim',([-0.150 0.151]), 'TickDir', 'out', 'FontSize',22); % analysis window size
        xlabel('Time (s)'); ylabel ('Firing rate (spk/s)');
        vline(0, 'k-');
        box off
        
         case 'sorted_colormap_sacc_selected'
        % get selected
        
        indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc; 
        indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc; 
        indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup;
        indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup; 
        
        shared_exc = intersect(indx_exc_pro, indx_exc_anti); shared_sup = intersect(indx_sup_pro, indx_sup_anti);
        shared_exc_pro = ismember(indx_exc_pro, shared_exc); shared_exc_anti = ismember(indx_exc_anti, shared_exc); 
        shared_sup_pro = ismember(indx_sup_pro, shared_sup); shared_sup_anti = ismember(indx_sup_anti, shared_sup); 
 
        win = [-0.151 0.251];
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=win(1) & t<=win(2));
        
        for j=1:length(indx_exc_pro), r_pro_exc(j,:) = units(indx_exc_pro(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        for j=1:length(indx_exc_anti), r_anti_exc(j,:) = units(indx_exc_anti(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        
        for j=1:length(indx_sup_pro), r_pro_sup(j,:) = units(indx_sup_pro(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        for j=1:length(indx_sup_anti), r_anti_sup(j,:) = units(indx_sup_anti(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
       
        %% 
        % pro exc
        [maxRates,pos_max_pro] = max(r_pro_exc, [], 2);
        [~,indx_max_pro] = sort(pos_max_pro);
        r_pro_norm = r_pro_exc./repmat(maxRates,[1 size(r_pro_exc,2)]);
        r_pro_sorted_exc = r_pro_norm(indx_max_pro,:); [~, max_pos_pro] =  max(r_pro_sorted_exc,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_pro_sorted_exc,1),r_pro_sorted_exc, [0 1]); hold on;
        scatter(t_win(max_pos_pro),1:size(r_pro_sorted_exc,1),40,'k','filled');
        scatter(t_win(max_pos_pro(shared_exc_pro)),find(shared_exc_pro), 40, 'c', 'filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Pro both ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_pro),1:size(r_pro_exc,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_exc_pro)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro exc');
        
        % anti exc
        [maxRates_anti,pos_max_anti] = max(r_anti_exc, [], 2);
        [~,indx_max_anti] = sort(pos_max_anti);
        r_anti_norm = r_anti_exc./repmat(maxRates_anti,[1 size(r_anti_exc,2)]);
        r_anti_sorted_exc = r_anti_norm(indx_max_anti,:); [~, max_pos_anti] =  max(r_anti_sorted_exc,[],2);
        % sort based on prosaccade
        % r_anti_pro_sorted_exc = r_anti_norm(indx_max_pro,:);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_anti_sorted_exc,1),r_anti_sorted_exc, [0 1]); hold on;
        scatter(t_win(max_pos_anti),1:size(r_anti_sorted_exc,1),40,'k','filled');
        scatter(t_win(max_pos_anti(shared_exc_anti)),find(shared_exc_anti), 40, 'c', 'filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Anti both ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_anti),1:size(r_anti_exc,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_exc_anti)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti exc');
        
        % plot anti with pro sorting
%         B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
%         figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
%         imagesc(t_win,1:size(r_anti_pro_sorted_exc,1),r_anti_pro_sorted_exc, [0 1]); hold on;
%         scatter(t_win(max_pos_pro),1:size(r_anti_pro_sorted_exc,1),40,'k','filled');
%         set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
%         vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Anti both Pro sorted ' recArea])
        
        % plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(max_pos_pro), 'm', 'LineWidth',2);
        plot(t_win(max_pos_anti), 'b', 'LineWidth',2);
        set(gca,'ylim',[0 0.250],'xlim', [1 length(max_pos_anti)], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti exc ' recArea]);
        axis square;
        
        % plot without sorting
%         figure; hold on;
%         plot(t_win(pos_max_pro), '.m', 'LineWidth',2, 'MarkerSize', 30);
%         plot(t_win(pos_max_anti), '.b', 'LineWidth',2, 'MarkerSize', 30);
%         plot(t_win(pos_max_pro), 'm', 'LineWidth',2, 'MarkerSize', 30);
%         plot(t_win(pos_max_anti), 'b', 'LineWidth',2, 'MarkerSize', 30);
%         set(gca,'ylim',[0 0.250],'xlim', [1 14], 'TickDir', 'out', 'FontSize', 18); box off;
%         ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti exc ' recArea]);
%         axis square;
        
       % plot scatter for shared ones
        pro = t_win(max_pos_pro(shared_exc_pro)); anti = t_win(max_pos_anti(shared_exc_anti)); 
        figure; hold on;
        plot(pro,anti, '.k','MarkerSize', 30); hold on; plot([0 0.25],[0 0.25])
        set(gca, 'xlim', [0 0.25],'xTick', [0 0.25], 'ylim', [0 0.25], 'yTick', [0 0.25], 'TickDir', 'out', 'FontSize', 22, 'xlabel', []); box off;
        xlabel('pro'); ylabel('anti'); title(['exc ' recArea]); axis square

        % compute difference between both - zero no change
%         figure; hold on;
%         plot(t_win(max_pos_pro)-t_win(max_pos_anti), 'k','LineWidth',2); %% mean(t(max_pos_pro)-t(max_pos_anti)) % for difference between both
%         set(gca,'ylim',[0 0.250],'xlim', [0 17], 'TickDir', 'out', 'FontSize', 18); box off;
%         axis square; hline(0); vline(68);
        
        %% pro sup
        [maxRates_sup_pro,pos_max] = max(r_pro_sup, [], 2); [~,pos_min] = min(r_pro_sup, [], 2);
        [~,indx_max] = sort(pos_max); [~,indx_min] = sort(pos_min);
        r_pro_norm = r_pro_sup./repmat(maxRates_sup_pro,[1 size(r_pro_sup,2)]);
        r_pro_sorted_sup = r_pro_norm(indx_min,:); [~, min_pos_pro] =  min(r_pro_sorted_sup,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_pro_sorted_sup,1),r_pro_sorted_sup, [0 1]); hold on;
        scatter(t_win(min_pos_pro),1:size(r_pro_sorted_sup,1),40,'k','filled');
        scatter(t_win(min_pos_pro(shared_sup_pro)),find(shared_sup_pro), 40, 'c', 'filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted Pro sup ' recArea])
       
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(min_pos_pro),1:size(r_pro_sup,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_min)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro sup');
        
        % anti sup
        [maxRates,anti_max] = max(r_anti_sup, [], 2); [~,pos_min] = min(r_anti_sup, [], 2);
        [~,indx_max] = sort(pos_max); [~,indx_min] = sort(pos_min);
        r_anti_norm = r_anti_sup./repmat(maxRates,[1 size(r_anti_sup,2)]);
        r_anti_sorted_sup = r_anti_norm(indx_min,:); [~, min_pos_anti] =  min(r_anti_sorted_sup,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_anti_sorted_sup,1),r_anti_sorted_sup, [0 1]); hold on;
        scatter(t_win(min_pos_anti),1:size(r_anti_sorted_sup,1),40,'k','filled');
        scatter(t_win(min_pos_anti(shared_sup_anti)),find(shared_sup_anti), 40, 'c', 'filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted Anti sup ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(min_pos_anti),1:size(r_anti_sup,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_sup_anti)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti sup');
        
        % plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(min_pos_pro), 'm', 'LineWidth',2);
        plot(t_win(min_pos_anti), 'b', 'LineWidth',2);
        set(gca,'ylim',[0 0.250],'xlim', [1 length(min_pos_pro)], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti sup ' recArea]);
        axis square;
        
        % plot scatter for shared ones
        figure; hold on;
        plot(t_win(min_pos_pro(shared_sup_pro)),t_win(min_pos_anti(shared_sup_anti)), '.k','MarkerSize', 30); hold on; plot([0 0.25],[0 0.25])
        set(gca, 'xlim', [0 0.25],'xTick', [0 0.25], 'ylim', [0 0.25], 'yTick', [0 0.25], 'TickDir', 'out', 'FontSize', 22, 'xlabel', []); box off;
        xlabel('pro'); ylabel('anti'); title(['sup ' recArea]); axis square
        
         case 'sorted_colormap_sacc_selected_both_lateral'
        % get selected
        
%         indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc; 
%         indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup;
%         indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
%         indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
        
        indx_exc_both = [pop.indx_sel.(recArea).sacc.both.exc pop.indx_sel.(recArea).sacc.both.mixed]; % Both lateral
        indx_sup_both = [pop.indx_sel.(recArea).sacc.both.sup pop.indx_sel.(recArea).sacc.both.mixed]; % Both lateral
        
        win = [-0.151 0.251];
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=win(1) & t<=win(2));
        
        
        for j=1:length(indx_exc_both), r_pro_exc(j,:) = units(indx_exc_both(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end % psth
        for j=1:length(indx_exc_both(1:end-2)), r_anti_exc(j,:) = units(indx_exc_both(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end % no 33 or 152
   
        
        for j=1:length(indx_sup_both(1:end-2)), r_pro_sup(j,:) = units(indx_sup_both(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end % psth no 33 or 152
        for j=1:length(indx_sup_both), r_anti_sup(j,:) = units(indx_sup_both(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end % psth

        % pro exc
        [maxRates,pos_max_pro] = max(r_pro_exc, [], 2);
        [~,indx_max_pro] = sort(pos_max_pro);
        r_pro_norm = r_pro_exc./repmat(maxRates,[1 size(r_pro_exc,2)]);
        r_pro_sorted_exc = r_pro_norm(indx_max_pro,:); [~, max_pos_pro] =  max(r_pro_sorted_exc,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_pro_sorted_exc,1),r_pro_sorted_exc, [0 1]); hold on;
        scatter(t_win(max_pos_pro),1:size(r_pro_sorted_exc,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Pro both ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_pro),1:size(r_pro_exc,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_exc_both)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro exc');
        
        % anti exc
        [maxRates_anti,pos_max_anti] = max(r_anti_exc, [], 2);
        [~,indx_max_anti] = sort(pos_max_anti);
        r_anti_norm = r_anti_exc./repmat(maxRates_anti,[1 size(r_anti_exc,2)]);
        r_anti_sorted_exc = r_anti_norm(indx_max_anti,:); [~, max_pos_anti] =  max(r_anti_sorted_exc,[],2);
        % sort based on prosaccade
        % r_anti_pro_sorted_exc = r_anti_norm(indx_max_pro,:);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_anti_sorted_exc,1),r_anti_sorted_exc, [0 1]); hold on;
        scatter(t_win(max_pos_anti),1:size(r_anti_sorted_exc,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Anti both ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_anti),1:size(r_anti_exc,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_exc_both)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti exc');
        
        % plot anti with pro sorting
%         B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
%         figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
%         imagesc(t_win,1:size(r_anti_pro_sorted_exc,1),r_anti_pro_sorted_exc, [0 1]); hold on;
%         scatter(t_win(max_pos_pro),1:size(r_anti_pro_sorted_exc,1),40,'k','filled');
%         set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
%         vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Anti both Pro sorted ' recArea])
        
        % plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(max_pos_pro(1:12)), 'm', 'LineWidth',2);
        plot(t_win(max_pos_anti), 'b', 'LineWidth',2);
        set(gca,'ylim',[0 0.250],'xlim', [1 length(max_pos_pro(1:12))], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti exc ' recArea]);
        axis square;
        
        % plot without sorting
        figure; hold on;
        plot(t_win(pos_max_pro(1:12)), '.m', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_anti), '.b', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_pro(1:12)), 'm', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_anti), 'b', 'LineWidth',2, 'MarkerSize', 30);
        set(gca,'ylim',[0 0.250],'xlim', [1 14], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti exc ' recArea]);
        axis square;
        
        %scatter
        figure; hold on;
        plot(t_win(pos_max_pro(1:12)),t_win(pos_max_anti), '.k','MarkerSize', 30); hold on; plot([0 0.25],[0 0.25], 'k')
        set(gca, 'xlim', [0 0.25],'xTick', [0 0.25], 'ylim', [0 0.25], 'yTick', [0 0.25], 'TickDir', 'out', 'FontSize', 22, 'xlabel', []); box off;
        xlabel('pro'); ylabel('anti'); axis square

        % compute difference between both - zero no change
%         figure; hold on;
%         plot(t_win(max_pos_pro)-t_win(max_pos_anti), 'k','LineWidth',2); %% mean(t(max_pos_pro)-t(max_pos_anti)) % for difference between both
%         set(gca,'ylim',[0 0.250],'xlim', [0 17], 'TickDir', 'out', 'FontSize', 18); box off;
%         axis square; hline(0); vline(68);
        
        %% pro sup
        [maxRates_sup_pro,pos_max] = max(r_pro_sup, [], 2); [~,pos_min] = min(r_pro_sup, [], 2);
        [~,indx_max] = sort(pos_max); [~,indx_min] = sort(pos_min);
        r_pro_norm = r_pro_sup./repmat(maxRates_sup_pro,[1 size(r_pro_sup,2)]);
        r_pro_sorted_sup = r_pro_norm(indx_min,:); [~, min_pos_pro] =  min(r_pro_sorted_sup,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_pro_sorted_sup,1),r_pro_sorted_sup, [0 1]); hold on;
        scatter(t_win(min_pos_pro),1:size(r_pro_sorted_sup,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted Pro sup ' recArea])
       
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(min_pos_pro),1:size(r_pro_sup,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_sup_both)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro sup');
        
        % anti sup
        [maxRates,anti_max] = max(r_anti_sup, [], 2); [~,pos_min] = min(r_anti_sup, [], 2);
        [~,indx_max] = sort(pos_max); [~,indx_min] = sort(pos_min);
        r_anti_norm = r_anti_sup./repmat(maxRates,[1 size(r_anti_sup,2)]);
        r_anti_sorted_sup = r_anti_norm(indx_min,:); [~, min_pos_anti] =  min(r_anti_sorted_sup,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_anti_sorted_sup,1),r_anti_sorted_sup, [0 1]); hold on;
        scatter(t_win(min_pos_anti),1:size(r_anti_sorted_sup,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted Anti sup ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(min_pos_anti),1:size(r_anti_sup,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_sup_both)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti sup');
        
        % plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(min_pos_pro), 'm', 'LineWidth',2);
        plot(t_win(min_pos_anti(1:9)), 'b', 'LineWidth',2);
        set(gca,'ylim',[0 0.250],'xlim', [1 length(min_pos_pro)], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti sup ' recArea]);
        axis square;
        
         %scatter
        figure; hold on;
        plot(t_win(min_pos_pro),t_win(min_pos_anti(1:9)), '.k','MarkerSize', 30); hold on; plot([0 0.25],[0 0.25])
        set(gca, 'xlim', [0 0.25],'xTick', [0 0.25], 'ylim', [0 0.25], 'yTick', [0 0.25], 'TickDir', 'out', 'FontSize', 22, 'xlabel', []); box off;
        xlabel('pro'); ylabel('anti'); axis square
        
         case 'sorted_colormap_sacc_selected_both'
        % get selected
        
        indx_exc = pop.indx_sel.(recArea).sacc.both.exc;
        indx_sup = pop.indx_sel.(recArea).sacc.both.sup;
        
        win = [-0.151 0.251];
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=win(1) & t<=win(2));
        
        
        for j=1:length(indx_exc)
            r_pro_exc(j,:) = units(indx_exc(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2));
            r_anti_exc(j,:) = units(indx_exc(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2));
        end
        
        for j=1:length(indx_sup)
            r_pro_sup(j,:) = units(indx_sup(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2));
            r_anti_sup(j,:) = units(indx_sup(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2));
        end


        % pro exc
        [maxRates,pos_max_pro] = max(r_pro_exc, [], 2);
        [~,indx_max_pro] = sort(pos_max_pro);
        r_pro_norm = r_pro_exc./repmat(maxRates,[1 size(r_pro_exc,2)]);
        r_pro_sorted_exc = r_pro_norm(indx_max_pro,:); [~, max_pos_pro] =  max(r_pro_sorted_exc,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_pro_sorted_exc,1),r_pro_sorted_exc, [0 1]); hold on;
        scatter(t_win(max_pos_pro),1:size(r_pro_sorted_exc,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Pro both ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_pro),1:size(r_pro_exc,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_exc)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro exc');
        
        % anti exc
        [maxRates_anti,pos_max_anti] = max(r_anti_exc, [], 2);
        [~,indx_max_anti] = sort(pos_max_anti);
        r_anti_norm = r_anti_exc./repmat(maxRates_anti,[1 size(r_anti_exc,2)]);
        r_anti_sorted_exc = r_anti_norm(indx_max_anti,:); [~, max_pos_anti] =  max(r_anti_sorted_exc,[],2);
        % sort based on prosaccade
        % r_anti_pro_sorted_exc = r_anti_norm(indx_max_pro,:);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_anti_sorted_exc,1),r_anti_sorted_exc, [0 1]); hold on;
        scatter(t_win(max_pos_anti),1:size(r_anti_sorted_exc,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Anti both ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_anti),1:size(r_anti_exc,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_exc)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti exc');
        
        % plot anti with pro sorting
%         B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
%         figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
%         imagesc(t_win,1:size(r_anti_pro_sorted_exc,1),r_anti_pro_sorted_exc, [0 1]); hold on;
%         scatter(t_win(max_pos_pro),1:size(r_anti_pro_sorted_exc,1),40,'k','filled');
%         set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
%         vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Anti both Pro sorted ' recArea])
%         
        % plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(max_pos_pro), 'm', 'LineWidth',2);
        plot(t_win(max_pos_anti), 'b', 'LineWidth',2);
        set(gca,'ylim',[0 0.250],'xlim', [1 length(max_pos_anti)], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti exc ' recArea]);
        axis square;
        
        % plot without sorting
        figure; hold on;
        plot(t_win(pos_max_pro), '.m', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_anti), '.b', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_pro), 'm', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_anti), 'b', 'LineWidth',2, 'MarkerSize', 30);
        set(gca,'ylim',[0 0.250],'xlim', [1 14], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti exc ' recArea]);
        axis square;
        
        %scatter
        figure; hold on;
        plot(t_win(max_pos_pro),t_win(max_pos_anti), '.k','MarkerSize', 30); hold on; plot([0 0.25],[0 0.25])
        set(gca, 'xlim', [0 0.25],'xTick', [0 0.25], 'ylim', [0 0.25], 'yTick', [0 0.25], 'TickDir', 'out', 'FontSize', 22, 'xlabel', []); box off;
        xlabel('pro'); ylabel('anti'); axis square

        % compute difference between both - zero no change
%         figure; hold on;
%         plot(t_win(max_pos_pro)-t_win(max_pos_anti), 'k','LineWidth',2); %% mean(t(max_pos_pro)-t(max_pos_anti)) % for difference between both
%         set(gca,'ylim',[0 0.250],'xlim', [0 17], 'TickDir', 'out', 'FontSize', 18); box off;
%         axis square; hline(0); vline(68);
        
        %% pro sup
        [maxRates_sup_pro,pos_max] = max(r_pro_sup, [], 2); [~,pos_min] = min(r_pro_sup, [], 2);
        [~,indx_max] = sort(pos_max); [~,indx_min] = sort(pos_min);
        r_pro_norm = r_pro_sup./repmat(maxRates_sup_pro,[1 size(r_pro_sup,2)]);
        r_pro_sorted_sup = r_pro_norm(indx_min,:); [~, min_pos_pro] =  min(r_pro_sorted_sup,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_pro_sorted_sup,1),r_pro_sorted_sup, [0 1]); hold on;
        scatter(t_win(min_pos_pro),1:size(r_pro_sorted_sup,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted Pro sup ' recArea])
       
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(min_pos_pro),1:size(r_pro_sup,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_sup)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro sup');
        
        % anti sup
        [maxRates,anti_max] = max(r_anti_sup, [], 2); [~,pos_min] = min(r_anti_sup, [], 2);
        [~,indx_max] = sort(pos_max); [~,indx_min] = sort(pos_min);
        r_anti_norm = r_anti_sup./repmat(maxRates,[1 size(r_anti_sup,2)]);
        r_anti_sorted_sup = r_anti_norm(indx_min,:); [~, min_pos_anti] =  min(r_anti_sorted_sup,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_anti_sorted_sup,1),r_anti_sorted_sup, [0 1]); hold on;
        scatter(t_win(min_pos_anti),1:size(r_anti_sorted_sup,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc Sorted Anti sup ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(min_pos_anti),1:size(r_anti_sup,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_sup)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti sup');
        
        % plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(min_pos_pro), 'm', 'LineWidth',2);
        plot(t_win(min_pos_anti), 'b', 'LineWidth',2);
        set(gca,'ylim',[0 0.250],'xlim', [1 length(min_pos_pro)], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti sup ' recArea]);
        axis square;
        
         %scatter
        figure; hold on;
        plot(t_win(min_pos_pro),t_win(min_pos_anti), '.k','MarkerSize', 30); hold on; plot([0 0.25],[0 0.25])
        set(gca, 'xlim', [0 0.25],'xTick', [0 0.25], 'ylim', [0 0.25], 'yTick', [0 0.25], 'TickDir', 'out', 'FontSize', 22, 'xlabel', []); box off;
        xlabel('pro'); ylabel('anti'); axis square
        
        case 'sorted_colormap_sacc_selected_both_instr'
        % get selected
        
        indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc;
        indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
        indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup;
        indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
        
        shared_exc = intersect(indx_exc_pro, indx_exc_anti); shared_sup = intersect(indx_sup_pro, indx_sup_anti);
        shared_exc_pro = ismember(indx_exc_pro, shared_exc); shared_exc_anti = ismember(indx_exc_anti, shared_exc);
        shared_sup_pro = ismember(indx_sup_pro, shared_sup); shared_sup_anti = ismember(indx_sup_anti, shared_sup);
        
%         indx_exc = pop.indx_sel.(recArea).instr_back.both.exc;
%         indx_sup = pop.indx_sel.(recArea).instr_back.both.sup;
        
        win = [-0.301 0.1];
        t_win = units(1).pro.neural.instrDir.ts_pst;
        
        
        for j=1:length(indx_exc_pro), r_pro_exc(j,:) = units(indx_exc_pro(j)).pro.neural.instrDir.rate_pst(t_win>=win(1) & t_win<=win(2)); end
        for j=1:length(indx_exc_anti), r_anti_exc(j,:) = units(indx_exc_anti(j)).anti.neural.instrDir.rate_pst(t_win>=win(1) & t_win<=win(2)); end
        
        for j=1:length(indx_sup_pro), r_pro_sup(j,:) = units(indx_sup_pro(j)).pro.neural.instrDir.rate_pst(t_win>=win(1) & t_win<=win(2));end
        for j=1:length(indx_sup_anti), r_anti_sup(j,:) = units(indx_sup_anti(j)).anti.neural.instrDir.rate_pst(t_win>=win(1) & t_win<=win(2)); end


        % pro exc
        [maxRates,pos_max_pro] = max(r_pro_exc, [], 2);
        [~,indx_max_pro] = sort(pos_max_pro);
        r_pro_norm = r_pro_exc./repmat(maxRates,[1 size(r_pro_exc,2)]);
        r_pro_sorted_exc = r_pro_norm(indx_max_pro,:); [~, max_pos_pro] =  max(r_pro_sorted_exc,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_pro_sorted_exc,1),r_pro_sorted_exc, [0 1]); hold on;
        scatter(t_win(max_pos_pro),1:size(r_pro_sorted_exc,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Instr selected Pro both ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_pro),1:size(r_pro_exc,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_exc_pro)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro exc');
        
        % anti exc
        [maxRates_anti,pos_max_anti] = max(r_anti_exc, [], 2);
        [~,indx_max_anti] = sort(pos_max_anti);
        r_anti_norm = r_anti_exc./repmat(maxRates_anti,[1 size(r_anti_exc,2)]);
        r_anti_sorted_exc = r_anti_norm(indx_max_anti,:); [~, max_pos_anti] =  max(r_anti_sorted_exc,[],2);
        % sort based on prosaccade
        % r_anti_pro_sorted_exc = r_anti_norm(indx_max_pro,:);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_anti_sorted_exc,1),r_anti_sorted_exc, [0 1]); hold on;
        scatter(t_win(max_pos_anti),1:size(r_anti_sorted_exc,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Instr selected Anti both ' recArea])
        %anti
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(max_pos_anti),1:size(r_anti_exc,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_exc_anti)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti exc');
        
        % plot anti with pro sorting
%         B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
%         figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
%         imagesc(t_win,1:size(r_anti_pro_sorted_exc,1),r_anti_pro_sorted_exc, [0 1]); hold on;
%         scatter(t_win(max_pos_pro),1:size(r_anti_pro_sorted_exc,1),40,'k','filled');
%         set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
%         vline(0, '--k'); ylabel('Neuron'); box off; title(['Sacc selected Anti both Pro sorted ' recArea])
%         
        % plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(max_pos_pro), 'm', 'LineWidth',2);
        plot(t_win(max_pos_anti), 'b', 'LineWidth',2);
        set(gca,'ylim',[win(1) win(2)],'xlim', [1 length(max_pos_anti)], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti exc ' recArea]);
        axis square;
        
        % plot without sorting
        figure; hold on;
        plot(t_win(pos_max_pro), '.m', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_anti), '.b', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_pro), 'm', 'LineWidth',2, 'MarkerSize', 30);
        plot(t_win(pos_max_anti), 'b', 'LineWidth',2, 'MarkerSize', 30);
        set(gca,'ylim',[win(1) win(2)],'xlim', [1 14], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti exc ' recArea]);
        axis square;
        
        %scatter
        figure; hold on;
        plot(t_win(max_pos_pro),t_win(max_pos_anti), '.k','MarkerSize', 30); hold on; plot([0 0.25],[0 0.25])
        set(gca, 'xlim', [win(1) win(2)],'xTick', [win(1) win(2)], 'ylim', [win(1) win(2)], 'yTick', [win(1) win(2)], 'TickDir', 'out', 'FontSize', 22, 'xlabel', []); box off;
        xlabel('pro'); ylabel('anti'); axis square

        % compute difference between both - zero no change
%         figure; hold on;
%         plot(t_win(max_pos_pro)-t_win(max_pos_anti), 'k','LineWidth',2); %% mean(t(max_pos_pro)-t(max_pos_anti)) % for difference between both
%         set(gca,'ylim',[0 0.250],'xlim', [0 17], 'TickDir', 'out', 'FontSize', 18); box off;
%         axis square; hline(0); vline(68);
        
        %% pro sup
        [maxRates_sup_pro,pos_max] = max(r_pro_sup, [], 2); [~,pos_min] = min(r_pro_sup, [], 2);
        [~,indx_max] = sort(pos_max); [~,indx_min] = sort(pos_min);
        r_pro_norm = r_pro_sup./repmat(maxRates_sup_pro,[1 size(r_pro_sup,2)]);
        r_pro_sorted_sup = r_pro_norm(indx_min,:); [~, min_pos_pro] =  min(r_pro_sorted_sup,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[1 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_pro_sorted_sup,1),r_pro_sorted_sup, [0 1]); hold on;
        scatter(t_win(min_pos_pro),1:size(r_pro_sorted_sup,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Instr Sorted Pro sup ' recArea])
       
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(min_pos_pro),1:size(r_pro_sup,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_sup_pro)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Pro sup');
        
        % anti sup
        [maxRates,anti_max] = max(r_anti_sup, [], 2); [~,pos_min] = min(r_anti_sup, [], 2);
        [~,indx_max] = sort(pos_max); [~,indx_min] = sort(pos_min);
        r_anti_norm = r_anti_sup./repmat(maxRates,[1 size(r_anti_sup,2)]);
        r_anti_sorted_sup = r_anti_norm(indx_min,:); [~, min_pos_anti] =  min(r_anti_sorted_sup,[],2);
        
        % colormap sorted
        B = makeColorMap([1 1 1],[0 0 1],[0 0 0], 100);
        figure; set(gcf,'Position',[100 200 300 300]); axes('DataAspectRatio',[1 1 1]); colormap(flip(B)); %colormap(flipud(pink)); 
        imagesc(t_win,1:size(r_anti_sorted_sup,1),r_anti_sorted_sup, [0 1]); hold on;
        scatter(t_win(min_pos_anti),1:size(r_anti_sorted_sup,1),40,'k','filled');
        set(gca,'xlim',[win(1) win(2)],'yTick',[], 'YTickLabel', [], 'TickDir', 'out', 'FontSize', 18); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Instr Sorted Anti sup ' recArea])
        %pro
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(min_pos_anti),1:size(r_anti_sup,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','k','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 length(indx_sup_anti)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        title('Anti sup');
        
        % plot pro and anti sequences on top of each other
        figure; hold on;
        plot(t_win(min_pos_pro), 'm', 'LineWidth',2);
        plot(t_win(min_pos_anti), 'b', 'LineWidth',2);
        set(gca,'ylim',[win(1) win(2)],'xlim', [1 length(min_pos_pro)], 'TickDir', 'out', 'FontSize', 18); box off;
        ylabel('time (s)'); xlabel('neuron'); axis square; title(['Pro Anti sup ' recArea]);
        axis square;
        
         %scatter
        figure; hold on;
        plot(t_win(min_pos_pro),t_win(min_pos_anti), '.k','MarkerSize', 30); hold on; plot([0 0.25],[0 0.25])
        set(gca, 'xlim', [win(1) win(2)],'xTick', [win(1) win(2)], 'ylim', [win(1) win(2)], 'yTick', [win(1) win(2)], 'TickDir', 'out', 'FontSize', 22, 'xlabel', []); box off;
        xlabel('pro'); ylabel('anti'); axis square
        
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
        
        
    case 'sorted_colormap_sacc_rand' 
        for cellNum = 1:length(units)
            indx(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        win = [-0.151 0.151]; 
        indx = find(indx);
        nunits_area = 1:length(indx);
        t = units(1).pro.neural.sacc.ts_pst; % time
        t_win = t(t>win(1) & t<win(2));
        
        %% original
        r_pro = []; r_anti = [];indx_max=[];
        for j=1:length(indx)
            r_pro(j,:) = units(indx(j)).pro.neural.sacc.rate_pst(1,t>win(1) & t<win(2)); % psth
            r_anti(j,:) = units(indx(j)).anti.neural.sacc.rate_pst(1,t>win(1) & t<win(2)); % psth
        end
        % pro
        [maxRates,pos_max_pro] = max(r_pro, [], 2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_pro),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Orig ' recArea]);
        
        % comparison
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_pro = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_pro(cell)) >=-0.151 & t_win(pos_max_pro(cell)) <= -0.075;
                win_indx_pro(cell,1) = 1;
            elseif t_win(pos_max_pro(cell)) > -0.075 & t_win(pos_max_pro(cell)) <= 0;
                win_indx_pro(cell,3) = 1;
            elseif t_win(pos_max_pro(cell)) > 0 & t_win(pos_max_pro(cell))<= 0.075;
                win_indx_pro(cell,5) = 1;
            else
                win_indx_pro(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_pro(cell)) >=-0.151 & t_win(pos_max_pro(cell)) <= 0;
                comparison_indx_pro(cell,1) = 1;
            else
                comparison_indx_pro(cell,3) = 1;
            end
        end
        
        
        % anti
        [maxRates,pos_max_anti] = max(r_anti, [], 2);
        
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_anti),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Orig ' recArea]);
        
        % get flag for position
        win_indx_anti = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_anti(cell)) >=-0.151 & t_win(pos_max_anti(cell)) <= -0.075;
                win_indx_anti(cell,1) = 1;
            elseif t_win(pos_max_anti(cell)) > -0.075 & t_win(pos_max_anti(cell)) <= 0;
                win_indx_anti(cell,3) = 1;
            elseif t_win(pos_max_anti(cell)) > 0 & t_win(pos_max_anti(cell))<= 0.075;
                win_indx_anti(cell,5) = 1;
            else
                win_indx_anti(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_anti(cell)) >=-0.151 & t_win(pos_max_anti(cell)) <= 0;
                comparison_indx_anti(cell,1) = 1;
            else
                comparison_indx_anti(cell,3) = 1;
            end
        end
        
        
        %% one
        for j=1:length(indx)
            r_pro(j,:) = units(indx(j)).pro.neural.sacc.rate_pst_rand(1,t>-0.151 & t<0.151); % psth
            r_anti(j,:) = units(indx(j)).anti.neural.sacc.rate_pst_rand(1,t>-0.151 & t<0.151); % psth
        end
        
        % pro
        [maxRates,pos_max_pro_1] = max(r_pro, [], 2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_pro_1),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand1 ' recArea]);
        
        % get flag for position
        win_indx_pro_1 = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_pro_1(cell)) >=-0.151 & t_win(pos_max_pro_1(cell)) <= -0.075;
                win_indx_pro_1(cell,2) = 1;
            elseif t_win(pos_max_pro_1(cell)) > -0.075 & t_win(pos_max_pro_1(cell)) <= 0;
                win_indx_pro_1(cell,4) = 1;
            elseif t_win(pos_max_pro_1(cell)) > 0 & t_win(pos_max_pro_1(cell))<= 0.075;
                win_indx_pro_1(cell,6) = 1;
            else
                win_indx_pro_1(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_pro_1(cell)) >=-0.151 & t_win(pos_max_pro_1(cell)) <= 0;
                comparison_indx_pro_1(cell,2) = 1;
            else
                comparison_indx_pro_1(cell,4) = 1;
            end
        end
        
        % anti
        [maxRates,pos_max_anti_1] = max(r_anti, [], 2);
        
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_anti_1),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand1 ' recArea]);
        
        % get flag for position
        win_indx_anti_1 = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_anti_1(cell)) >=-0.151 & t_win(pos_max_anti_1(cell)) <= -0.075;
                win_indx_anti_1(cell,2) = 1;
            elseif t_win(pos_max_anti_1(cell)) > -0.075 & t_win(pos_max_anti_1(cell)) <= 0;
                win_indx_anti_1(cell,4) = 1;
            elseif t_win(pos_max_anti_1(cell)) > 0 & t_win(pos_max_anti_1(cell))<= 0.075;
                win_indx_anti_1(cell,6) = 1;
            else
                win_indx_anti_1(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_anti_1(cell)) >=-0.151 & t_win(pos_max_anti_1(cell)) <= 0;
                comparison_indx_anti_1(cell,2) = 1;
            else
                comparison_indx_anti_1(cell,4) = 1;
            end
        end
        
        
        %% two
        for j=1:length(indx)
            r_pro(j,:) = units(indx(j)).pro.neural.sacc.rate_pst_rand(2,t>-0.151 & t<0.151); % psth
            r_anti(j,:) = units(indx(j)).anti.neural.sacc.rate_pst_rand(2,t>-0.151 & t<0.151); % psth
        end
        % pro
        [maxRates,pos_max_pro_2] = max(r_pro, [], 2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_pro_2),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand2 ' recArea]);
        
        % get flag for position
        win_indx_pro_2 = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_pro_2(cell)) >=-0.151 & t_win(pos_max_pro_2(cell)) <= -0.075;
                win_indx_pro_2(cell,2) = 1;
            elseif t_win(pos_max_pro_2(cell)) > -0.075 & t_win(pos_max_pro_2(cell)) <= 0;
                win_indx_pro_2(cell,4) = 1;
            elseif t_win(pos_max_pro_2(cell)) > 0 & t_win(pos_max_pro_2(cell))<= 0.075;
                win_indx_pro_2(cell,6) = 1;
            else
                win_indx_pro_2(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_pro_2(cell)) >=-0.151 & t_win(pos_max_pro_2(cell)) <= 0;
                comparison_indx_pro_2(cell,2) = 1;
            else
                comparison_indx_pro_2(cell,4) = 1;
            end
        end
        
        % anti
        [maxRates,pos_max_anti_2] = max(r_anti, [], 2);
        r_anti_norm = r_anti./repmat(maxRates,[1 size(r_anti,2)]);
        r_anti_sorted = r_anti_norm(pos_max_anti_2,:); [~, max_pos_anti] =  max(r_anti_sorted,[],2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_anti_2),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand2 ' recArea]);
        
        % get flag for position
        win_indx_anti_2 = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_anti_2(cell)) >=-0.151 & t_win(pos_max_anti_2(cell)) <= -0.075;
                win_indx_anti_2(cell,2) = 1;
            elseif t_win(pos_max_anti_2(cell)) > -0.075 & t_win(pos_max_anti_2(cell)) <= 0;
                win_indx_anti_2(cell,4) = 1;
            elseif t_win(pos_max_anti_2(cell)) > 0 & t_win(pos_max_anti_2(cell))<= 0.075;
                win_indx_anti_2(cell,6) = 1;
            else
                win_indx_anti_2(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_anti_2(cell)) >=-0.151 & t_win(pos_max_anti_2(cell)) <= 0;
                comparison_indx_anti_2(cell,2) = 1;
            else
                comparison_indx_anti_2(cell,4) = 1;
            end
        end
        
        
        
        %% three
        for j=1:length(indx)
            r_pro(j,:) = units(indx(j)).pro.neural.sacc.rate_pst_rand(3,t>-0.151 & t<0.151); % psth
            r_anti(j,:) = units(indx(j)).anti.neural.sacc.rate_pst_rand(3,t>-0.151 & t<0.151); % psth
        end
        % pro
        [maxRates,pos_max_pro_3] = max(r_pro, [], 2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_pro_3),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand3 ' recArea]);
        
        % get flag for position
        win_indx_pro_3 = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_pro_3(cell)) >=-0.151 & t_win(pos_max_pro_3(cell)) <= -0.075;
                win_indx_pro_3(cell,2) = 1;
            elseif t_win(pos_max_pro_3(cell)) > -0.075 & t_win(pos_max_pro_3(cell)) <= 0;
                win_indx_pro_3(cell,4) = 1;
            elseif t_win(pos_max_pro_3(cell)) > 0 & t_win(pos_max_pro_3(cell))<= 0.075;
                win_indx_pro_3(cell,6) = 1;
            else
                win_indx_pro_3(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_pro_3(cell)) >=-0.151 & t_win(pos_max_pro_3(cell)) <= 0;
                comparison_indx_pro_3(cell,2) = 1;
            else
                comparison_indx_pro_3(cell,4) = 1;
            end
        end
        
        
        % anti
        [maxRates,pos_max_anti_3] = max(r_anti, [], 2);
        r_anti_norm = r_anti./repmat(maxRates,[1 size(r_anti,2)]);
        r_anti_sorted = r_anti_norm(pos_max_anti_3,:); [~, max_pos_anti] =  max(r_anti_sorted,[],2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_anti_3),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_area(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand3 ' recArea]);
        
        
        % get flag for position
        win_indx_anti_3 = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_anti_3(cell)) >=-0.151 & t_win(pos_max_anti_3(cell)) <= -0.075;
                win_indx_anti_3(cell,2) = 1;
            elseif t_win(pos_max_anti_3(cell)) > -0.075 & t_win(pos_max_anti_3(cell)) <= 0;
                win_indx_anti_3(cell,4) = 1;
            elseif t_win(pos_max_anti_3(cell)) > 0 & t_win(pos_max_anti_3(cell))<= 0.075;
                win_indx_anti_3(cell,6) = 1;
            else
                win_indx_anti_3(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_anti_3(cell)) >=-0.151 & t_win(pos_max_anti_3(cell)) <= 0;
                comparison_indx_anti_3(cell,2) = 1;
            else
                comparison_indx_anti_3(cell,4) = 1;
            end
        end
        
        
        %% Compute cosine similarity index between 4 windows for pro
        
        for k = 1:length(win_indx_pro(1,:))
            if k==8
                break
            else
                cos_sim_indx_pro(k) = sum(win_indx_pro(:,k).*win_indx_pro_3(:,k+1))/(sqrt(sum(win_indx_pro(:,k))*sum(win_indx_pro_3(:,k+1))));
            end
        end
        
        grand_cos_pro(1) = mean(cos_sim_indx_pro(:,[1 3 5 7]));
        
        % plot
        figure; hold on;
        plot(cos_sim_indx_pro, '.m', 'MarkerSize', 20);
        set(gca, 'xlim', [0 8], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        %% Compute cosine similarity index between 4 windows for anti
        
        for k = 1:length(win_indx_anti(1,:))
            if k==8
                break
            else
                cos_sim_indx_anti(k) = sum(win_indx_anti(:,k).*win_indx_anti_3(:,k+1))/(sqrt(sum(win_indx_anti(:,k))*sum(win_indx_anti_3(:,k+1))));
            end
        end
        
        grand_cos_anti(1) = mean(cos_sim_indx_anti(:,[1 3 5 7]));
        
        % plot
        figure; hold on;
        plot(cos_sim_indx_anti, '.b', 'MarkerSize', 20);
        set(gca, 'xlim', [0 8], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        
        
        %% Compute cosine similarity index between 2 windows
        
        for k = 1:length(comparison_indx_pro(1,:))
            if k==4
                break
            else
                cos_sim_indx(k) = sum(comparison_indx_pro(:,k).*comparison_indx_pro_1(:,k+1))/(sqrt(sum(comparison_indx_pro(:,k))*sum(comparison_indx_pro_1(:,k+1))));
            end
        end
        
        grand_cos = mean(cos_sim_indx(:,[1 3]));
        
        %% compute cosine similarity index for 100 iterations and plot.
        
        % % pro
        % original
        r_pro = []; indx_max=[];
        for j=1:length(indx)
            r_pro(j,:) = units(indx(j)).pro.neural.sacc.rate_pst_win; % psth
        end
        % pro
        [maxRates,pos_max_pro] = max(r_pro, [], 2);
        
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_pro = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_pro(cell)) >=-0.151 & t_win(pos_max_pro(cell)) <= -0.075;
                win_indx_pro(cell,1) = 1;
            elseif t_win(pos_max_pro(cell)) > -0.075 & t_win(pos_max_pro(cell)) <= 0;
                win_indx_pro(cell,3) = 1;
            elseif t_win(pos_max_pro(cell)) > 0 & t_win(pos_max_pro(cell))<= 0.075;
                win_indx_pro(cell,5) = 1;
            else
                win_indx_pro(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_pro(cell)) >=-0.151 & t_win(pos_max_pro(cell)) <= 0;
                comparison_indx_pro(cell,1) = 1;
            else
                comparison_indx_pro(cell,3) = 1;
            end
        end
        
        % rand
        for iter = 1:length(units(indx(1)).pro.neural.sacc.rate_pst_rand(:,1))
            
            for j=1:length(indx)
                r_pro_rnd(j,:) = units(indx(j)).pro.neural.sacc.rate_pst_rand(iter,t>-0.151 & t<0.151); % psth
            end
            % pro
            [maxRates_rnd,pos_max_pro_rnd] = max(r_pro_rnd, [], 2);
            
            % get flag for position
            win_indx_pro_rnd = zeros(length(indx),8);
            for cell = 1:length(indx)
                if t_win(pos_max_pro_rnd(cell)) >=-0.151 & t_win(pos_max_pro_rnd(cell)) <= -0.075;
                    win_indx_pro_rnd(cell,2) = 1;
                elseif t_win(pos_max_pro_rnd(cell)) > -0.075 & t_win(pos_max_pro_rnd(cell)) <= 0;
                    win_indx_pro_rnd(cell,4) = 1;
                elseif t_win(pos_max_pro_rnd(cell)) > 0 & t_win(pos_max_pro_rnd(cell))<= 0.075;
                    win_indx_pro_rnd(cell,6) = 1;
                else
                    win_indx_pro_rnd(cell,8) = 1;
                end
            end
            
            % two window comparison
            for cell = 1:length(indx)
                if t_win(pos_max_pro_rnd(cell)) >=-0.151 & t_win(pos_max_pro_rnd(cell)) <= 0;
                    comparison_indx_pro_rnd(cell,2) = 1;
                else
                    comparison_indx_pro_rnd(cell,4) = 1;
                end
            end
            
            % compute cosine similarity index
            for k = 1:length(win_indx_pro(1,:))
                if k==8
                    break
                else
                    cos_sim_indx_pro(iter,k) = sum(win_indx_pro(:,k).*win_indx_pro_rnd(:,k+1))/(sqrt(sum(win_indx_pro(:,k))*sum(win_indx_pro_rnd(:,k+1))));
                end
            end
            
        end
        
        grand_cos_pro_rnd(1,:) = mean(cos_sim_indx_pro(:,[1 3 5 7]));
        grand_cos_pro_rnd_std(1,:) = std(cos_sim_indx_pro(:,[1 3 5 7]));
        
        
        % compute similarity index for two windows
        for k = 1:length(comparison_indx_pro(1,:))
            if k==4
                break
            else
                cos_sim_indx_pro_two_win(k) = sum(comparison_indx_pro(:,k).*comparison_indx_pro_rnd(:,k+1))/(sqrt(sum(comparison_indx_pro(:,k))*sum(comparison_indx_pro_rnd(:,k+1))));
            end
        end
        
        grand_cos_pro_two_Win = mean(cos_sim_indx_pro_two_win(:,[1 3]));
        
        % % anti
        % original
        r_anti = []; indx_max=[];
        for j=1:length(indx)
            r_anti(j,:) = units(indx(j)).anti.neural.sacc.rate_pst_win; % psth
        end
        % pro
        [maxRates,pos_max_anti] = max(r_anti, [], 2);
        
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_anti = zeros(length(indx),8);
        for cell = 1:length(indx)
            if t_win(pos_max_anti(cell)) >=-0.151 & t_win(pos_max_anti(cell)) <= -0.075;
                win_indx_anti(cell,1) = 1;
            elseif t_win(pos_max_anti(cell)) > -0.075 & t_win(pos_max_anti(cell)) <= 0;
                win_indx_anti(cell,3) = 1;
            elseif t_win(pos_max_anti(cell)) > 0 & t_win(pos_max_anti(cell))<= 0.075;
                win_indx_anti(cell,5) = 1;
            else
                win_indx_anti(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx)
            if t_win(pos_max_anti(cell)) >=-0.151 & t_win(pos_max_anti(cell)) <= 0;
                comparison_indx_anti(cell,1) = 1;
            else
                comparison_indx_anti(cell,3) = 1;
            end
        end
        
        % rand
        for iter = 1:length(units(indx(1)).anti.neural.sacc.rate_pst_rand(:,1))
            
            for j=1:length(indx)
                r_anti_rnd(j,:) = units(indx(j)).anti.neural.sacc.rate_pst_rand(iter,t>-0.151 & t<0.151); % psth
            end
            % pro
            [maxRates_rnd,pos_max_anti_rnd] = max(r_anti_rnd, [], 2);
            
            % get flag for position
            win_indx_pro_rnd = zeros(length(indx),8);
            for cell = 1:length(indx)
                if t_win(pos_max_anti_rnd(cell)) >=-0.151 & t_win(pos_max_anti_rnd(cell)) <= -0.075;
                    win_indx_anti_rnd(cell,2) = 1;
                elseif t_win(pos_max_anti_rnd(cell)) > -0.075 & t_win(pos_max_anti_rnd(cell)) <= 0;
                    win_indx_anti_rnd(cell,4) = 1;
                elseif t_win(pos_max_anti_rnd(cell)) > 0 & t_win(pos_max_anti_rnd(cell))<= 0.075;
                    win_indx_anti_rnd(cell,6) = 1;
                else
                    win_indx_anti_rnd(cell,8) = 1;
                end
            end
            
            % two window comparison
            for cell = 1:length(indx)
                if t_win(pos_max_anti_rnd(cell)) >=-0.151 & t_win(pos_max_anti_rnd(cell)) <= 0;
                    comparison_indx_anti_rnd(cell,2) = 1;
                else
                    comparison_indx_anti_rnd(cell,4) = 1;
                end
            end
            
            
            % compute cosine similarity index
            for k = 1:length(win_indx_anti(1,:))
                if k==8
                    break
                else
                    cos_sim_indx_anti(iter,k) = sum(win_indx_anti(:,k).*win_indx_anti_rnd(:,k+1))/(sqrt(sum(win_indx_anti(:,k))*sum(win_indx_anti_rnd(:,k+1))));
                end
            end
            
        end
        
        grand_cos_anti_rnd(1,:) = nanmean(cos_sim_indx_anti(:,[1 3 5 7]));
        grand_cos_anti_rnd_std(1,:) = nanstd(cos_sim_indx_anti(:,[1 3 5 7]));
        
        % compute similarity index for two windows
        for k = 1:length(comparison_indx_anti(1,:))
            if k==4
                break
            else
                cos_sim_indx_anti_two_win(k) = sum(comparison_indx_anti(:,k).*comparison_indx_anti_rnd(:,k+1))/(sqrt(sum(comparison_indx_anti(:,k))*sum(comparison_indx_anti_rnd(:,k+1))));
            end
        end
        
        grand_cos_anti_two_Win = mean(cos_sim_indx_anti_two_win(:,[1 3]));
        
        %% plot
        
        % pro
        figure; hold on;
        errorbar(1:4,grand_cos_pro_rnd,grand_cos_pro_rnd_std, '.m', 'MarkerSize', 20, 'CapSize', 0);
        set(gca, 'xlim', [0 4], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        % anti
        figure; hold on;
        errorbar(1:4,grand_cos_anti_rnd,grand_cos_anti_rnd_std, '.b', 'MarkerSize', 20, 'CapSize', 0);
        set(gca, 'xlim', [0 4], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
    case 'sorted_colormap_sacc_rand_selected_exc'
        
        indx_pro = pop.indx_sel.(recArea).sacc.all.pro.exc; 
        indx_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
        
        win = [-0.151 0.151]; 
        nunits_pro = 1:length(indx_pro); nunits_anti = 1:length(indx_anti);
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=win(1) & t<=win(2));
        
        %% original
        r_pro = []; r_anti = [];indx_max=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        % pro
        [maxRates,pos_max_pro] = max(r_pro, [], 2);
        % plot
        subplot(2,4,1); hold on; %set(gcf,'Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]); hold on; 
        scatterhist(t_win(pos_max_pro),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 nunits_pro(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Orig ' recArea]); axis square; 
        
        % comparison
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_pro = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro(cell)) >=-0.151 & t_win(pos_max_pro(cell)) <= -0.075;
                win_indx_pro(cell,1) = 1;
            elseif t_win(pos_max_pro(cell)) > -0.075 & t_win(pos_max_pro(cell)) <= 0;
                win_indx_pro(cell,3) = 1;
            elseif t_win(pos_max_pro(cell)) > 0 & t_win(pos_max_pro(cell))<= 0.075;
                win_indx_pro(cell,5) = 1;
            else
                win_indx_pro(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro(cell)) >=-0.151 & t_win(pos_max_pro(cell)) <= 0;
                comparison_indx_pro(cell,1) = 1;
            else
                comparison_indx_pro(cell,3) = 1;
            end
        end
        
        
        % anti
        [maxRates,pos_max_anti] = max(r_anti, [], 2);
        
        % plot
        subplot(2,4,5); hold on; %set(gcf,'Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]); hold on; 
        scatterhist(t_win(pos_max_anti),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_anti(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Orig ' recArea]); axis square
        
        % get flag for position
        win_indx_anti = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti(cell)) >=-0.151 & t_win(pos_max_anti(cell)) <= -0.075;
                win_indx_anti(cell,1) = 1;
            elseif t_win(pos_max_anti(cell)) > -0.075 & t_win(pos_max_anti(cell)) <= 0;
                win_indx_anti(cell,3) = 1;
            elseif t_win(pos_max_anti(cell)) > 0 & t_win(pos_max_anti(cell))<= 0.075;
                win_indx_anti(cell,5) = 1;
            else
                win_indx_anti(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti(cell)) >=-0.151 & t_win(pos_max_anti(cell)) <= 0;
                comparison_indx_anti(cell,1) = 1;
            else
                comparison_indx_anti(cell,3) = 1;
            end
        end
        
        
        %% one
        r_pro = []; r_anti = [];indx_max=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst_rand(1,t>=win(1) & t<=win(2)); end
        for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst_rand(1,t>=win(1) & t<=win(2)); end
        % pro
        [maxRates,pos_max_pro_1] = max(r_pro, [], 2);
        % plot
        subplot(1,4,2); set(gcf,'Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]); hold on; 
        % figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_pro_1),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_pro(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand1 ' recArea]); axis square; 
        
        % get flag for position
        win_indx_pro_1 = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro_1(cell)) >=-0.151 & t_win(pos_max_pro_1(cell)) <= -0.075;
                win_indx_pro_1(cell,2) = 1;
            elseif t_win(pos_max_pro_1(cell)) > -0.075 & t_win(pos_max_pro_1(cell)) <= 0;
                win_indx_pro_1(cell,4) = 1;
            elseif t_win(pos_max_pro_1(cell)) > 0 & t_win(pos_max_pro_1(cell))<= 0.075;
                win_indx_pro_1(cell,6) = 1;
            else
                win_indx_pro_1(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro_1(cell)) >=-0.151 & t_win(pos_max_pro_1(cell)) <= 0;
                comparison_indx_pro_1(cell,2) = 1;
            else
                comparison_indx_pro_1(cell,4) = 1;
            end
        end
        
        % anti
        [maxRates,pos_max_anti_1] = max(r_anti, [], 2);
        
        % plot
        subplot(2,4,2);set(gcf,'Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]); hold on; 
        % figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_anti_1),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_anti(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand1 ' recArea]); axis square; 
        
        % get flag for position
        win_indx_anti_1 = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti_1(cell)) >=-0.151 & t_win(pos_max_anti_1(cell)) <= -0.075;
                win_indx_anti_1(cell,2) = 1;
            elseif t_win(pos_max_anti_1(cell)) > -0.075 & t_win(pos_max_anti_1(cell)) <= 0;
                win_indx_anti_1(cell,4) = 1;
            elseif t_win(pos_max_anti_1(cell)) > 0 & t_win(pos_max_anti_1(cell))<= 0.075;
                win_indx_anti_1(cell,6) = 1;
            else
                win_indx_anti_1(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti_1(cell)) >=-0.151 & t_win(pos_max_anti_1(cell)) <= 0;
                comparison_indx_anti_1(cell,2) = 1;
            else
                comparison_indx_anti_1(cell,4) = 1;
            end
        end
        
        
        %% two
        r_pro = []; r_anti = [];indx_max=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst_rand(2,t>=win(1) & t<=win(2)); end
        for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst_rand(2,t>=win(1) & t<=win(2)); end
        % pro
        [maxRates,pos_max_pro_2] = max(r_pro, [], 2);
        % plot
        subplot(1,4,3);
        set(gcf,'Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_pro_2),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_pro(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand2 ' recArea]); axis square; 
        
        % get flag for position
        win_indx_pro_2 = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro_2(cell)) >=-0.151 & t_win(pos_max_pro_2(cell)) <= -0.075;
                win_indx_pro_2(cell,2) = 1;
            elseif t_win(pos_max_pro_2(cell)) > -0.075 & t_win(pos_max_pro_2(cell)) <= 0;
                win_indx_pro_2(cell,4) = 1;
            elseif t_win(pos_max_pro_2(cell)) > 0 & t_win(pos_max_pro_2(cell))<= 0.075;
                win_indx_pro_2(cell,6) = 1;
            else
                win_indx_pro_2(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro_2(cell)) >=-0.151 & t_win(pos_max_pro_2(cell)) <= 0;
                comparison_indx_pro_2(cell,2) = 1;
            else
                comparison_indx_pro_2(cell,4) = 1;
            end
        end
        
        % anti
        [maxRates,pos_max_anti_2] = max(r_anti, [], 2);

        % plot
        subplot(2,4,3);
        set(gcf,'Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_anti_2),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_anti(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand2 ' recArea]); axis square; 
        
        % get flag for position
        win_indx_anti_2 = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti_2(cell)) >=-0.151 & t_win(pos_max_anti_2(cell)) <= -0.075;
                win_indx_anti_2(cell,2) = 1;
            elseif t_win(pos_max_anti_2(cell)) > -0.075 & t_win(pos_max_anti_2(cell)) <= 0;
                win_indx_anti_2(cell,4) = 1;
            elseif t_win(pos_max_anti_2(cell)) > 0 & t_win(pos_max_anti_2(cell))<= 0.075;
                win_indx_anti_2(cell,6) = 1;
            else
                win_indx_anti_2(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti_2(cell)) >=-0.151 & t_win(pos_max_anti_2(cell)) <= 0;
                comparison_indx_anti_2(cell,2) = 1;
            else
                comparison_indx_anti_2(cell,4) = 1;
            end
        end
        
        
        
        %% three
        r_pro = []; r_anti = [];indx_max=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst_rand(3,t>=win(1) & t<=win(2)); end
        for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst_rand(3,t>=win(1) & t<=win(2)); end

        % pro
        [maxRates,pos_max_pro_3] = max(r_pro, [], 2);
        % plot
        subplot(1,4,4);
        set(gcf,'Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_pro_3),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_pro(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand3 ' recArea]);axis square; 
        
        % get flag for position
        win_indx_pro_3 = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro_3(cell)) >=-0.151 & t_win(pos_max_pro_3(cell)) <= -0.075
                win_indx_pro_3(cell,2) = 1;
            elseif t_win(pos_max_pro_3(cell)) > -0.075 & t_win(pos_max_pro_3(cell)) <= 0
                win_indx_pro_3(cell,4) = 1;
            elseif t_win(pos_max_pro_3(cell)) > 0 & t_win(pos_max_pro_3(cell))<= 0.075
                win_indx_pro_3(cell,6) = 1;
            else
                win_indx_pro_3(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro_3(cell)) >=-0.151 & t_win(pos_max_pro_3(cell)) <= 0
                comparison_indx_pro_3(cell,2) = 1;
            else
                comparison_indx_pro_3(cell,4) = 1;
            end
        end
        
        
        % anti
        [maxRates,pos_max_anti_3] = max(r_anti, [], 2);

        % plot
        subplot(2,4,4);
        set(gcf,'Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_max_anti_3),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_anti(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand3 ' recArea]);axis square; 
        
        
        % get flag for position
        win_indx_anti_3 = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti_3(cell)) >=-0.151 & t_win(pos_max_anti_3(cell)) <= -0.075
                win_indx_anti_3(cell,2) = 1;
            elseif t_win(pos_max_anti_3(cell)) > -0.075 & t_win(pos_max_anti_3(cell)) <= 0
                win_indx_anti_3(cell,4) = 1;
            elseif t_win(pos_max_anti_3(cell)) > 0 & t_win(pos_max_anti_3(cell))<= 0.075
                win_indx_anti_3(cell,6) = 1;
            else
                win_indx_anti_3(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti_3(cell)) >=-0.151 & t_win(pos_max_anti_3(cell)) <= 0
                comparison_indx_anti_3(cell,2) = 1;
            else
                comparison_indx_anti_3(cell,4) = 1;
            end
        end
        
        
        %% Compute cosine similarity index between 4 windows for pro
        
        for k = 1:length(win_indx_pro(1,:))
            if k==8
                break
            else
                cos_sim_indx_pro(k) = sum(win_indx_pro(:,k).*win_indx_pro_3(:,k+1))/(sqrt(sum(win_indx_pro(:,k))*sum(win_indx_pro_3(:,k+1))));
            end
        end
        
        grand_cos_pro(1) = nanmean(cos_sim_indx_pro(:,[1 3 5 7]));
        
        % plot
        figure; hold on;
        plot(cos_sim_indx_pro, '.m', 'MarkerSize', 20);
        set(gca, 'xlim', [0 8], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        %% Compute cosine similarity index between 4 windows for anti
        
        for k = 1:length(win_indx_anti(1,:))
            if k==8
                break
            else
                cos_sim_indx_anti(k) = sum(win_indx_anti(:,k).*win_indx_anti_3(:,k+1))/(sqrt(sum(win_indx_anti(:,k))*sum(win_indx_anti_3(:,k+1))));
            end
        end
        
        % grand_cos_anti(1) = mean(cos_sim_indx_anti(:,[1 3 5 7]));
        grand_cos_anti(1) = mean(cos_sim_indx_anti(:,[5 7]));
        
        % plot
        figure; hold on;
        plot(cos_sim_indx_anti, '.b', 'MarkerSize', 20);
        set(gca, 'xlim', [0 8], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        
        
        %% Compute cosine similarity index between 2 windows
        
%         for k = 1:length(comparison_indx_pro(1,:))
%             if k==4
%                 break
%             else
%                 cos_sim_indx(k) = sum(comparison_indx_pro(:,k).*comparison_indx_pro_1(:,k+1))/(sqrt(sum(comparison_indx_pro(:,k))*sum(comparison_indx_pro_1(:,k+1))));
%             end
%         end
%         
%         grand_cos = mean(cos_sim_indx(:,[1 3]));
        
        %% compute cosine similarity index for 100 iterations and plot.
        
        % % pro
        % original
        r_pro = []; indx_max=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        % pro
        [maxRates,pos_max_pro] = max(r_pro, [], 2);
        
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_pro = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro(cell)) >=-0.151 & t_win(pos_max_pro(cell)) <= -0.075
                win_indx_pro(cell,1) = 1;
            elseif t_win(pos_max_pro(cell)) > -0.075 & t_win(pos_max_pro(cell)) <= 0
                win_indx_pro(cell,3) = 1;
            elseif t_win(pos_max_pro(cell)) > 0 & t_win(pos_max_pro(cell))<= 0.075
                win_indx_pro(cell,5) = 1;
            else
                win_indx_pro(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_max_pro(cell)) >=-0.151 & t_win(pos_max_pro(cell)) <= 0
                comparison_indx_pro(cell,1) = 1;
            else
                comparison_indx_pro(cell,3) = 1;
            end
        end
        
        % rand
        for iter = 1:length(units(indx_pro(1)).pro.neural.sacc.rate_pst_rand(:,1))
            
             for j=1:length(indx_pro), r_pro_rnd(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst_rand(iter,t>=win(1) & t<=win(2)); end
            % pro
            [maxRates_rnd,pos_max_pro_rnd] = max(r_pro_rnd, [], 2);
            
            % get flag for position
            win_indx_pro_rnd = zeros(length(indx_pro),8);
            for cell = 1:length(indx_pro)
                if t_win(pos_max_pro_rnd(cell)) >=-0.151 & t_win(pos_max_pro_rnd(cell)) <= -0.075
                    win_indx_pro_rnd(cell,2) = 1;
                elseif t_win(pos_max_pro_rnd(cell)) > -0.075 & t_win(pos_max_pro_rnd(cell)) <= 0
                    win_indx_pro_rnd(cell,4) = 1;
                elseif t_win(pos_max_pro_rnd(cell)) > 0 & t_win(pos_max_pro_rnd(cell))<= 0.075
                    win_indx_pro_rnd(cell,6) = 1;
                else
                    win_indx_pro_rnd(cell,8) = 1;
                end
            end
            
            % two window comparison
            for cell = 1:length(indx_pro)
                if t_win(pos_max_pro_rnd(cell)) >=win(1) & t_win(pos_max_pro_rnd(cell)) <= 0
                    comparison_indx_pro_rnd(cell,2) = 1;
                else
                    comparison_indx_pro_rnd(cell,4) = 1;
                end
            end
            
            % compute cosine similarity index
            for k = 1:length(win_indx_pro(1,:))
                if k==8
                    break
                else
                    cos_sim_indx_pro(iter,k) = sum(win_indx_pro(:,k).*win_indx_pro_rnd(:,k+1))/(sqrt(sum(win_indx_pro(:,k))*sum(win_indx_pro_rnd(:,k+1))));
                end
            end
            
        end
        
                grand_cos_pro_rnd(1,:) = nanmean(cos_sim_indx_pro(:,[1 3 5 7]));
                grand_cos_pro_rnd_std(1,:) = nanstd(cos_sim_indx_pro(:,[1 3 5 7]));
        
        % compute similarity index for two windows
        for k = 1:length(comparison_indx_pro(1,:))
            if k==4
                break
            else
                cos_sim_indx_pro_two_win(k) = sum(comparison_indx_pro(:,k).*comparison_indx_pro_rnd(:,k+1))/(sqrt(sum(comparison_indx_pro(:,k))*sum(comparison_indx_pro_rnd(:,k+1))));
            end
        end
        
        grand_cos_pro_two_Win = mean(cos_sim_indx_pro_two_win(:,[1 3]));
        
        % % anti
        % original
        r_anti = []; indx_max=[];
         for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        % pro
        [maxRates,pos_max_anti] = max(r_anti, [], 2);
        
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_anti = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti(cell)) >=-0.151 & t_win(pos_max_anti(cell)) <= -0.075
                win_indx_anti(cell,1) = 1;
            elseif t_win(pos_max_anti(cell)) > -0.075 & t_win(pos_max_anti(cell)) <= 0
                win_indx_anti(cell,3) = 1;
            elseif t_win(pos_max_anti(cell)) > 0 & t_win(pos_max_anti(cell))<= 0.075
                win_indx_anti(cell,5) = 1;
            else
                win_indx_anti(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_max_anti(cell)) >=-0.151 & t_win(pos_max_anti(cell)) <= 0
                comparison_indx_anti(cell,1) = 1;
            else
                comparison_indx_anti(cell,3) = 1;
            end
        end
        
        % rand
        for iter = 1:length(units(indx_anti(1)).anti.neural.sacc.rate_pst_rand(:,1))
            
            for j=1:length(indx_anti), r_anti_rnd(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst_rand(iter,t>=win(1) & t<=win(2)); end
            % pro
            [maxRates_rnd,pos_max_anti_rnd] = max(r_anti_rnd, [], 2);
            
            % get flag for position
            win_indx_pro_rnd = zeros(length(indx_anti),8);
            for cell = 1:length(indx_anti)
                if t_win(pos_max_anti_rnd(cell)) >=-0.151 & t_win(pos_max_anti_rnd(cell)) <= -0.075
                    win_indx_anti_rnd(cell,2) = 1;
                elseif t_win(pos_max_anti_rnd(cell)) > -0.075 & t_win(pos_max_anti_rnd(cell)) <= 0
                    win_indx_anti_rnd(cell,4) = 1;
                elseif t_win(pos_max_anti_rnd(cell)) > 0 & t_win(pos_max_anti_rnd(cell))<= 0.075
                    win_indx_anti_rnd(cell,6) = 1;
                else
                    win_indx_anti_rnd(cell,8) = 1;
                end
            end
            
            % two window comparison
            for cell = 1:length(indx_anti)
                if t_win(pos_max_anti_rnd(cell)) >=-0.151 & t_win(pos_max_anti_rnd(cell)) <= 0
                    comparison_indx_anti_rnd(cell,2) = 1;
                else
                    comparison_indx_anti_rnd(cell,4) = 1;
                end
            end
            
            
            % compute cosine similarity index
            for k = 1:length(win_indx_anti(1,:))
                if k==8
                    break
                else
                    cos_sim_indx_anti(iter,k) = sum(win_indx_anti(:,k).*win_indx_anti_rnd(:,k+1))/(sqrt(sum(win_indx_anti(:,k))*sum(win_indx_anti_rnd(:,k+1))));
                end
            end
            
        end
        
        grand_cos_anti_rnd(1,:) = nanmean(cos_sim_indx_anti(:,[1 3 5 7]));
        grand_cos_anti_rnd_std(1,:) = nanstd(cos_sim_indx_anti(:,[1 3 5 7]));
        
        % compute similarity index for two windows
        for k = 1:length(comparison_indx_anti(1,:))
            if k==4
                break
            else
                cos_sim_indx_anti_two_win(k) = sum(comparison_indx_anti(:,k).*comparison_indx_anti_rnd(:,k+1))/(sqrt(sum(comparison_indx_anti(:,k))*sum(comparison_indx_anti_rnd(:,k+1))));
            end
        end
        
        grand_cos_anti_two_Win = mean(cos_sim_indx_anti_two_win(:,[1 3]));
        
        %% plot
        
        % pro
        figure; hold on;
        errorbar(1:4,grand_cos_pro_rnd,grand_cos_pro_rnd_std, '.m', 'MarkerSize', 20, 'CapSize', 0);
        set(gca, 'xlim', [0 4], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        % anti
        figure; hold on;
        errorbar(1:4,grand_cos_anti_rnd,grand_cos_anti_rnd_std, '.b', 'MarkerSize', 20, 'CapSize', 0);
        set(gca, 'xlim', [0 4], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
         case 'sorted_colormap_sacc_rand_selected_sup'
        
        indx_pro = pop.indx_sel.(recArea).sacc.all.pro.sup; 
        indx_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
        
        win = [-0.151 0.151]; 
        nunits_pro = 1:length(indx_pro); nunits_anti = 1:length(indx_anti);
        t = units(1).pro.neural.sacc.ts_pst; t_win = t(t>=win(1) & t<=win(2));
        
        %% original
        r_pro = []; r_anti = [];indx_min=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        % pro
        [minRates,pos_min_pro] = min(r_pro, [], 2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_min_pro),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[win(1) win(2)],'ylim',[0 nunits_pro(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Orig ' recArea]);
        
        % comparison
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_pro = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro(cell)) >=-0.151 & t_win(pos_min_pro(cell)) <= -0.075
                win_indx_pro(cell,1) = 1;
            elseif t_win(pos_min_pro(cell)) > -0.075 & t_win(pos_min_pro(cell)) <= 0
                win_indx_pro(cell,3) = 1;
            elseif t_win(pos_min_pro(cell)) > 0 & t_win(pos_min_pro(cell))<= 0.075
                win_indx_pro(cell,5) = 1;
            else
                win_indx_pro(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro(cell)) >=-0.151 & t_win(pos_min_pro(cell)) <= 0
                comparison_indx_pro(cell,1) = 1;
            else
                comparison_indx_pro(cell,3) = 1;
            end
        end
        
        
        % anti
        [minRates,pos_min_anti] = min(r_anti, [], 2);
        
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_min_anti),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_anti(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Orig ' recArea]);
        
        % get flag for position
        win_indx_anti = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti(cell)) >=-0.151 & t_win(pos_min_anti(cell)) <= -0.075
                win_indx_anti(cell,1) = 1;
            elseif t_win(pos_min_anti(cell)) > -0.075 & t_win(pos_min_anti(cell)) <= 0
                win_indx_anti(cell,3) = 1;
            elseif t_win(pos_min_anti(cell)) > 0 & t_win(pos_min_anti(cell))<= 0.075
                win_indx_anti(cell,5) = 1;
            else
                win_indx_anti(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti(cell)) >=-0.151 & t_win(pos_min_anti(cell)) <= 0
                comparison_indx_anti(cell,1) = 1;
            else
                comparison_indx_anti(cell,3) = 1;
            end
        end
        
        
        %% one
        r_pro = []; r_anti = [];indx_min=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst_rand(1,t>=win(1) & t<=win(2)); end
        for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst_rand(1,t>=win(1) & t<=win(2)); end
        % pro
        [minRates,pos_min_pro_1] = min(r_pro, [], 2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_min_pro_1),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_pro(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand1 ' recArea]);
        
        % get flag for position
        win_indx_pro_1 = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro_1(cell)) >=-0.151 & t_win(pos_min_pro_1(cell)) <= -0.075
                win_indx_pro_1(cell,2) = 1;
            elseif t_win(pos_min_pro_1(cell)) > -0.075 & t_win(pos_min_pro_1(cell)) <= 0
                win_indx_pro_1(cell,4) = 1;
            elseif t_win(pos_min_pro_1(cell)) > 0 & t_win(pos_min_pro_1(cell))<= 0.075
                win_indx_pro_1(cell,6) = 1;
            else
                win_indx_pro_1(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro_1(cell)) >=-0.151 & t_win(pos_min_pro_1(cell)) <= 0;
                comparison_indx_pro_1(cell,2) = 1;
            else
                comparison_indx_pro_1(cell,4) = 1;
            end
        end
        
        % anti
        [minRates,pos_min_anti_1] = min(r_anti, [], 2);
        
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_min_anti_1),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_anti(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand1 ' recArea]);
        
        % get flag for position
        win_indx_anti_1 = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti_1(cell)) >=-0.151 & t_win(pos_min_anti_1(cell)) <= -0.075
                win_indx_anti_1(cell,2) = 1;
            elseif t_win(pos_min_anti_1(cell)) > -0.075 & t_win(pos_min_anti_1(cell)) <= 0
                win_indx_anti_1(cell,4) = 1;
            elseif t_win(pos_min_anti_1(cell)) > 0 & t_win(pos_min_anti_1(cell))<= 0.075
                win_indx_anti_1(cell,6) = 1;
            else
                win_indx_anti_1(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti_1(cell)) >=-0.151 & t_win(pos_min_anti_1(cell)) <= 0
                comparison_indx_anti_1(cell,2) = 1;
            else
                comparison_indx_anti_1(cell,4) = 1;
            end
        end
        
        
        %% two
        r_pro = []; r_anti = [];indx_max=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst_rand(2,t>=win(1) & t<=win(2)); end
        for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst_rand(2,t>=win(1) & t<=win(2)); end
        % pro
        [minRates,pos_min_pro_2] = min(r_pro, [], 2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_min_pro_2),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_pro(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand2 ' recArea]);
        
        % get flag for position
        win_indx_pro_2 = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro_2(cell)) >=-0.151 & t_win(pos_min_pro_2(cell)) <= -0.075
                win_indx_pro_2(cell,2) = 1;
            elseif t_win(pos_min_pro_2(cell)) > -0.075 & t_win(pos_min_pro_2(cell)) <= 0
                win_indx_pro_2(cell,4) = 1;
            elseif t_win(pos_min_pro_2(cell)) > 0 & t_win(pos_min_pro_2(cell))<= 0.075
                win_indx_pro_2(cell,6) = 1;
            else
                win_indx_pro_2(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro_2(cell)) >=-0.151 & t_win(pos_min_pro_2(cell)) <= 0
                comparison_indx_pro_2(cell,2) = 1;
            else
                comparison_indx_pro_2(cell,4) = 1;
            end
        end
        
        % anti
        [minRates,pos_min_anti_2] = min(r_anti, [], 2);

        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_min_anti_2),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_anti(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand2 ' recArea]);
        
        % get flag for position
        win_indx_anti_2 = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti_2(cell)) >=-0.151 & t_win(pos_min_anti_2(cell)) <= -0.075
                win_indx_anti_2(cell,2) = 1;
            elseif t_win(pos_min_anti_2(cell)) > -0.075 & t_win(pos_min_anti_2(cell)) <= 0
                win_indx_anti_2(cell,4) = 1;
            elseif t_win(pos_min_anti_2(cell)) > 0 & t_win(pos_min_anti_2(cell))<= 0.075
                win_indx_anti_2(cell,6) = 1;
            else
                win_indx_anti_2(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti_2(cell)) >=-0.151 & t_win(pos_min_anti_2(cell)) <= 0
                comparison_indx_anti_2(cell,2) = 1;
            else
                comparison_indx_anti_2(cell,4) = 1;
            end
        end
        
        
        
        %% three
        r_pro = []; r_anti = [];indx_max=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst_rand(3,t>=win(1) & t<=win(2)); end
        for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst_rand(3,t>=win(1) & t<=win(2)); end

        % pro
        [minRates,pos_min_pro_3] = min(r_pro, [], 2);
        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_min_pro_3),1:size(r_pro,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','m','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_pro(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        vline(0, '--k'); ylabel('Neuron'); box off; title(['Rand3 ' recArea]);
        
        % get flag for position
        win_indx_pro_3 = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro_3(cell)) >=-0.151 & t_win(pos_min_pro_3(cell)) <= -0.075
                win_indx_pro_3(cell,2) = 1;
            elseif t_win(pos_min_pro_3(cell)) > -0.075 & t_win(pos_min_pro_3(cell)) <= 0
                win_indx_pro_3(cell,4) = 1;
            elseif t_win(pos_min_pro_3(cell)) > 0 & t_win(pos_min_pro_3(cell))<= 0.075
                win_indx_pro_3(cell,6) = 1;
            else
                win_indx_pro_3(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro_3(cell)) >=-0.151 & t_win(pos_min_pro_3(cell)) <= 0
                comparison_indx_pro_3(cell,2) = 1;
            else
                comparison_indx_pro_3(cell,4) = 1;
            end
        end
        
        
        % anti
        [minRates,pos_min_anti_3] = min(r_anti, [], 2);

        % plot
        figure('Position', [719 545 346 420]); axes('DataAspectRatio',[1 1 1]);
        scatterhist(t_win(pos_min_anti_3),1:size(r_anti,1), 'Kernel', 'off', 'Marker', '.', 'MarkerSize',12, 'Color','b','NBins',31);
        set(gca,'xlim',[-0.15 0.15],'ylim',[0 nunits_anti(end)], 'TickDir', 'out', 'FontSize', 22, 'XGrid', 'on','YGrid', 'on', 'xlabel', []); box off;
        ylabel('Neuron'); box off; title(['Rand3 ' recArea]);
        
        
        % get flag for position
        win_indx_anti_3 = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti_3(cell)) >=-0.151 & t_win(pos_min_anti_3(cell)) <= -0.075
                win_indx_anti_3(cell,2) = 1;
            elseif t_win(pos_min_anti_3(cell)) > -0.075 & t_win(pos_min_anti_3(cell)) <= 0
                win_indx_anti_3(cell,4) = 1;
            elseif t_win(pos_min_anti_3(cell)) > 0 & t_win(pos_min_anti_3(cell))<= 0.075
                win_indx_anti_3(cell,6) = 1;
            else
                win_indx_anti_3(cell,8) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti_3(cell)) >=-0.151 & t_win(pos_min_anti_3(cell)) <= 0
                comparison_indx_anti_3(cell,2) = 1;
            else
                comparison_indx_anti_3(cell,4) = 1;
            end
        end
        
        
        %% Compute cosine similarity index between 4 windows for pro
        
        for k = 1:length(win_indx_pro(1,:))
            if k==8
                break
            else
                cos_sim_indx_pro(k) = sum(win_indx_pro(:,k).*win_indx_pro_3(:,k+1))/(sqrt(sum(win_indx_pro(:,k))*sum(win_indx_pro_3(:,k+1))));
            end
        end
        
        grand_cos_pro(1) = nanmean(cos_sim_indx_pro(:,[1 3 5 7]));
        
        % plot
        figure; hold on;
        plot(cos_sim_indx_pro, '.m', 'MarkerSize', 20);
        set(gca, 'xlim', [0 8], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        %% Compute cosine similarity index between 4 windows for anti
        
        for k = 1:length(win_indx_anti(1,:))
            if k==8
                break
            else
                cos_sim_indx_anti(k) = sum(win_indx_anti(:,k).*win_indx_anti_3(:,k+1))/(sqrt(sum(win_indx_anti(:,k))*sum(win_indx_anti_3(:,k+1))));
            end
        end
        
        % grand_cos_anti(1) = mean(cos_sim_indx_anti(:,[1 3 5 7]));
        grand_cos_anti(1) = mean(cos_sim_indx_anti(:,[5 7]));
        
        % plot
        figure; hold on;
        plot(cos_sim_indx_anti, '.b', 'MarkerSize', 20);
        set(gca, 'xlim', [0 8], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        
        
        %% Compute cosine similarity index between 2 windows
        
%         for k = 1:length(comparison_indx_pro(1,:))
%             if k==4
%                 break
%             else
%                 cos_sim_indx(k) = sum(comparison_indx_pro(:,k).*comparison_indx_pro_1(:,k+1))/(sqrt(sum(comparison_indx_pro(:,k))*sum(comparison_indx_pro_1(:,k+1))));
%             end
%         end
%         
%         grand_cos = mean(cos_sim_indx(:,[1 3]));
        
        %% compute cosine similarity index for 100 iterations and plot.
        
        % % pro
        % original
        r_pro = []; indx_min=[];
        for j=1:length(indx_pro), r_pro(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        % pro
        [minRates,pos_min_pro] = min(r_pro, [], 2);
        
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_pro = zeros(length(indx_pro),8);
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro(cell)) >=-0.151 & t_win(pos_min_pro(cell)) <= -0.075
                win_indx_pro(cell,1) = 1;
            elseif t_win(pos_min_pro(cell)) > -0.075 & t_win(pos_min_pro(cell)) <= 0
                win_indx_pro(cell,3) = 1;
            elseif t_win(pos_min_pro(cell)) > 0 & t_win(pos_min_pro(cell))<= 0.075
                win_indx_pro(cell,5) = 1;
            else
                win_indx_pro(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_pro)
            if t_win(pos_min_pro(cell)) >=-0.151 & t_win(pos_min_pro(cell)) <= 0
                comparison_indx_pro(cell,1) = 1;
            else
                comparison_indx_pro(cell,3) = 1;
            end
        end
        
        % rand
        for iter = 1:length(units(indx_pro(1)).pro.neural.sacc.rate_pst_rand(:,1))
            
             for j=1:length(indx_pro), r_pro_rnd(j,:) = units(indx_pro(j)).pro.neural.sacc.rate_pst_rand(iter,t>=win(1) & t<=win(2)); end
            % pro
            [minRates_rnd,pos_min_pro_rnd] = min(r_pro_rnd, [], 2);
            
            % get flag for position
            win_indx_pro_rnd = zeros(length(indx_pro),8);
            for cell = 1:length(indx_pro)
                if t_win(pos_min_pro_rnd(cell)) >=-0.151 & t_win(pos_min_pro_rnd(cell)) <= -0.075
                    win_indx_pro_rnd(cell,2) = 1;
                elseif t_win(pos_min_pro_rnd(cell)) > -0.075 & t_win(pos_min_pro_rnd(cell)) <= 0
                    win_indx_pro_rnd(cell,4) = 1;
                elseif t_win(pos_min_pro_rnd(cell)) > 0 & t_win(pos_min_pro_rnd(cell))<= 0.075
                    win_indx_pro_rnd(cell,6) = 1;
                else
                    win_indx_pro_rnd(cell,8) = 1;
                end
            end
            
            % two window comparison
            for cell = 1:length(indx_pro)
                if t_win(pos_min_pro_rnd(cell)) >=win(1) & t_win(pos_min_pro_rnd(cell)) <= 0
                    comparison_indx_pro_rnd(cell,2) = 1;
                else
                    comparison_indx_pro_rnd(cell,4) = 1;
                end
            end
            
            % compute cosine similarity index
            for k = 1:length(win_indx_pro(1,:))
                if k==8
                    break
                else
                    cos_sim_indx_pro(iter,k) = sum(win_indx_pro(:,k).*win_indx_pro_rnd(:,k+1))/(sqrt(sum(win_indx_pro(:,k))*sum(win_indx_pro_rnd(:,k+1))));
                end
            end
            
        end
        
                grand_cos_pro_rnd(1,:) = nanmean(cos_sim_indx_pro(:,[1 3 5 7]));
                grand_cos_pro_rnd_std(1,:) = nanstd(cos_sim_indx_pro(:,[1 3 5 7]));
        
        % compute similarity index for two windows
        for k = 1:length(comparison_indx_pro(1,:))
            if k==4
                break
            else
                cos_sim_indx_pro_two_win(k) = sum(comparison_indx_pro(:,k).*comparison_indx_pro_rnd(:,k+1))/(sqrt(sum(comparison_indx_pro(:,k))*sum(comparison_indx_pro_rnd(:,k+1))));
            end
        end
        
        grand_cos_pro_two_Win = mean(cos_sim_indx_pro_two_win(:,[1 3]));
        
        % % anti
        % original
        r_anti = []; indx_min=[];
         for j=1:length(indx_anti), r_anti(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst(t>=win(1) & t<=win(2)); end
        % pro
        [minRates,pos_min_anti] = min(r_anti, [], 2);
        
        % w1 = -150 ms to -75 ms ; w2 = -75 ms to 0 ; w3 = 0 to 75 ms ; w4 = 75 to 150 ms;
        % get flag for position
        win_indx_anti = zeros(length(indx_anti),8);
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti(cell)) >=-0.151 & t_win(pos_min_anti(cell)) <= -0.075;
                win_indx_anti(cell,1) = 1;
            elseif t_win(pos_min_anti(cell)) > -0.075 & t_win(pos_min_anti(cell)) <= 0;
                win_indx_anti(cell,3) = 1;
            elseif t_win(pos_min_anti(cell)) > 0 & t_win(pos_min_anti(cell))<= 0.075;
                win_indx_anti(cell,5) = 1;
            else
                win_indx_anti(cell,7) = 1;
            end
        end
        
        % two window comparison
        for cell = 1:length(indx_anti)
            if t_win(pos_min_anti(cell)) >=-0.151 & t_win(pos_min_anti(cell)) <= 0;
                comparison_indx_anti(cell,1) = 1;
            else
                comparison_indx_anti(cell,3) = 1;
            end
        end
        
        % rand
        for iter = 1:length(units(indx_anti(1)).anti.neural.sacc.rate_pst_rand(:,1))
            
            for j=1:length(indx_anti), r_anti_rnd(j,:) = units(indx_anti(j)).anti.neural.sacc.rate_pst_rand(iter,t>=win(1) & t<=win(2)); end
            % pro
            [minRates_rnd,pos_min_anti_rnd] = min(r_anti_rnd, [], 2);
            
            % get flag for position
            win_indx_pro_rnd = zeros(length(indx_anti),8);
            for cell = 1:length(indx_anti)
                if t_win(pos_min_anti_rnd(cell)) >=-0.151 & t_win(pos_min_anti_rnd(cell)) <= -0.075;
                    win_indx_anti_rnd(cell,2) = 1;
                elseif t_win(pos_min_anti_rnd(cell)) > -0.075 & t_win(pos_min_anti_rnd(cell)) <= 0;
                    win_indx_anti_rnd(cell,4) = 1;
                elseif t_win(pos_min_anti_rnd(cell)) > 0 & t_win(pos_min_anti_rnd(cell))<= 0.075;
                    win_indx_anti_rnd(cell,6) = 1;
                else
                    win_indx_anti_rnd(cell,8) = 1;
                end
            end
            
            % two window comparison
            for cell = 1:length(indx_anti)
                if t_win(pos_min_anti_rnd(cell)) >=-0.151 & t_win(pos_min_anti_rnd(cell)) <= 0;
                    comparison_indx_anti_rnd(cell,2) = 1;
                else
                    comparison_indx_anti_rnd(cell,4) = 1;
                end
            end
            
            
            % compute cosine similarity index
            for k = 1:length(win_indx_anti(1,:))
                if k==8
                    break
                else
                    cos_sim_indx_anti(iter,k) = sum(win_indx_anti(:,k).*win_indx_anti_rnd(:,k+1))/(sqrt(sum(win_indx_anti(:,k))*sum(win_indx_anti_rnd(:,k+1))));
                end
            end
            
        end
        
        grand_cos_anti_rnd(1,:) = nanmean(cos_sim_indx_anti(:,[1 3 5 7]));
        grand_cos_anti_rnd_std(1,:) = nanstd(cos_sim_indx_anti(:,[1 3 5 7]));
        
        % compute similarity index for two windows
        for k = 1:length(comparison_indx_anti(1,:))
            if k==4
                break
            else
                cos_sim_indx_anti_two_win(k) = sum(comparison_indx_anti(:,k).*comparison_indx_anti_rnd(:,k+1))/(sqrt(sum(comparison_indx_anti(:,k))*sum(comparison_indx_anti_rnd(:,k+1))));
            end
        end
        
        grand_cos_anti_two_Win = mean(cos_sim_indx_anti_two_win(:,[1 3]));
        
        %% plot
        
        % pro
        figure; hold on;
        errorbar(1:4,grand_cos_pro_rnd,grand_cos_pro_rnd_std, '.m', 'MarkerSize', 20, 'CapSize', 0);
        set(gca, 'xlim', [0 4], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
        % anti
        figure; hold on;
        errorbar(1:4,grand_cos_anti_rnd,grand_cos_anti_rnd_std, '.b', 'MarkerSize', 20, 'CapSize', 0);
        set(gca, 'xlim', [0 4], 'xTick', [], 'ylim',[0 1], 'yTick', [], 'TickDir', 'out','FontSize',20); box off
        
    case 'latencies'
        indx_exc_pro_omv = pop.indx_sel.vermis.sacc.all.pro.exc; indx_exc_pro_lat = pop.indx_sel.lateral.sacc.all.pro.exc; 
        indx_exc_anti_omv = pop.indx_sel.vermis.sacc.all.anti.exc; indx_exc_anti_lat = pop.indx_sel.lateral.sacc.all.anti.exc;
        indx_sup_pro_omv = pop.indx_sel.vermis.sacc.all.pro.sup; indx_sup_pro_lat = pop.indx_sel.lateral.sacc.all.pro.sup;
        indx_sup_anti_omv = pop.indx_sel.vermis.sacc.all.anti.sup; indx_sup_anti_lat = pop.indx_sel.lateral.sacc.all.anti.sup;
        
        % omv
        for j=1:length(indx_exc_pro_omv), pro_exc_omv(j,:) = units(indx_exc_pro_omv(j)).pro.neural.sacc.peak_resp_time; end
        for j=1:length(indx_exc_anti_omv), anti_exc_omv(j,:) = units(indx_exc_anti_omv(j)).anti.neural.sacc.peak_resp_time; end 
        
        for j=1:length(indx_sup_pro_omv), pro_sup_omv(j,:) = units(indx_sup_pro_omv(j)).pro.neural.sacc.peak_resp_time; end  % the psths are flipped in extractWholeNeuronR
        for j=1:length(indx_sup_anti_omv), anti_sup_omv(j,:) = units(indx_sup_anti_omv(j)).anti.neural.sacc.peak_resp_time; end  % the psths are flipped in extractWholeNeuronR
        
        % lateral
        for j=1:length(indx_exc_pro_lat), pro_exc_lat(j,:) = units(indx_exc_pro_lat(j)).pro.neural.sacc.peak_resp_time; end
        for j=1:length(indx_exc_anti_lat), anti_exc_lat(j,:) = units(indx_exc_anti_lat(j)).anti.neural.sacc.peak_resp_time; end
        
        for j=1:length(indx_sup_pro_lat), pro_sup_lat(j,:) = units(indx_sup_pro_lat(j)).pro.neural.sacc.peak_resp_time; end  % the psths are flipped in extractWholeNeuronR
        for j=1:length(indx_sup_anti_lat), anti_sup_lat(j,:) = units(indx_sup_anti_lat(j)).anti.neural.sacc.peak_resp_time; end % the psths are flipped in extractWholeNeuronR
        
        % plot exc medial vs lateral Anti         
%         figure; boxplot(pro_exc_omv); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); xlabel('Medial'); box off; axis square;
%         figure; boxplot(pro_exc_lat); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); xlabel('Lateral'); box off; axis square;
        figure; hold on;
        plot(1,pro_exc_omv,'.k','MarkerSize',12); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); box off; axis square; plot(1,mean(pro_exc_omv),'sr','MarkerSize',12);
        plot(2,pro_exc_lat,'.k','MarkerSize',12); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); box off; axis square; plot(2,mean(pro_exc_lat),'sr','MarkerSize',12);
        ranksum(pro_exc_omv, pro_exc_lat)
        
        % plot sup medial vs lateral Anti
        %         figure; boxplot(pro_sup_omv); set (gca, 'yLim',[-0.02 0.15], 'TickDir', 'out','FontSize', 22); xlabel('Medial');box off; axis square;
        %         figure; boxplot(pro_sup_lat); set (gca, 'yLim',[-0.02 0.15], 'TickDir', 'out','FontSize', 22); xlabel('Lateral'); box off; axis square;
        plot(3,pro_sup_omv,'.k','MarkerSize',12); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); box off; axis square; plot(3,mean(pro_sup_omv),'sr','MarkerSize',12);
        plot(4,pro_sup_lat,'.k','MarkerSize',12); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); box off; axis square; plot(4,mean(pro_sup_lat),'sr','MarkerSize',12);
        ranksum(pro_sup_omv, pro_sup_lat)
        
        % plot exc medial vs lateral Anti
        %         figure; boxplot(anti_exc_omv); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); xlabel('Medial');box off; axis square;
        %         figure; boxplot(anti_exc_lat); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); xlabel('Lateral');box off; axis square;
        plot(5,anti_exc_omv,'.k','MarkerSize',12); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); box off; axis square; plot(5,mean(anti_exc_omv),'sr','MarkerSize',12);
        plot(6,anti_exc_lat,'.k','MarkerSize',12); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); box off; axis square; plot(6,mean(anti_exc_lat),'sr','MarkerSize',12);
        ranksum(anti_exc_omv, anti_exc_lat)
        
        % plot sup medial vs lateral Anti
        %         figure; boxplot(anti_sup_omv); set (gca, 'yLim',[-0.02 0.15], 'TickDir', 'out','FontSize', 22); xlabel('Medial');box off; axis square;
        %         figure; boxplot(anti_sup_lat); set (gca, 'yLim',[-0.02 0.15], 'TickDir', 'out','FontSize', 22); xlabel('Lateral');box off; axis square;
        plot(7,anti_sup_omv,'.k','MarkerSize',12); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); box off; axis square; plot(7,mean(anti_sup_omv),'sr','MarkerSize',12);
        plot(8,anti_sup_lat,'.k','MarkerSize',12); set (gca, 'yLim',[0 0.16], 'TickDir', 'out','FontSize', 22); box off; axis square; plot(8,mean(anti_sup_lat),'sr','MarkerSize',12);
        ranksum(anti_sup_omv, anti_sup_lat)
        
        ranksum(anti_exc_lat,anti_sup_lat)
        ranksum([pro_sup_omv ; anti_sup_omv],[pro_sup_lat ; anti_sup_lat])
        ranksum(pro_exc_omv,pro_sup_omv) ; ranksum(anti_exc_omv,anti_sup_omv) 
        ranksum(pro_exc_lat,pro_sup_lat) ; ranksum(anti_exc_lat,anti_sup_lat) 
        
    case 'mag_sensitivity_per_area' 
        % If calculating medial and lateral separately

% Where do the files originate from...
% mag_sens_pro = pop.stats.sacc.pro.mag_sensitivity(indx_area);
% mag_sens_anti = pop.stats.sacc.anti.mag_sensitivity(indx_area);
% save('mag_sens_omv', 'mag_sens_pro', 'mag_sens_anti');

medial = load('mag_sens_omv')
lateral = load('mag_sens_lat')

medial_mag = [medial.mag_sens_pro medial.mag_sens_anti]; 
lateral_mag = [lateral.mag_sens_pro lateral.mag_sens_anti];

[h,p] = ttest2(medial.mag_sens_pro,lateral.mag_sens_pro, 'Tail', 'left')
[h,p] = ttest2(medial.mag_sens_anti,lateral.mag_sens_anti, 'Tail', 'left')
        


        
        %% Check indexing
%         fprintf(['        >>> loading ' recArea ' cells <<< \n']);
%         for cellNum = 1:length(units)
%             indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
%         end
%         indx_area = find(indx_area);
%         
%         indx_exc_pro = pop.indx_sel.(recArea).sacc.all.pro.exc; indx_exc_anti = pop.indx_sel.(recArea).sacc.all.anti.exc;
%         indx_sup_pro = pop.indx_sel.(recArea).sacc.all.pro.sup; indx_sup_anti = pop.indx_sel.(recArea).sacc.all.anti.sup;
%         indx_pro = [indx_exc_pro indx_sup_pro];
%         indx_anti = [indx_exc_anti indx_sup_anti];
%         
%         % plot
%         figure; hold on; 
%         plot(pop.stats.sacc.pro.mag_sensitivity(indx_area),pop.stats.sacc.anti.mag_sensitivity(indx_area), '.k', 'MarkerSize', 18)
%         set(gca,'FontSize', 22, 'TickDir', 'out'); xlabel('Mag sensitivity Pro (spks/sec)'); ylabel('Mag sensitivity Anti (spks/sec)')
%         xlim([0 250]); ylim([0 250]);
%         plot([0 250],[0 250],'k'); axis square;
%         % mark the ones significantly different
%         plot(pop.stats.sacc.pro.mag_sensitivity(indx_pro),pop.stats.sacc.anti.mag_sensitivity(indx_pro), '.c', 'MarkerSize', 18);
%         plot(pop.stats.sacc.pro.mag_sensitivity(indx_anti),pop.stats.sacc.anti.mag_sensitivity(indx_anti), '.c', 'MarkerSize', 18);
%         title(['Magnitude sensitivity ' recArea]); axis square;
%         
%         ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_area)', pop.stats.sacc.anti.mag_sensitivity(indx_area)')
%         ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_pro)', pop.stats.sacc.anti.mag_sensitivity(indx_pro)')
%         ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_anti)', pop.stats.sacc.anti.mag_sensitivity(indx_anti)')
%         
%         kstest2(pop.stats.sacc.pro.mag_sensitivity(indx_area)', pop.stats.sacc.anti.mag_sensitivity(indx_area)')
%         kstest2(pop.stats.sacc.pro.mag_sensitivity(indx_pro)', pop.stats.sacc.anti.mag_sensitivity(indx_pro)')
%         kstest2(pop.stats.sacc.pro.mag_sensitivity(indx_anti)', pop.stats.sacc.anti.mag_sensitivity(indx_anti)')
%         
%         mag_sens_pro = pop.stats.sacc.pro.mag_sensitivity(indx_area);
%         mag_sens_anti = pop.stats.sacc.anti.mag_sensitivity(indx_area);
%         
%         save('mag_sens_omv', 'mag_sens_pro', 'mag_sens_anti');
%         medial = load('mag_sens_omv')
%         lateral = load('mag_sens_lat')
        
    case 'mag_sensitivity'

        indx_exc_pro_omv = pop.indx_sel.vermis.sacc.all.pro.exc; indx_exc_anti_omv = pop.indx_sel.vermis.sacc.all.anti.exc;
        indx_sup_pro_omv = pop.indx_sel.vermis.sacc.all.pro.sup; indx_sup_anti_omv = pop.indx_sel.vermis.sacc.all.anti.sup;
 
        indx_pro_omv = [indx_exc_pro_omv indx_sup_pro_omv];
        indx_anti_omv = [indx_exc_anti_omv indx_sup_anti_omv];
        indx_omv = [indx_pro_omv indx_anti_omv];
        
        
        indx_exc_pro_lat = pop.indx_sel.lateral.sacc.all.pro.exc; indx_exc_anti_lat = pop.indx_sel.lateral.sacc.all.anti.exc;
        indx_sup_pro_lat = pop.indx_sel.lateral.sacc.all.pro.sup; indx_sup_anti_lat = pop.indx_sel.lateral.sacc.all.anti.sup;
 
        indx_pro_lat = [indx_exc_pro_lat indx_sup_pro_lat];
        indx_anti_lat = [indx_exc_anti_lat indx_sup_anti_lat];
        indx_lat = [indx_pro_lat indx_anti_lat]; 
        
        % plot
        figure; hold on;
        plot(pop.stats.sacc.pro.mag_sensitivity(indx_omv),pop.stats.sacc.anti.mag_sensitivity(indx_omv), '.k', 'MarkerSize', 18)
        set(gca,'FontSize', 22, 'TickDir', 'out'); xlabel('Mag sensitivity Pro (spks/sec/deg)'); ylabel('Mag sensitivity Anti (spks/sec)')
        xlim([0 250]); ylim([0 250]);
        plot([0 250],[0 250],'k'); axis square;
        % diagonal histogram 
        figure;
        [hscatter,hbar,ax] = scatterDiagHist(pop.stats.sacc.pro.mag_sensitivity(indx_omv),pop.stats.sacc.anti.mag_sensitivity(indx_omv), 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 200 400], 'yTick', [0 200 400]); %xlim([0 400]); ylim([0 400]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square; 
        
        figure; hold on;
        plot(pop.stats.sacc.pro.mag_sensitivity(indx_lat),pop.stats.sacc.anti.mag_sensitivity(indx_lat), '.k', 'MarkerSize', 18)
        set(gca,'FontSize', 22, 'TickDir', 'out'); xlabel('Mag sensitivity Pro (spks/sec/deg)'); ylabel('Mag sensitivity Anti (spks/sec)');
        xlim([0 250]); ylim([0 250]);
        plot([0 250],[0 250],'k'); axis square;
        % diagonal histogram
        figure;
        [hscatter,hbar,ax] = scatterDiagHist(pop.stats.sacc.pro.mag_sensitivity(indx_lat),pop.stats.sacc.anti.mag_sensitivity(indx_lat), 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 200 400], 'yTick', [0 200 400]); xlim([0 400]); ylim([0 400]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square; 
        
        % comparison pro vs anti per area
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_pro_omv)', pop.stats.sacc.anti.mag_sensitivity(indx_pro_omv)')
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_anti_omv)', pop.stats.sacc.anti.mag_sensitivity(indx_anti_omv)')
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_pro_lat)', pop.stats.sacc.anti.mag_sensitivity(indx_pro_lat)')
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_anti_lat)', pop.stats.sacc.anti.mag_sensitivity(indx_anti_lat)')
        
        % Fac vs Supp per area
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_omv)',pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_omv)')
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_exc_anti_omv)',pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_omv)')
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_lat)',pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_lat)')
        ranksum(pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_lat)',pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_lat)')
        
        % Fac vs Supp per area one tailed
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_omv)',pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_omv)','tail', 'right')
        ranksum(pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_omv)',pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_omv)','tail', 'right')
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_lat)',pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_lat)','tail', 'right')
        ranksum(pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_lat)',pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_lat)','tail', 'right')
        
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_omv)', pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_lat)','tail', 'left')
        ranksum(pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_omv)', pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_lat)','tail', 'left')
        ranksum(pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_omv)', pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_lat)','tail', 'left')
        ranksum(pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_omv)', pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_lat)','tail', 'left')
        
        % Medial Pro vs anti (Fac+Sup)
        ranksum([pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_omv)'],[pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_omv)'],'tail', 'left')

        % Lateral Pro vs anti (Fac+Sup)
        ranksum([pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_lat)' ; pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_lat)'],[pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_lat)' ; pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_lat)'],'tail', 'right')
        
        % Medial+Lat Fac vs Medial+Lat Supp
        %Pro
        ranksum([pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_lat)'],[pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_lat)'])
        ranksum([pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_lat)'],[pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_lat)'],'tail', 'right')
        %Anti
        ranksum([pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_lat)'],[pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_lat)'])
        ranksum([pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_lat)'],[pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_lat)'],'tail', 'right')
        
        % Fac vs Supp (Fac Pro+anti medial/lat    vs Supp  Pro+anti medial/lat)
        ranksum([pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_lat)' ; pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_lat)'],...
            [pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_lat)' ; pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_lat)']) 
        ranksum([pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_exc_pro_lat)' ; pop.stats.sacc.anti.mag_sensitivity(indx_exc_anti_lat)'],...
            [pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_omv)' ; pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_omv)' ; pop.stats.sacc.pro.mag_sensitivity(indx_sup_pro_lat)' ; pop.stats.sacc.anti.mag_sensitivity(indx_sup_anti_lat)'],'tail', 'right')

    case 'regress_amp_pv'
        indx_exc_pro_omv = pop.indx_sel.vermis.sacc.all.pro.exc; indx_exc_anti_omv = pop.indx_sel.vermis.sacc.all.anti.exc;
        indx_sup_pro_omv = pop.indx_sel.vermis.sacc.all.pro.sup; indx_sup_anti_omv = pop.indx_sel.vermis.sacc.all.anti.sup;
        
        indx_pro_omv = [indx_exc_pro_omv indx_sup_pro_omv];
        indx_anti_omv = [indx_exc_anti_omv indx_sup_anti_omv];
        indx_omv = [indx_pro_omv indx_anti_omv];
        
        indx_exc_pro_lat = pop.indx_sel.lateral.sacc.all.pro.exc; indx_exc_anti_lat = pop.indx_sel.lateral.sacc.all.anti.exc;
        indx_sup_pro_lat = pop.indx_sel.lateral.sacc.all.pro.sup; indx_sup_anti_lat = pop.indx_sel.lateral.sacc.all.anti.sup;
        
        indx_pro_lat = [indx_exc_pro_lat indx_sup_pro_lat];
        indx_anti_lat = [indx_exc_anti_lat indx_sup_anti_lat];
        indx_lat = [indx_pro_lat indx_anti_lat];
        
        pro_amp_corr = [pop.stats.sacc.regress.coeff_pro_amp_corrected]; pro_amp_corr(3) = NaN; % removed this one outlier
        anti_amp_corr = [pop.stats.sacc.regress.coeff_anti_amp_corrected]; 
        pro_pv_corr = [pop.stats.sacc.regress.coeff_pro_pv_corrected]; 
        anti_pv_corr = [pop.stats.sacc.regress.coeff_anti_pv_corrected]; anti_pv_corr(46) = NaN;
        
        % amp
        for i = 1:length(indx_omv), uncorr_pro_amp_omv(i) = units(indx_omv(i)).stats.sacc.regress.coeff_pro_amp(2); end 
        for i = 1:length(indx_omv), uncorr_anti_amp_omv(i) = units(indx_omv(i)).stats.sacc.regress.coeff_anti_amp(2); end 
        for i = 1:length(indx_lat), uncorr_pro_amp_lat(i) = units(indx_lat(i)).stats.sacc.regress.coeff_pro_amp(2); end 
        for i = 1:length(indx_lat), uncorr_anti_amp_lat(i) = units(indx_lat(i)).stats.sacc.regress.coeff_anti_amp(2); end 
        
        % peak vel 
        for i = 1:length(indx_omv), uncorr_pro_pv_omv(i) = units(indx_omv(i)).stats.sacc.regress.coeff_pro_peakVel(2); end 
        for i = 1:length(indx_omv), uncorr_anti_pv_omv(i) = units(indx_omv(i)).stats.sacc.regress.coeff_anti_peakVel(2); end 
        for i = 1:length(indx_lat), uncorr_pro_pv_lat(i) = units(indx_lat(i)).stats.sacc.regress.coeff_pro_peakVel(2); end 
        for i = 1:length(indx_lat), uncorr_anti_pv_lat(i) = units(indx_lat(i)).stats.sacc.regress.coeff_anti_peakVel(2); end 
       
        % plot omv
        figure; hold on;
        plot(pro_amp_corr(indx_omv),anti_amp_corr(indx_omv), '.k', 'MarkerSize', 18)
        set(gca,'FontSize', 22, 'TickDir', 'out'); xlabel('Pro (corrected)'); ylabel('Anti (corrected)')
        xlim([0 2]); ylim([0 2]);
        plot([0 2],[0 2],'k'); axis square;
        [p,h] = kstest2(pro_amp_corr(indx_omv)',anti_amp_corr(indx_omv)')
        % diagonal histogram amp
        figure;
        [hscatter,hbar,ax] = scatterDiagHist(pro_amp_corr(indx_omv),anti_amp_corr(indx_omv), 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 1 2], 'yTick', [0 1 2]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square;
        xlabel('Pro'); ylabel('Anti'); title('Amplitude - Medial cb corrected')
        % diagonal histogram pv
        figure;
        [hscatter,hbar,ax] = scatterDiagHist(pro_pv_corr(indx_omv),anti_pv_corr(indx_omv), 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 1 2], 'yTick', [0 1 2]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square;
        xlabel('Pro'); ylabel('Anti'); title('Peak vel - Medial cb corrected')
        [p,h] = kstest2(pro_pv_corr(indx_omv)',anti_pv_corr(indx_omv)')
        % uncorrected amp
        figure; 
        [hscatter,hbar,ax] = scatterDiagHist(uncorr_pro_amp_omv,uncorr_anti_amp_omv, 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 100 200], 'yTick', [0 100 200]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square;
        xlabel('Pro'); ylabel('Anti'); title('Amplitude - Medial cb uncorrected')
        [p,h] = kstest2(uncorr_pro_amp_omv',uncorr_anti_amp_omv')
        % uncorrected pv
        figure; 
        [hscatter,hbar,ax] = scatterDiagHist(uncorr_pro_pv_omv,uncorr_anti_pv_omv, 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 100 200], 'yTick', [0 100 200]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square;
        xlabel('Pro'); ylabel('Anti'); title('Peak vel - Medial cb uncorrected')
        
        % plot lat
        figure; hold on;
        plot(pro_amp_corr(indx_lat),anti_amp_corr(indx_lat), '.k', 'MarkerSize', 18)
        set(gca,'FontSize', 22, 'TickDir', 'out'); xlabel('Mag sensitivity Pro (spks/sec/deg)'); ylabel('Mag sensitivity Anti (spks/sec)');
        xlim([0 2]); ylim([0 2]);
        plot([0 2],[0 2],'k'); axis square;
        [p,h] = kstest2(pro_amp_corr(indx_lat)',anti_amp_corr(indx_lat)')
        % diagonal histogram
        figure;
        [hscatter,hbar,ax] = scatterDiagHist(pro_amp_corr(indx_lat),anti_amp_corr(indx_lat), 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 1 2], 'yTick', [0 1 2]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square;
        xlabel('Pro'); ylabel('Anti'); title('Amplitude - Lateral cb corrected')
        % diagonal histogram pv
        figure;
        [hscatter,hbar,ax] = scatterDiagHist(pro_pv_corr(indx_lat),anti_pv_corr(indx_lat), 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 1 2], 'yTick', [0 1 2]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square;
        xlabel('Pro'); ylabel('Anti'); title('Peak vel - Lateral cb corrected')
        [p,h] = kstest2(pro_pv_corr(indx_lat)',anti_pv_corr(indx_lat)')
        % uncorrected amp
        figure; 
        [hscatter,hbar,ax] = scatterDiagHist(uncorr_pro_amp_lat,uncorr_anti_amp_lat, 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 100 200], 'yTick', [0 100 200], 'xlim', [0 200], 'ylim', [0 200]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square;
        xlabel('Pro'); ylabel('Anti'); title('Amplitude - Lateral cb uncorrected')
        [p,h] = kstest2(uncorr_pro_amp_lat',uncorr_anti_amp_lat')
        % uncorrected pv
        figure; 
        [hscatter,hbar,ax] = scatterDiagHist(uncorr_pro_pv_lat,uncorr_anti_pv_lat, 50);
        set(gca,'FontSize', 22, 'TickDir', 'out', 'xTick', [0 100 200], 'yTick', [0 100 200]);
        hscatter.Marker = '.'; hscatter.SizeData = 400; hscatter.MarkerFaceColor = 'k'; hscatter.MarkerEdgeColor = 'k';
        hbar.FaceColor = [0 0 0]; hbar.EdgeColor = [0 0 0]; axis square;
        xlabel('Pro'); ylabel('Anti'); title('Peak vel - Lateral cb uncorrected')
        [p,h] = kstest2(uncorr_pro_pv_lat',uncorr_anti_pv_lat')
    case 'rec_location'
        
        load rec_locations.mat
        %         cnt = 1;
        %         for i = 1:length(neuronList)
        %             if ~isempty(neuronList(i).gridLocation)
        %                 loc(cnt,:) = neuronList(i).gridLocation;
        %                 cnt=cnt+1;
        %             end
        %         end
        %
        %         [Cu,~,ic] = unique(loc, 'rows');
        %         count_mtx = accumarray(ic,1);
        
        for cellNum = 1:length(units)
            indx(cellNum) = strcmp(units(cellNum).area, recArea);
        end
        indx = find(indx);
        nunits_area = 1:length(indx);
        
        for j=1:length(indx)
            %if ~isempty(units(indx(j)).coord.loc)
            %area_loc(j,:) = units(indx(j)).coord.loc; % grid coordinate location from units file
            area_loc(j,:) = loc(indx(j),:); % grid coordinate location from rec_locations file. Third column is repetition number
            % depth(j) = units(indx(j)).coord.depth; % depth of electrode
            
            if area_loc(j,1) < 15 && strcmp(recArea,'vermis')
                area_loc(j,1) = area_loc(j,1)+11; % angled insertion lateral offset
            end
            %end
        end
        % plot
        figure; hold on;
        for ii = 1:length(area_loc)
            plot(area_loc(ii,1),area_loc(ii,2), '.k', 'MarkerSize', ((area_loc(ii,3))+10)*3);
        end
        % plot(loc(:,1),loc(:,2), '.k', 'MarkerSize', loc(:,3));
        title('vermis')
        
        
    case 'x_square'
        [h,p, chi2stat,df] = prop_test([9 12], [90 72] , true) % Exc pro omv vs lat
        [h,p, chi2stat,df] = prop_test([10 10], [90 72] , true) % Sup pro omv vs lat
        [h,p, chi2stat,df] = prop_test([10 13], [90 72] , true) % Exc anti omv vs lat
        [h,p, chi2stat,df] = prop_test([6 11], [90 72] , true) % Sup anti omv lat
        
    case 'max_sorted_response'
        
        omv = load('colormap_sorted_fr_omv.mat');
        lat = load('colormap_sorted_fr_lat.mat');
        ranksum([omv.pro_exc'; omv.anti_exc'],[lat.pro_exc'; lat.anti_exc'])
        ranksum([omv.pro_sup'; omv.anti_sup'],[lat.pro_sup'; lat.anti_sup'])
        
        ranksum(omv.pro_exc',lat.pro_exc')
        ranksum(omv.anti_exc',lat.anti_exc')
        ranksum(omv.pro_sup',lat.pro_sup')
        ranksum(omv.anti_sup',lat.anti_sup')
end
