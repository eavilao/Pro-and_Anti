function trial_hist= trial_history(units,recArea)
%         n-2 n-1 n   only sequential correct trials!!
%      1   pro pro pro       1 1 1
%      2   anti anti anti    0 0 0
%
%      3   pro pro anti      1 1 0
%      4   anti anti pro     0 0 1
%
%      5   pro anti pro      1 0 1
%      6   anti pro anti     0 1 0
%
%      7   pro anti anti     1 0 0
%      8   anti pro pro      0 1 1

% output trial_hist contains a array of numbers that correspond to the
% different options shown above

pro_anti_encoder    =  [ 1 1 1;  0 0 0 ; 1 1 0; 0 0 1; 1 0 1; 0 1 0; 1 0 0; 0 1 1];
proLegend           =   [10 8 15 14 4 11 5 1];


% extract neural responses for those trials
default_prs_pro_anti; % load prs
instr=[]; sacc=[]; 

% get indx of cells of the rec area
fprintf(['        >>> loading ' recArea ' cells <<< \n']);
for cell = 1:length(units)
    n_area(cell) = strcmp(units(cell).area, recArea);
    n_area_signif(cell) = units(cell).stats.sacc.flags.proVsAnti_sacc_ks_nspk;
    exc(cell) = units(cell).pro.neural.exc;
    sup(cell) = units(cell).pro.neural.sup; 
    monk(cell) = strcmp(units(cell).monk, 'mickey');
end

indx_area = find(n_area); % all
%indx_area = find(n_area & n_area_signif); % significantly diff
%indx_area = find(n_area & n_area_signif & monk); % significantly diff and monk


for cellNum = 1:length(indx_area)
    trial_hist          = zeros(size(units(indx_area(cellNum)).trial.behav,2),1) ; % with that don't have have three consesequtive correct trials are 0
    
    % instr -- initialize in case there are missing conditions in a cell
    instr(cellNum).units_t_hist.nspk_1 = []; instr(cellNum).units_t_hist.nspk_2 = []; instr(cellNum).units_t_hist.nspk_3 = []; instr(cellNum).units_t_hist.nspk_4 = [];
    instr(cellNum).units_t_hist.nspk_5 = []; instr(cellNum).units_t_hist.nspk_6 = []; instr(cellNum).units_t_hist.nspk_7 = []; instr(cellNum).units_t_hist.nspk_8 = [];
    
    instr(cellNum).units_t_hist.nspk_1_min1 = []; instr(cellNum).units_t_hist.nspk_2_min1 = []; instr(cellNum).units_t_hist.nspk_3_min1 = []; instr(cellNum).units_t_hist.nspk_4_min1 = [];
    instr(cellNum).units_t_hist.nspk_5_min1 = []; instr(cellNum).units_t_hist.nspk_6_min1 = []; instr(cellNum).units_t_hist.nspk_7_min1 = []; instr(cellNum).units_t_hist.nspk_8_min1 = [];
    
    instr(cellNum).units_t_hist.nspk_1_min2 = []; instr(cellNum).units_t_hist.nspk_2_min2 = []; instr(cellNum).units_t_hist.nspk_3_min2 = []; instr(cellNum).units_t_hist.nspk_4_min2 = [];
    instr(cellNum).units_t_hist.nspk_5_min2 = []; instr(cellNum).units_t_hist.nspk_6_min2 = []; instr(cellNum).units_t_hist.nspk_7_min2 = []; instr(cellNum).units_t_hist.nspk_8_min2 = [];
    
    instr(cellNum).units_t_hist.rate_pst = NaN(8,301); instr(cellNum).units_t_hist.rate_pst_min1 = NaN(8,301); instr(cellNum).units_t_hist.rate_pst_min2 = NaN(8,301); 
    
    % sacc
    sacc(cellNum).units_t_hist.nspk_1 = []; sacc(cellNum).units_t_hist.nspk_2 = []; sacc(cellNum).units_t_hist.nspk_3 = []; sacc(cellNum).units_t_hist.nspk_4 = [];
    sacc(cellNum).units_t_hist.nspk_5 = []; sacc(cellNum).units_t_hist.nspk_6 = []; sacc(cellNum).units_t_hist.nspk_7 = []; sacc(cellNum).units_t_hist.nspk_8 = [];
    
    sacc(cellNum).units_t_hist.nspk_1_min1 = []; sacc(cellNum).units_t_hist.nspk_2_min1 = []; sacc(cellNum).units_t_hist.nspk_3_min1 = []; sacc(cellNum).units_t_hist.nspk_4_min1 = [];
    sacc(cellNum).units_t_hist.nspk_5_min1 = []; sacc(cellNum).units_t_hist.nspk_6_min1 = []; sacc(cellNum).units_t_hist.nspk_7_min1 = []; sacc(cellNum).units_t_hist.nspk_8_min1 = [];
    
    sacc(cellNum).units_t_hist.nspk_1_min2 = []; sacc(cellNum).units_t_hist.nspk_2_min2 = []; sacc(cellNum).units_t_hist.nspk_3_min2 = []; sacc(cellNum).units_t_hist.nspk_4_min2 = [];
    sacc(cellNum).units_t_hist.nspk_5_min2 = []; sacc(cellNum).units_t_hist.nspk_6_min2 = []; sacc(cellNum).units_t_hist.nspk_7_min2 = []; sacc(cellNum).units_t_hist.nspk_8_min2 = [];
    
    sacc(cellNum).units_t_hist.rate_pst = NaN(8,111); sacc(cellNum).units_t_hist.rate_pst_min1 = NaN(8,111); sacc(cellNum).units_t_hist.rate_pst_min2 = NaN(8,111); 
    
    % starting from trial 3 cause trial 1 and 2 don't have a history
    for i = 3:size(units(indx_area(cellNum)).trial.behav,2)
        
        
        these_trials = units(indx_area(cellNum)).trial.behav(i-2:i);
        if sum([these_trials.correctResponse]) == 6 % all 3 trials have correct responses
            
            code(:,1) = ismember(units(indx_area(cellNum)).trial.behav(i-2).conditionCode,proLegend(1,:));
            code(:,2) = ismember(units(indx_area(cellNum)).trial.behav(i-1 ).conditionCode,proLegend(1,:));
            code(:,3) = ismember(units(indx_area(cellNum)).trial.behav( i).conditionCode,proLegend(1,:));
            
            
            trial_hist(i,:) = find(ismember(pro_anti_encoder,code,'rows'));
            
            
        end
    end
    
    
    for ii = 1:length(trial_hist)
        if trial_hist(ii) ~= 0
            id = 'SS'; analyse_sacc_align=0;
            cond = trial_hist(ii); instr(cellNum).unique_cond = unique(trial_hist);
            t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
            t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
            
            %%
            if cond == 1
                % instr
                clear cond_trl
                %nspk
                cond_trl = find(trial_hist==cond); if length(cond_trl) < 9, continue; end
                
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS>prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS < prs.instruction_win(2);
                    t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
                    tspk_instr_1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS < prs.instruction_win(2);
                    tspk_instr_1_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_1_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_1_min1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS < prs.instruction_win(2);
                    tspk_instr_1_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_1_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_1_min2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [instr(cellNum).units_t_hist.rate_pst(cond,:), instr(cellNum).units_t_hist.ts] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                instr(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);
                
                % extract n-1 and n-2 trials
                [instr(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [instr(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                % sacc
                analyse_sacc_align=1;
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc < prs.instruction_win(2);
                    t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
                    tspk_sacc_1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_1_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_1_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_1_min1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_1_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_1_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_1_min2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [sacc(cellNum).units_t_hist.rate_pst(cond,:), sacc(cellNum).units_t_hist.ts] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                sacc(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);

                % extract n-1 and n-2 trials
                [sacc(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [sacc(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                analyse_sacc_align=0;
                
                %%
            elseif cond == 2
                clear cond_trl
                %nspk
                cond_trl = find(trial_hist==cond); if length(cond_trl) < 9, continue; end
                
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS>prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS < prs.instruction_win(2);
                    t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
                    tspk_instr_2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS < prs.instruction_win(2);
                    tspk_instr_2_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_2_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_2_min1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS < prs.instruction_win(2);
                    tspk_instr_2_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_2_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_2_min2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % psth
                [instr(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==2)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);
                
                % extract n-1 and n-2 trials
                [instr(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [instr(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                % sacc
                analyse_sacc_align=1;
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc < prs.instruction_win(2);
                    t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
                    tspk_sacc_2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_2_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_2_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_2_min1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_2_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_2_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_2_min2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [sacc(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                sacc(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);

                % extract n-1 and n-2 trials
                [sacc(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [sacc(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                analyse_sacc_align=0;
                
                %%
            elseif cond == 3
                clear cond_trl
                % nspk
                cond_trl = find(trial_hist==cond); if length(cond_trl) < 9, continue; end
                
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS>prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS < prs.instruction_win(2);
                    t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
                    tspk_instr_3{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_3(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_3{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS < prs.instruction_win(2);
                    tspk_instr_3_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_3_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_3_min1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS < prs.instruction_win(2);
                    tspk_instr_3_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_3_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_3_min2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % psth
                [instr(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==3)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);
               
                % extract n-1 and n-2 trials
                [instr(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [instr(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                % sacc
                analyse_sacc_align=1;
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc < prs.instruction_win(2);
                    t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
                    tspk_sacc_3{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_3(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_3{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_3_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_3_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_3_min1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_3_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_3_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_3_min2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                
                %psth
                [sacc(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                sacc(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);

                % extract n-1 and n-2 trials
                [sacc(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [sacc(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                analyse_sacc_align=0;
                
                %%
            elseif cond == 4
                clear cond_trl
                %nspk
                cond_trl = find(trial_hist==cond); if length(cond_trl) < 9, continue; end
                
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS>prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS < prs.instruction_win(2);
                    t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
                    tspk_instr_4{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_4(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_4{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS < prs.instruction_win(2);
                    tspk_instr_4_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_4_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_4_min1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS < prs.instruction_win(2);
                    tspk_instr_4_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_4_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_4_min2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [instr(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==4)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);
                
                % extract n-1 and n-2 trials
                [instr(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [instr(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                % sacc
                analyse_sacc_align=1;
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc < prs.instruction_win(2);
                    t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
                    tspk_sacc_4{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_4(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_4{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_4_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_4_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_4_min1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_4_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_4_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_4_min2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [sacc(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                sacc(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);

                % extract n-1 and n-2 trials
                [sacc(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [sacc(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                analyse_sacc_align=0;
                
                %%
            elseif cond == 5
                clear cond_trl
                % nspk
                cond_trl = find(trial_hist==cond); if length(cond_trl) < 9, continue; end
                
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS>prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS < prs.instruction_win(2);
                    t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
                    tspk_instr_5{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_5(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_5{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS < prs.instruction_win(2);
                    tspk_instr_5_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_5_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_5_min1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS < prs.instruction_win(2);
                    tspk_instr_5_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_5_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_5_min2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % psth
                [instr(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);
                
                % extract n-1 and n-2 trials
                [instr(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [instr(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                % sacc
                analyse_sacc_align=1;
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc < prs.instruction_win(2);
                    t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
                    tspk_sacc_5{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_5(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_5{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_5_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_5_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_5_min1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_5_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_5_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_5_min2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [sacc(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                sacc(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);

                % extract n-1 and n-2 trials
                [sacc(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [sacc(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                analyse_sacc_align=0;
                
                %%
            elseif cond == 6
                clear cond_trl
                % nspk
                cond_trl = find(trial_hist==cond); if length(cond_trl) < 9, continue; end
                
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS>prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS < prs.instruction_win(2);
                    t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
                    tspk_instr_6{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_6(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_6{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min1
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS < prs.instruction_win(2);
                    tspk_instr_6_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_6_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_6_min1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS < prs.instruction_win(2);
                    tspk_instr_6_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_6_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_6_min2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % psth
                [instr(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);
                
                % extract n-1 and n-2 trials
                [instr(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [instr(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                 % sacc
                analyse_sacc_align=1;
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc < prs.instruction_win(2);
                    t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
                    tspk_sacc_6{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_6(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_6{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end

                 % min1
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_6_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_6_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_6_min1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_6_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_6_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_6_min2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [sacc(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                sacc(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);

                % extract n-1 and n-2 trials
                [sacc(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [sacc(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                analyse_sacc_align=0;
                
                %%
            elseif cond == 7
                clear cond_trl
                %nspk
                cond_trl = find(trial_hist==cond); if length(cond_trl) < 9, continue; end
                
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS>prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS < prs.instruction_win(2);
                    t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
                    tspk_instr_7{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_7(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_7{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                 % min1
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS < prs.instruction_win(2);
                    tspk_instr_7_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_7_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_7_min1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS < prs.instruction_win(2);
                    tspk_instr_7_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_7_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_7_min2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [instr(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);
            
                % extract n-1 and n-2 trials
                [instr(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [instr(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                 % sacc
                analyse_sacc_align=1;
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc < prs.instruction_win(2);
                    t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
                    tspk_sacc_7{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_7(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_7{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                 % min1
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_7_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_7_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_7_min1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_7_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_7_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_7_min2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [sacc(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                sacc(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);

                % extract n-1 and n-2 trials
                [sacc(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [sacc(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                analyse_sacc_align=0;
                
                %%
            else
                clear cond_trl
                % nspk
                cond_trl = find(trial_hist==cond); if length(cond_trl) < 9, continue; end
                
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS>prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS < prs.instruction_win(2);
                    t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
                    tspk_instr_8{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_8(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_8{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                 % min1
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS < prs.instruction_win(2);
                    tspk_instr_8_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_8_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_8_min1{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_instr
                for j = 1:length(cond_trl)
                    tspk_instr= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS < prs.instruction_win(2);
                    tspk_instr_8_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS(tspk_instr);
                    [instr(cellNum).units_t_hist.nspk_8_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_instr_8_min2{j,:},timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [instr(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);
              
                % extract n-1 and n-2 trials
                [instr(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [instr(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
                instr(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(instr(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                
                 % sacc
                analyse_sacc_align=1;
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc < prs.instruction_win(2);
                    t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
                    tspk_sacc_8{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_8(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_8{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                 % min1
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_8_min1{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-1).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_8_min1(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_8_min1{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                % min2
                clear tspk_sacc
                for j = 1:length(cond_trl)
                    tspk_sacc= units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc > prs.instruction_win(1) & units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc < prs.instruction_win(2);
                    tspk_sacc_8_min2{j,:} = units(indx_area(cellNum)).trial.neural(cond_trl(j)-2).tspk_SS_align_sacc(tspk_sacc);
                    [sacc(cellNum).units_t_hist.nspk_8_min2(j,:),~,~] = Spiketimes2RateTrial(tspk_sacc_8_min2{j,:},timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                end
                
                %psth
                [sacc(cellNum).units_t_hist.rate_pst(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
                sacc(cellNum).units_t_hist.rate_pst(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst(cond,:),prs.binwidth,prs.tsmooth);

                % extract n-1 and n-2 trials
                [sacc(cellNum).units_t_hist.rate_pst_min1(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-1),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min1(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min1(cond,:),prs.binwidth,prs.tsmooth);
                
                [sacc(cellNum).units_t_hist.rate_pst_min2(cond,:), ~] = Spiketimes2Rate(units(indx_area(cellNum)).trial.neural(find(trial_hist==cond)-2),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
                sacc(cellNum).units_t_hist.rate_pst_min2(cond,:) = smooth_pst(sacc(cellNum).units_t_hist.rate_pst_min2(cond,:),prs.binwidth,prs.tsmooth);
                analyse_sacc_align=0;
                
            end
           
            
        end
    end
end


%% plot

% instr
ts=instr(3).units_t_hist.ts;
for cellNum = 1:length(instr)
    cond_1_abs(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(1,:) - instr(cellNum).units_t_hist.rate_pst_min1(1,:));
    cond_2_abs(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(2,:) - instr(cellNum).units_t_hist.rate_pst_min1(2,:));
    cond_3_abs(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(3,:) - instr(cellNum).units_t_hist.rate_pst_min1(3,:));
    cond_4_abs(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(4,:) - instr(cellNum).units_t_hist.rate_pst_min1(4,:));
    cond_5_abs(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(5,:) - instr(cellNum).units_t_hist.rate_pst_min1(5,:));
    cond_6_abs(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(6,:) - instr(cellNum).units_t_hist.rate_pst_min1(6,:));
    cond_7_abs(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(7,:) - instr(cellNum).units_t_hist.rate_pst_min1(7,:));
    cond_8_abs(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(8,:) - instr(cellNum).units_t_hist.rate_pst_min1(8,:));
end
conds = [nanmean(cond_1_abs)' nanmean(cond_2_abs)' nanmean(cond_3_abs)' nanmean(cond_4_abs)' nanmean(cond_5_abs)' nanmean(cond_6_abs)' nanmean(cond_7_abs)' nanmean(cond_8_abs)'];

% raw plot
 figure; hold on; 
 plot(ts,nanmean(cond_1_abs)', 'Color', [0 0 0.125]); plot(ts,nanmean(cond_2_abs)','Color', [0 0 0.250]); plot(ts,nanmean(cond_3_abs)','Color', [0 0 0.375]); plot(ts,nanmean(cond_4_abs)','Color', [0 0 0.5]); 
 plot(ts,nanmean(cond_5_abs)','Color', [0 0 0.625]); plot(ts,nanmean(cond_6_abs)','Color', [0 0 0.750]); plot(ts,nanmean(cond_7_abs)','Color', [0 0 0.875]); plot(ts,nanmean(cond_8_abs),'Color', [0 0 1]); 
 set(gca,'xlim', [-0 0.355], 'FontSize', 22, 'TickDir', 'out');  title('Instr min1')
 

%  figure; hold on; 
%  plot(ts,cond_1_abs', 'Color', [0 0 0.125]); plot(ts,cond_2_abs','Color', [0 0 0.250]); plot(ts,cond_3_abs','Color', [0 0 0.375]); plot(ts,cond_4_abs','Color', [0 0 0.5]); 
%  plot(ts,cond_5_abs','Color', [0 0 0.625]); plot(ts,cond_6_abs','Color', [0 0 0.750]); plot(ts,cond_7_abs','Color', [0 0 0.875]); plot(ts,cond_8_abs,'Color', [0 0 1]); 
%  set(gca,'xlim', [-0 0.355], 'FontSize', 22, 'TickDir', 'out');  title('Instr min1')

 
 
 % plot colormap
 r = [nanmean(cond_1_abs) ; nanmean(cond_2_abs) ; nanmean(cond_3_abs) ; nanmean(cond_4_abs) ; nanmean(cond_5_abs) ; nanmean(cond_6_abs) ; nanmean(cond_7_abs) ; nanmean(cond_8_abs)];
 B = goodcolormap('bwr'); ncond = 8; 
 figure; set(gcf,'Position',[100 200 300 300]);
 hold on; colormap(winter);
 imagesc(ts,1:ncond,r,[0,15]); 
 set(gca,'xlim',[0 0.350],'ylim',[1 ncond(end)],'TickDir','Out','Fontsize',18);
 xlabel('Time (s)'); ylabel('cond');
 title('Instr min1')
 
 
 % plot normalized
%  figure; hold on; 
%  plot(ts,c_1_norm); plot(ts,c_2_norm); plot(ts,c_3_norm); plot(ts,c_4_norm); 
%  plot(ts,c_5_norm); plot(ts,c_6_norm); plot(ts,c_7_norm); plot(ts,c_8_norm); 
%  set(gca,'xlim', [-0 0.350]); vline(0,'--k')
 

% instr nspk
for cellNum = 1:length(instr)
   cond_1_nspk(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_1-instr(cellNum).units_t_hist.nspk_1_min1)); 
   cond_2_nspk(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_2-instr(cellNum).units_t_hist.nspk_2_min1)); 
   cond_3_nspk(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_3-instr(cellNum).units_t_hist.nspk_3_min1)); 
   cond_4_nspk(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_4-instr(cellNum).units_t_hist.nspk_4_min1)); 
   cond_5_nspk(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_5-instr(cellNum).units_t_hist.nspk_5_min1)); 
   cond_6_nspk(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_6-instr(cellNum).units_t_hist.nspk_6_min1)); 
   cond_7_nspk(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_7-instr(cellNum).units_t_hist.nspk_7_min1)); 
   cond_8_nspk(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_8-instr(cellNum).units_t_hist.nspk_8_min1)); 
end

conds_nspk1 = [cond_1_nspk' cond_2_nspk' cond_3_nspk' cond_4_nspk' cond_5_nspk' cond_6_nspk' cond_7_nspk' cond_8_nspk'];
[p_instr_min1,~,stats] = kruskalwallis(conds_nspk1);

figure; hold on;
plot(1,cond_1_nspk, '.k','MarkerSize', 10); plot(1,nanmean(cond_1_nspk), '.g', 'MarkerSize', 30);
plot(2,cond_2_nspk, '.k','MarkerSize', 10); plot(2,nanmean(cond_2_nspk), '.g', 'MarkerSize', 30);
plot(3,cond_3_nspk, '.k','MarkerSize', 10); plot(3,nanmean(cond_3_nspk), '.g', 'MarkerSize', 30);
plot(4,cond_4_nspk, '.k','MarkerSize', 10); plot(4,nanmean(cond_4_nspk), '.g', 'MarkerSize', 30);
plot(5,cond_5_nspk, '.k','MarkerSize', 10); plot(5,nanmean(cond_5_nspk), '.g', 'MarkerSize', 30);
plot(6,cond_6_nspk, '.k','MarkerSize', 10); plot(6,nanmean(cond_6_nspk), '.g', 'MarkerSize', 30);
plot(7,cond_7_nspk, '.k','MarkerSize', 10); plot(7,nanmean(cond_7_nspk), '.g', 'MarkerSize', 30);
plot(8,cond_8_nspk, '.k','MarkerSize', 10); plot(8,nanmean(cond_8_nspk), '.g', 'MarkerSize', 30);
set(gca, 'ylim', [0 40],'xlim',[0.5 8.5], 'TickDir', 'out', 'FontSize', 22); box off
title('instr min1')

blue       = [0.2235    0.3255    0.6431];
magenta    = [0.4588    0.1725    0.3922];
pro = ['\color{magenta}' 'pro '];
anti = ['\color{blue}' 'anti '];

 xlabels = {[pro pro pro]  ;...
             [anti anti anti];...
             [pro pro anti] ;...
             [anti anti pro];...
             [pro anti pro] ;...
             [anti pro anti];...
             [pro anti anti];...
             [anti pro pro]};
 ax=gca;
 ax.XTick = [1:8]; %
 ax.XTickLabel = (xlabels(1:8)); % with these labels
 ax.XTickLabelRotation = 45;
 
 % plot for individual cells
% p=numSubplots(length(instr));
%  for cellNum = 1:length(instr)
%      
%      cond_1 = abs(instr(cellNum).units_t_hist.nspk_1-instr(cellNum).units_t_hist.nspk_1_min1);
%      cond_2 = abs(instr(cellNum).units_t_hist.nspk_2-instr(cellNum).units_t_hist.nspk_2_min1);
%      cond_3 = abs(instr(cellNum).units_t_hist.nspk_3-instr(cellNum).units_t_hist.nspk_3_min1);
%      cond_4 = abs(instr(cellNum).units_t_hist.nspk_4-instr(cellNum).units_t_hist.nspk_4_min1);
%      cond_5 = abs(instr(cellNum).units_t_hist.nspk_5-instr(cellNum).units_t_hist.nspk_5_min1);
%      cond_6 = abs(instr(cellNum).units_t_hist.nspk_6-instr(cellNum).units_t_hist.nspk_6_min1);
%      cond_7 = abs(instr(cellNum).units_t_hist.nspk_7-instr(cellNum).units_t_hist.nspk_7_min1);
%      cond_8 = abs(instr(cellNum).units_t_hist.nspk_8-instr(cellNum).units_t_hist.nspk_8_min1);
%      
%      
%      figure(1); subplot(p(1),p(2),cellNum); hold on; 
%      if ~isempty(cond_1), plot(1,cond_1, '.k','MarkerSize', 10); plot(1,nanmean(cond_1), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_2), plot(2,cond_2, '.k','MarkerSize', 10); plot(2,nanmean(cond_2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_3), plot(3,cond_3, '.k','MarkerSize', 10); plot(3,nanmean(cond_3), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_4), plot(4,cond_4, '.k','MarkerSize', 10); plot(4,nanmean(cond_4), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_5), plot(5,cond_5, '.k','MarkerSize', 10); plot(5,nanmean(cond_5), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_6), plot(6,cond_6, '.k','MarkerSize', 10); plot(6,nanmean(cond_6), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_7), plot(7,cond_7, '.k','MarkerSize', 10); plot(7,nanmean(cond_7), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_8), plot(8,cond_8, '.k','MarkerSize', 10); plot(8,nanmean(cond_8), '.g', 'MarkerSize', 30);...
%      set(gca,'ylim', [0 max([cond_1 ; cond_2 ; cond_3 ; cond_4 ; cond_5 ; cond_6 ; cond_7 ; cond_8])], 'xlim',[0.5 8.5], 'TickDir', 'out', 'FontSize', 14); box off; end
%      
%      title(['min1 cell ' num2str(cellNum)]);
%  end
%      ax=gca;
%      ax.XTick = [1:8]; %
%      ax.XTickLabel = (xlabels(1:8)); % with these labels
%      ax.XTickLabelRotation = 45;
 

% min2
for cellNum = 1:length(instr)
  
    cond_1_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(1,:) - instr(cellNum).units_t_hist.rate_pst_min2(1,:));
    cond_2_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(2,:) - instr(cellNum).units_t_hist.rate_pst_min2(2,:));
    cond_3_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(3,:) - instr(cellNum).units_t_hist.rate_pst_min2(3,:));
    cond_4_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(4,:) - instr(cellNum).units_t_hist.rate_pst_min2(4,:));
    cond_5_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(5,:) - instr(cellNum).units_t_hist.rate_pst_min2(5,:));
    cond_6_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(6,:) - instr(cellNum).units_t_hist.rate_pst_min2(6,:));
    cond_7_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(7,:) - instr(cellNum).units_t_hist.rate_pst_min2(7,:));
    cond_8_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(8,:) - instr(cellNum).units_t_hist.rate_pst_min2(8,:));
end
% raw plot
figure; hold on; 
 plot(ts,nanmean(cond_1_abs_min2)', 'Color', [0 0 0.125]); plot(ts,nanmean(cond_2_abs_min2)','Color', [0 0 0.250]); plot(ts,nanmean(cond_3_abs_min2)','Color', [0 0 0.375]); plot(ts,nanmean(cond_4_abs_min2)','Color', [0 0 0.5]); 
 plot(ts,nanmean(cond_5_abs_min2)','Color', [0 0 0.625]); plot(ts,nanmean(cond_6_abs_min2)','Color', [0 0 0.750]); plot(ts,nanmean(cond_7_abs_min2)','Color', [0 0 0.875]); plot(ts,nanmean(cond_8_abs_min2),'Color', [0 0 1]); 
 set(gca,'xlim', [-0 0.355], 'FontSize', 22, 'TickDir', 'out'); title('instr min2')
 
%  figure; hold on; 
%  plot(ts,cond_1_abs_min2', 'Color', [0 0 0.125]); plot(ts,cond_2_abs_min2','Color', [0 0 0.250]); plot(ts,cond_3_abs_min2','Color', [0 0 0.375]); plot(ts,cond_4_abs_min2','Color', [0 0 0.5]); 
%  plot(ts,cond_5_abs_min2','Color', [0 0 0.625]); plot(ts,cond_6_abs_min2','Color', [0 0 0.750]); plot(ts,cond_7_abs_min2','Color', [0 0 0.875]); plot(ts,cond_8_abs_min2,'Color', [0 0 1]); 
%  set(gca,'xlim', [-0 0.355], 'FontSize', 22, 'TickDir', 'out'); title('instr min2')
 
 % plot colormap
 r_min2 = [nanmean(cond_1_abs_min2) ; nanmean(cond_2_abs_min2) ; nanmean(cond_3_abs_min2) ; nanmean(cond_4_abs_min2) ; nanmean(cond_5_abs_min2) ; nanmean(cond_6_abs_min2) ; nanmean(cond_7_abs_min2) ; nanmean(cond_8_abs_min2)];
 B = goodcolormap('bwr'); ncond = 8; 
 figure; set(gcf,'Position',[100 200 300 300]);
 hold on; colormap(winter);
 imagesc(ts,1:ncond,r_min2,[0,15]);
 set(gca,'xlim',[0 0.350],'ylim',[1 ncond(end)],'TickDir','Out','Fontsize',18);
 xlabel('Time (s)'); ylabel('cond');
 title('Instr min2')
 
 
 % instr nspk min2
for cellNum = 1:length(instr)
   cond_1_nspk2(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_1-instr(cellNum).units_t_hist.nspk_1_min2)); 
   cond_2_nspk2(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_2-instr(cellNum).units_t_hist.nspk_2_min2)); 
   cond_3_nspk2(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_3-instr(cellNum).units_t_hist.nspk_3_min2)); 
   cond_4_nspk2(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_4-instr(cellNum).units_t_hist.nspk_4_min2)); 
   cond_5_nspk2(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_5-instr(cellNum).units_t_hist.nspk_5_min2)); 
   cond_6_nspk2(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_6-instr(cellNum).units_t_hist.nspk_6_min2)); 
   cond_7_nspk2(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_7-instr(cellNum).units_t_hist.nspk_7_min2)); 
   cond_8_nspk2(cellNum) = nanmean(abs(instr(cellNum).units_t_hist.nspk_8-instr(cellNum).units_t_hist.nspk_8_min2)); 
end

conds_nspk2 = [cond_1_nspk' cond_2_nspk' cond_3_nspk' cond_4_nspk' cond_5_nspk' cond_6_nspk' cond_7_nspk' cond_8_nspk'];
p_instr_min2 = kruskalwallis(conds_nspk2);

figure; hold on;
plot(1,cond_1_nspk2, '.k','MarkerSize', 10); plot(1,nanmean(cond_1_nspk2), '.g', 'MarkerSize', 30);
plot(2,cond_2_nspk2, '.k','MarkerSize', 10); plot(2,nanmean(cond_2_nspk2), '.g', 'MarkerSize', 30);
plot(3,cond_3_nspk2, '.k','MarkerSize', 10); plot(3,nanmean(cond_3_nspk2), '.g', 'MarkerSize', 30);
plot(4,cond_4_nspk2, '.k','MarkerSize', 10); plot(4,nanmean(cond_4_nspk2), '.g', 'MarkerSize', 30);
plot(5,cond_5_nspk2, '.k','MarkerSize', 10); plot(5,nanmean(cond_5_nspk2), '.g', 'MarkerSize', 30);
plot(6,cond_6_nspk2, '.k','MarkerSize', 10); plot(6,nanmean(cond_6_nspk2), '.g', 'MarkerSize', 30);
plot(7,cond_7_nspk2, '.k','MarkerSize', 10); plot(7,nanmean(cond_7_nspk2), '.g', 'MarkerSize', 30);
plot(8,cond_8_nspk2, '.k','MarkerSize', 10); plot(8,nanmean(cond_8_nspk2), '.g', 'MarkerSize', 30);
set(gca, 'ylim', [0 40],'xlim',[0.5 8.5], 'TickDir', 'out', 'FontSize', 22); box off
title('instr min2')

blue       = [0.2235    0.3255    0.6431];   
magenta    = [0.4588    0.1725    0.3922];
pro = ['\color{magenta}' 'pro '];
anti = ['\color{blue}' 'anti '];

 xlabels = {[pro pro pro]  ;...
             [anti anti anti];...
             [pro pro anti] ;...
             [anti anti pro];...
             [pro anti pro] ;...
             [anti pro anti];...
             [pro anti anti];...
             [anti pro pro]};
 ax=gca;
 ax.XTick = [1:8]; %
 ax.XTickLabel = (xlabels(1:8)); % with these labels
 ax.XTickLabelRotation = 45;
 
 % instr min2
%  p=numSubplots(length(instr));
%  for cellNum = 1:length(instr)
%      
%      cond_1_min2 = abs(instr(cellNum).units_t_hist.nspk_1-instr(cellNum).units_t_hist.nspk_1_min2);
%      cond_2_min2 = abs(instr(cellNum).units_t_hist.nspk_2-instr(cellNum).units_t_hist.nspk_2_min2);
%      cond_3_min2 = abs(instr(cellNum).units_t_hist.nspk_3-instr(cellNum).units_t_hist.nspk_3_min2);
%      cond_4_min2 = abs(instr(cellNum).units_t_hist.nspk_4-instr(cellNum).units_t_hist.nspk_4_min2);
%      cond_5_min2 = abs(instr(cellNum).units_t_hist.nspk_5-instr(cellNum).units_t_hist.nspk_5_min2);
%      cond_6_min2 = abs(instr(cellNum).units_t_hist.nspk_6-instr(cellNum).units_t_hist.nspk_6_min2);
%      cond_7_min2 = abs(instr(cellNum).units_t_hist.nspk_7-instr(cellNum).units_t_hist.nspk_7_min2);
%      cond_8_min2 = abs(instr(cellNum).units_t_hist.nspk_8-instr(cellNum).units_t_hist.nspk_8_min2);
%      
%      
%      figure(1); subplot(p(1),p(2),cellNum); hold on; 
%      if ~isempty(cond_1_min2), plot(1,cond_1_min2, '.k','MarkerSize', 10); plot(1,nanmean(cond_1_min2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_2_min2), plot(2,cond_2_min2, '.k','MarkerSize', 10); plot(2,nanmean(cond_2_min2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_3_min2), plot(3,cond_3_min2, '.k','MarkerSize', 10); plot(3,nanmean(cond_3_min2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_4_min2), plot(4,cond_4_min2, '.k','MarkerSize', 10); plot(4,nanmean(cond_4_min2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_5_min2), plot(5,cond_5_min2, '.k','MarkerSize', 10); plot(5,nanmean(cond_5_min2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_6_min2), plot(6,cond_6_min2, '.k','MarkerSize', 10); plot(6,nanmean(cond_6_min2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_7_min2), plot(7,cond_7_min2, '.k','MarkerSize', 10); plot(7,nanmean(cond_7_min2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_8_min2), plot(8,cond_8_min2, '.k','MarkerSize', 10); plot(8,nanmean(cond_8_min2), '.g', 'MarkerSize', 30); ...
%      set(gca, 'ylim', [0 max([cond_1_min2 ; cond_2_min2 ; cond_3_min2 ; cond_4_min2 ; cond_5_min2 ; cond_6_min2 ; cond_7_min2 ; cond_8_min2])], 'xlim',[0.5 8.5], 'TickDir', 'out', 'FontSize', 14); box off; end
%      title(['min2 cell ' num2str(cellNum)]);
%  end
%      ax=gca;
%      ax.XTick = [1:8]; %
%      ax.XTickLabel = (xlabels(1:8)); % with these labels
%      ax.XTickLabelRotation = 45;


%% sacc
% sacc(19)=[]; % leave out for now, deal with it later. 

ts=sacc(3).units_t_hist.ts;
for cellNum = 1:length(sacc)
    cond_1_abs_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(1,:) - sacc(cellNum).units_t_hist.rate_pst_min1(1,:));
    cond_2_abs_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(2,:) - sacc(cellNum).units_t_hist.rate_pst_min1(2,:));
    cond_3_abs_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(3,:) - sacc(cellNum).units_t_hist.rate_pst_min1(3,:));
    cond_4_abs_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(4,:) - sacc(cellNum).units_t_hist.rate_pst_min1(4,:));
    cond_5_abs_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(5,:) - sacc(cellNum).units_t_hist.rate_pst_min1(5,:));
    cond_6_abs_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(6,:) - sacc(cellNum).units_t_hist.rate_pst_min1(6,:));
    cond_7_abs_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(7,:) - sacc(cellNum).units_t_hist.rate_pst_min1(7,:));
    cond_8_abs_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(8,:) - sacc(cellNum).units_t_hist.rate_pst_min1(8,:));
end


% raw plot
 figure; hold on; 
 plot(ts,nanmean(cond_1_abs_sacc)', 'Color', [0 0 0.125]); plot(ts,nanmean(cond_2_abs_sacc)','Color', [0 0 0.250]); plot(ts,nanmean(cond_3_abs_sacc)','Color', [0 0 0.375]); plot(ts,nanmean(cond_4_abs_sacc)','Color', [0 0 0.5]); 
 plot(ts,nanmean(cond_5_abs_sacc)','Color', [0 0 0.625]); plot(ts,nanmean(cond_6_abs_sacc)','Color', [0 0 0.750]); plot(ts,nanmean(cond_7_abs_sacc)','Color', [0 0 0.875]); plot(ts,nanmean(cond_8_abs_sacc),'Color', [0 0 1]); 
 set(gca,'xlim', [-0.150 0.151], 'ylim', [0 15], 'FontSize', 22, 'TickDir', 'out'); vline(0,'--k'); 
 
%  figure; hold on; 
%  plot(ts,cond_1_abs_sacc', 'Color', [0 0 0.125]); plot(ts,cond_2_abs_sacc','Color', [0 0 0.250]); plot(ts,cond_3_abs_sacc','Color', [0 0 0.375]); plot(ts,cond_4_abs_sacc','Color', [0 0 0.5]); 
%  plot(ts,cond_5_abs_sacc','Color', [0 0 0.625]); plot(ts,cond_6_abs_sacc','Color', [0 0 0.750]); plot(ts,cond_7_abs_sacc','Color', [0 0 0.875]); plot(ts,cond_8_abs_sacc,'Color', [0 0 1]); 
%  set(gca,'xlim', [-0.150 0.151], 'ylim', [0 30], 'FontSize', 22, 'TickDir', 'out'); vline(0,'--k'); 
 
 % plot colormap
 r_sacc = [nanmean(cond_1_abs_sacc) ; nanmean(cond_2_abs_sacc) ; nanmean(cond_3_abs_sacc) ; nanmean(cond_4_abs_sacc) ; nanmean(cond_5_abs_sacc) ; nanmean(cond_6_abs_sacc) ; nanmean(cond_7_abs_sacc) ; nanmean(cond_8_abs_sacc)];
 figure; set(gcf,'Position',[100 200 300 300]);
 hold on; colormap(winter);
 imagesc(ts,1:ncond,r_sacc,[0,15]);
 set(gca,'xlim',[-0.150 0.151],'ylim',[1 ncond(end)],'TickDir','Out','Fontsize',18);
 xlabel('Time (s)'); ylabel('cond');
 title('Sacc min1')
 
 
% sacc nspk
for cellNum = 1:length(sacc)
   cond_1_nspk_sacc(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_1-sacc(cellNum).units_t_hist.nspk_1_min1)); 
   cond_2_nspk_sacc(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_2-sacc(cellNum).units_t_hist.nspk_2_min1)); 
   cond_3_nspk_sacc(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_3-sacc(cellNum).units_t_hist.nspk_3_min1)); 
   cond_4_nspk_sacc(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_4-sacc(cellNum).units_t_hist.nspk_4_min1)); 
   cond_5_nspk_sacc(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_5-sacc(cellNum).units_t_hist.nspk_5_min1)); 
   cond_6_nspk_sacc(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_6-sacc(cellNum).units_t_hist.nspk_6_min1)); 
   cond_7_nspk_sacc(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_7-sacc(cellNum).units_t_hist.nspk_7_min1)); 
   cond_8_nspk_sacc(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_8-sacc(cellNum).units_t_hist.nspk_8_min1)); 
end

conds_nspk_sacc = [cond_1_nspk_sacc' cond_2_nspk_sacc' cond_3_nspk_sacc' cond_4_nspk_sacc' cond_5_nspk_sacc' cond_6_nspk_sacc' cond_7_nspk_sacc' cond_8_nspk_sacc'];
[p_sacc,~,stats_sacc] = kruskalwallis(conds_nspk_sacc);
c = multcompare(stats_sacc); 
%c = multcompare(stats_sacc, 'CType','bonferroni'); 


figure; hold on;
plot(1,cond_1_nspk_sacc, '.k','MarkerSize', 10); plot(1,nanmean(cond_1_nspk_sacc), '.g', 'MarkerSize', 30);
plot(2,cond_2_nspk_sacc, '.k','MarkerSize', 10); plot(2,nanmean(cond_2_nspk_sacc), '.g', 'MarkerSize', 30);
plot(3,cond_3_nspk_sacc, '.k','MarkerSize', 10); plot(3,nanmean(cond_3_nspk_sacc), '.g', 'MarkerSize', 30);
plot(4,cond_4_nspk_sacc, '.k','MarkerSize', 10); plot(4,nanmean(cond_4_nspk_sacc), '.g', 'MarkerSize', 30);
plot(5,cond_5_nspk_sacc, '.k','MarkerSize', 10); plot(5,nanmean(cond_5_nspk_sacc), '.g', 'MarkerSize', 30);
plot(6,cond_6_nspk_sacc, '.k','MarkerSize', 10); plot(6,nanmean(cond_6_nspk_sacc), '.g', 'MarkerSize', 30);
plot(7,cond_7_nspk_sacc, '.k','MarkerSize', 10); plot(7,nanmean(cond_7_nspk_sacc), '.g', 'MarkerSize', 30);
plot(8,cond_8_nspk_sacc, '.k','MarkerSize', 10); plot(8,nanmean(cond_8_nspk_sacc), '.g', 'MarkerSize', 30);
set(gca, 'ylim', [0 25],'xlim',[0.5 8.5], 'TickDir', 'out', 'FontSize', 22); box off
title('sacc min1')

 ax=gca;
 ax.XTick = [1:8]; %
 ax.XTickLabel = (xlabels(1:8)); % with these labels
 ax.XTickLabelRotation = 45;
 
 p=numSubplots(length(sacc));
 for cellNum = 1:length(sacc)
     
     cond_1 = abs(sacc(cellNum).units_t_hist.nspk_1-sacc(cellNum).units_t_hist.nspk_1_min1);
     cond_2 = abs(sacc(cellNum).units_t_hist.nspk_2-sacc(cellNum).units_t_hist.nspk_2_min1);
     cond_3 = abs(sacc(cellNum).units_t_hist.nspk_3-sacc(cellNum).units_t_hist.nspk_3_min1);
     cond_4 = abs(sacc(cellNum).units_t_hist.nspk_4-sacc(cellNum).units_t_hist.nspk_4_min1);
     cond_5 = abs(sacc(cellNum).units_t_hist.nspk_5-sacc(cellNum).units_t_hist.nspk_5_min1);
     cond_6 = abs(sacc(cellNum).units_t_hist.nspk_6-sacc(cellNum).units_t_hist.nspk_6_min1);
     cond_7 = abs(sacc(cellNum).units_t_hist.nspk_7-sacc(cellNum).units_t_hist.nspk_7_min1);
     cond_8 = abs(sacc(cellNum).units_t_hist.nspk_8-sacc(cellNum).units_t_hist.nspk_8_min1);
     
     
%      figure(1); subplot(p(1),p(2),cellNum); hold on; 
%      if ~isempty(cond_1), plot(1,cond_1, '.k','MarkerSize', 10); plot(1,nanmean(cond_1), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_2), plot(2,cond_2, '.k','MarkerSize', 10); plot(2,nanmean(cond_2), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_3), plot(3,cond_3, '.k','MarkerSize', 10); plot(3,nanmean(cond_3), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_4), plot(4,cond_4, '.k','MarkerSize', 10); plot(4,nanmean(cond_4), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_5), plot(5,cond_5, '.k','MarkerSize', 10); plot(5,nanmean(cond_5), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_6), plot(6,cond_6, '.k','MarkerSize', 10); plot(6,nanmean(cond_6), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_7), plot(7,cond_7, '.k','MarkerSize', 10); plot(7,nanmean(cond_7), '.g', 'MarkerSize', 30); end
%      if ~isempty(cond_8), plot(8,cond_8, '.k','MarkerSize', 10); plot(8,nanmean(cond_8), '.g', 'MarkerSize', 30);...
%       set(gca, 'ylim', [0 max([cond_1 ; cond_2 ; cond_3 ; cond_4 ; cond_5 ; cond_6 ; cond_7 ; cond_8])],'xlim',[0.5 8.5], 'TickDir', 'out', 'FontSize', 14); box off; end
%     
%      title(['min1 cell ' num2str(cellNum)]);
 end
     ax=gca;
     ax.XTick = [1:8]; %
     ax.XTickLabel = (xlabels(1:8)); % with these labels
     ax.XTickLabelRotation = 45;
 

% min2
for cellNum = 1:length(sacc)
  
    cond_1_abs_min2_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(1,:) - sacc(cellNum).units_t_hist.rate_pst_min2(1,:));
    cond_2_abs_min2_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(2,:) - sacc(cellNum).units_t_hist.rate_pst_min2(2,:));
    cond_3_abs_min2_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(3,:) - sacc(cellNum).units_t_hist.rate_pst_min2(3,:));
    cond_4_abs_min2_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(4,:) - sacc(cellNum).units_t_hist.rate_pst_min2(4,:));
    cond_5_abs_min2_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(5,:) - sacc(cellNum).units_t_hist.rate_pst_min2(5,:));
    cond_6_abs_min2_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(6,:) - sacc(cellNum).units_t_hist.rate_pst_min2(6,:));
    cond_7_abs_min2_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(7,:) - sacc(cellNum).units_t_hist.rate_pst_min2(7,:));
    cond_8_abs_min2_sacc(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(8,:) - sacc(cellNum).units_t_hist.rate_pst_min2(8,:));
end
% raw plot
 figure; hold on; 
 plot(ts,nanmean(cond_1_abs_min2_sacc)', 'Color', [0 0 0.125]); plot(ts,nanmean(cond_2_abs_min2_sacc)','Color', [0 0 0.250]); plot(ts,nanmean(cond_3_abs_min2_sacc)','Color', [0 0 0.375]); plot(ts,nanmean(cond_4_abs_min2_sacc)','Color', [0 0 0.5]); 
 plot(ts,nanmean(cond_5_abs_min2_sacc)','Color', [0 0 0.625]); plot(ts,nanmean(cond_6_abs_min2_sacc)','Color', [0 0 0.750]); plot(ts,nanmean(cond_7_abs_min2_sacc)','Color', [0 0 0.875]); plot(ts,nanmean(cond_8_abs_min2_sacc),'Color', [0 0 1]); 
 set(gca,'xlim', [-0.150 0.151], 'FontSize', 22, 'TickDir', 'out'); vline(0,'--k'); title('sacc min2'); ylim([0 15]);
 
%  figure; hold on; 
%  plot(ts,cond_1_abs_min2_sacc', 'Color', [0 0 0.125]); plot(ts,cond_2_abs_min2_sacc','Color', [0 0 0.250]); plot(ts,cond_3_abs_min2_sacc','Color', [0 0 0.375]); plot(ts,cond_4_abs_min2_sacc','Color', [0 0 0.5]); 
%  plot(ts,cond_5_abs_min2_sacc','Color', [0 0 0.625]); plot(ts,cond_6_abs_min2_sacc','Color', [0 0 0.750]); plot(ts,cond_7_abs_min2_sacc','Color', [0 0 0.875]); plot(ts,cond_8_abs_min2_sacc,'Color', [0 0 1]); 
%  set(gca,'xlim', [-0.150 0.151], 'FontSize', 22, 'TickDir', 'out'); vline(0,'--k'); title('sacc min2')
 
 % plot colormap
 r_sacc2 = [nanmean(cond_1_abs_min2_sacc) ; nanmean(cond_2_abs_min2_sacc) ; nanmean(cond_3_abs_min2_sacc) ; nanmean(cond_4_abs_min2_sacc) ; nanmean(cond_5_abs_min2_sacc) ; nanmean(cond_6_abs_min2_sacc) ; nanmean(cond_7_abs_min2_sacc) ; nanmean(cond_8_abs_min2_sacc)];
 figure; set(gcf,'Position',[100 200 300 300]);
 hold on; colormap(winter);
 imagesc(ts,1:ncond,r_sacc2,[0,15]);
 set(gca,'xlim',[-0.150 0.151],'ylim',[1 ncond(end)],'TickDir','Out','Fontsize',18);
 xlabel('Time (s)'); ylabel('cond');
 title('Sacc min2')
 
 
 % sacc nspk min2
for cellNum = 1:length(sacc)
   cond_1_nspk_sacc2(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_1-sacc(cellNum).units_t_hist.nspk_1_min2)); 
   cond_2_nspk_sacc2(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_2-sacc(cellNum).units_t_hist.nspk_2_min2)); 
   cond_3_nspk_sacc2(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_3-sacc(cellNum).units_t_hist.nspk_3_min2)); 
   cond_4_nspk_sacc2(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_4-sacc(cellNum).units_t_hist.nspk_4_min2)); 
   cond_5_nspk_sacc2(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_5-sacc(cellNum).units_t_hist.nspk_5_min2)); 
   cond_6_nspk_sacc2(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_6-sacc(cellNum).units_t_hist.nspk_6_min2)); 
   cond_7_nspk_sacc2(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_7-sacc(cellNum).units_t_hist.nspk_7_min2)); 
   cond_8_nspk_sacc2(cellNum) = nanmean(abs(sacc(cellNum).units_t_hist.nspk_8-sacc(cellNum).units_t_hist.nspk_8_min2)); 
end

conds_nspk_sacc2 = [cond_1_nspk_sacc2' cond_2_nspk_sacc2' cond_3_nspk_sacc2' cond_4_nspk_sacc2' cond_5_nspk_sacc2' cond_6_nspk_sacc2' cond_7_nspk_sacc2' cond_8_nspk_sacc2'];
[p_sacc2,~,stats_sacc2] = kruskalwallis(conds_nspk_sacc2);
c = multcompare(stats_sacc2);

figure; hold on;
plot(1,cond_1_nspk_sacc2, '.k', 'MarkerSize', 10); plot(1,nanmean(cond_1_nspk_sacc2), '.g', 'MarkerSize', 30);
plot(2,cond_2_nspk_sacc2, '.k', 'MarkerSize', 10); plot(2,nanmean(cond_2_nspk_sacc2), '.g', 'MarkerSize', 30);
plot(3,cond_3_nspk_sacc2, '.k', 'MarkerSize', 10); plot(3,nanmean(cond_3_nspk_sacc2), '.g', 'MarkerSize', 30);
plot(4,cond_4_nspk_sacc2, '.k', 'MarkerSize', 10); plot(4,nanmean(cond_4_nspk_sacc2), '.g', 'MarkerSize', 30);
plot(5,cond_5_nspk_sacc2, '.k', 'MarkerSize', 10); plot(5,nanmean(cond_5_nspk_sacc2), '.g', 'MarkerSize', 30);
plot(6,cond_6_nspk_sacc2, '.k', 'MarkerSize', 10); plot(6,nanmean(cond_6_nspk_sacc2), '.g', 'MarkerSize', 30);
plot(7,cond_7_nspk_sacc2, '.k', 'MarkerSize', 10); plot(7,nanmean(cond_7_nspk_sacc2), '.g', 'MarkerSize', 30);
plot(8,cond_8_nspk_sacc2, '.k', 'MarkerSize', 10); plot(8,nanmean(cond_8_nspk_sacc2), '.g', 'MarkerSize', 30);
hline(0,'--k');
set(gca, 'ylim', [0 25],'xlim',[0.5 8.5], 'TickDir', 'out', 'FontSize', 22); box off
title('sacc min2')


 ax=gca;
 ax.XTick = [1:8]; %
 ax.XTickLabel = (xlabels(1:8)); % with these labels
 ax.XTickLabelRotation = 45;

 p=numSubplots(length(sacc));
 for cellNum = 1:length(sacc)
     
     cond_1_min2 = abs(sacc(cellNum).units_t_hist.nspk_1-sacc(cellNum).units_t_hist.nspk_1_min2);
     cond_2_min2 = abs(sacc(cellNum).units_t_hist.nspk_2-sacc(cellNum).units_t_hist.nspk_2_min2);
     cond_3_min2 = abs(sacc(cellNum).units_t_hist.nspk_3-sacc(cellNum).units_t_hist.nspk_3_min2);
     cond_4_min2 = abs(sacc(cellNum).units_t_hist.nspk_4-sacc(cellNum).units_t_hist.nspk_4_min2);
     cond_5_min2 = abs(sacc(cellNum).units_t_hist.nspk_5-sacc(cellNum).units_t_hist.nspk_5_min2);
     cond_6_min2 = abs(sacc(cellNum).units_t_hist.nspk_6-sacc(cellNum).units_t_hist.nspk_6_min2);
     cond_7_min2 = abs(sacc(cellNum).units_t_hist.nspk_7-sacc(cellNum).units_t_hist.nspk_7_min2);
     cond_8_min2 = abs(sacc(cellNum).units_t_hist.nspk_8-sacc(cellNum).units_t_hist.nspk_8_min2);
     
     
     figure(1); subplot(p(1),p(2),cellNum); hold on; 
     if ~isempty(cond_1_min2), plot(1,cond_1_min2, '.k','MarkerSize', 10); plot(1,nanmean(cond_1_min2), '.g', 'MarkerSize', 30); end
     if ~isempty(cond_2_min2), plot(2,cond_2_min2, '.k','MarkerSize', 10); plot(2,nanmean(cond_2_min2), '.g', 'MarkerSize', 30); end
     if ~isempty(cond_3_min2), plot(3,cond_3_min2, '.k','MarkerSize', 10); plot(3,nanmean(cond_3_min2), '.g', 'MarkerSize', 30); end
     if ~isempty(cond_4_min2), plot(4,cond_4_min2, '.k','MarkerSize', 10); plot(4,nanmean(cond_4_min2), '.g', 'MarkerSize', 30); end
     if ~isempty(cond_5_min2), plot(5,cond_5_min2, '.k','MarkerSize', 10); plot(5,nanmean(cond_5_min2), '.g', 'MarkerSize', 30); end
     if ~isempty(cond_6_min2), plot(6,cond_6_min2, '.k','MarkerSize', 10); plot(6,nanmean(cond_6_min2), '.g', 'MarkerSize', 30); end
     if ~isempty(cond_7_min2), plot(7,cond_7_min2, '.k','MarkerSize', 10); plot(7,nanmean(cond_7_min2), '.g', 'MarkerSize', 30); end
     if ~isempty(cond_8_min2), plot(8,cond_8_min2, '.k','MarkerSize', 10); plot(8,nanmean(cond_8_min2), '.g', 'MarkerSize', 30); ...
     set(gca, 'ylim', [0 max([cond_1_min2 ; cond_2_min2 ; cond_3_min2 ; cond_4_min2 ; cond_5_min2 ; cond_6_min2 ; cond_7_min2 ; cond_8_min2])],'xlim',[0.5 8.5], 'TickDir', 'out', 'FontSize', 14); box off; end
     title(['min2 cell ' num2str(cellNum)]);
 end
     ax=gca;
     ax.XTick = [1:8]; %
     ax.XTickLabel = (xlabels(1:8)); % with these labels
     ax.XTickLabelRotation = 45;



end

