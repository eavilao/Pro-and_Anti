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
end
%indx_area = find(n_area); % all 
indx_area = find(n_area & n_area_signif); % significantly diff

for cellNum = 1:length(indx_area)
    trial_hist          = zeros(size(units(indx_area(cellNum)).trial.behav,2),1) ; % with that don't have have three consesequtive correct trials are 0
    
    % instr -- initialize in case there are missing conditions in a cell
    instrinstr(cellNum).units_t_hist.nspk_1 = []; instr(cellNum).units_t_hist.nspk_2 = []; instr(cellNum).units_t_hist.nspk_3 = []; instr(cellNum).units_t_hist.nspk_4 = [];
    instr(cellNum).units_t_hist.nspk_5 = []; instr(cellNum).units_t_hist.nspk_6 = []; instr(cellNum).units_t_hist.nspk_7 = []; instr(cellNum).units_t_hist.nspk_8 = [];
    
    instr(cellNum).units_t_hist.nspk_1_min1 = []; instr(cellNum).units_t_hist.nspk_2_min1 = []; instr(cellNum).units_t_hist.nspk_3_min1 = []; instr(cellNum).units_t_hist.nspk_4_min1 = [];
    instr(cellNum).units_t_hist.nspk_5_min1 = []; instr(cellNum).units_t_hist.nspk_6_min1 = []; instr(cellNum).units_t_hist.nspk_7_min1 = []; instr(cellNum).units_t_hist.nspk_8_min1 = [];
    
    instr(cellNum).units_t_hist.nspk_1_min2 = []; instr(cellNum).units_t_hist.nspk_2_min2 = []; instr(cellNum).units_t_hist.nspk_3_min2 = []; instr(cellNum).units_t_hist.nspk_4_min2 = [];
    instr(cellNum).units_t_hist.nspk_5_min2 = []; instr(cellNum).units_t_hist.nspk_6_min2 = []; instr(cellNum).units_t_hist.nspk_7_min2 = []; instr(cellNum).units_t_hist.nspk_8_min2 = [];
    
    % sacc
    sacc(cellNum).units_t_hist.nspk_1 = []; sacc(cellNum).units_t_hist.nspk_2 = []; sacc(cellNum).units_t_hist.nspk_3 = []; sacc(cellNum).units_t_hist.nspk_4 = [];
    sacc(cellNum).units_t_hist.nspk_5 = []; sacc(cellNum).units_t_hist.nspk_6 = []; sacc(cellNum).units_t_hist.nspk_7 = []; sacc(cellNum).units_t_hist.nspk_8 = [];
    
    sacc(cellNum).units_t_hist.nspk_1_min1 = []; sacc(cellNum).units_t_hist.nspk_2_min1 = []; sacc(cellNum).units_t_hist.nspk_3_min1 = []; sacc(cellNum).units_t_hist.nspk_4_min1 = [];
    sacc(cellNum).units_t_hist.nspk_5_min1 = []; sacc(cellNum).units_t_hist.nspk_6_min1 = []; sacc(cellNum).units_t_hist.nspk_7_min1 = []; sacc(cellNum).units_t_hist.nspk_8_min1 = [];
    
    sacc(cellNum).units_t_hist.nspk_1_min2 = []; sacc(cellNum).units_t_hist.nspk_2_min2 = []; sacc(cellNum).units_t_hist.nspk_3_min2 = []; sacc(cellNum).units_t_hist.nspk_4_min2 = [];
    sacc(cellNum).units_t_hist.nspk_5_min2 = []; sacc(cellNum).units_t_hist.nspk_6_min2 = []; sacc(cellNum).units_t_hist.nspk_7_min2 = []; sacc(cellNum).units_t_hist.nspk_8_min2 = [];
    
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
                cond_trl = find(trial_hist==cond);
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
                cond_trl = find(trial_hist==cond);
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
                cond_trl = find(trial_hist==cond);
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
                cond_trl = find(trial_hist==cond);
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
                cond_trl = find(trial_hist==cond);
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
                cond_trl = find(trial_hist==cond);
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
                cond_trl = find(trial_hist==cond);
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
                cond_trl = find(trial_hist==cond);
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
ts=instr(1).units_t_hist.ts; 
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

% normalize by max
c_1_norm = nanmean(cond_1_abs)./repmat(max(nanmean(cond_1_abs)),[1 size(nanmean(cond_1_abs),2)]);
c_2_norm = nanmean(cond_2_abs)./repmat(max(nanmean(cond_2_abs)),[1 size(nanmean(cond_2_abs),2)]);
c_3_norm = nanmean(cond_3_abs)./repmat(max(nanmean(cond_3_abs)),[1 size(nanmean(cond_3_abs),2)]);
c_4_norm = nanmean(cond_4_abs)./repmat(max(nanmean(cond_4_abs)),[1 size(nanmean(cond_4_abs),2)]);
c_5_norm = nanmean(cond_5_abs)./repmat(max(nanmean(cond_5_abs)),[1 size(nanmean(cond_5_abs),2)]);
c_6_norm = nanmean(cond_6_abs)./repmat(max(nanmean(cond_6_abs)),[1 size(nanmean(cond_6_abs),2)]);
c_7_norm = nanmean(cond_7_abs)./repmat(max(nanmean(cond_7_abs)),[1 size(nanmean(cond_7_abs),2)]);
c_8_norm = nanmean(cond_8_abs)./repmat(max(nanmean(cond_8_abs)),[1 size(nanmean(cond_8_abs),2)]);

% raw plot
 figure; hold on; 
 plot(ts,nanmean(cond_1_abs)'); plot(ts,nanmean(cond_2_abs)'); plot(ts,nanmean(cond_3_abs)'); plot(ts,nanmean(cond_4_abs)'); 
 plot(ts,nanmean(cond_5_abs)'); plot(ts,nanmean(cond_6_abs)'); plot(ts,nanmean(cond_7_abs)'); plot(ts,nanmean(cond_8_abs)'); 
 set(gca,'xlim', [-0.150 0.150]); vline(0,'--k')
 
 % plot normalized
 figure; hold on; 
 plot(ts,c_1_norm); plot(ts,c_2_norm); plot(ts,c_3_norm); plot(ts,c_4_norm); 
 plot(ts,c_5_norm); plot(ts,c_6_norm); plot(ts,c_7_norm); plot(ts,c_8_norm); 
 set(gca,'xlim', [-0.150 0.150]); vline(0,'--k')
 

% instr nspk
for cellNum = 1:length(instr)
   cond_1_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_1) - nanmean(instr(cellNum).units_t_hist.nspk_1_min1); 
   cond_2_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_2) - nanmean(instr(cellNum).units_t_hist.nspk_2_min1); 
   cond_3_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_3) - nanmean(instr(cellNum).units_t_hist.nspk_3_min1); 
   cond_4_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_4) - nanmean(instr(cellNum).units_t_hist.nspk_4_min1); 
   cond_5_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_5) - nanmean(instr(cellNum).units_t_hist.nspk_5_min1); 
   cond_6_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_6) - nanmean(instr(cellNum).units_t_hist.nspk_6_min1); 
   cond_7_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_7) - nanmean(instr(cellNum).units_t_hist.nspk_7_min1); 
   cond_8_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_8) - nanmean(instr(cellNum).units_t_hist.nspk_8_min1); 
end

figure; hold on;
plot(1,cond_1_nspk, '.k'); plot(1,nanmean(cond_1_nspk), '.g', 'MarkerSize', 20);
plot(2,cond_2_nspk, '.k'); plot(2,nanmean(cond_2_nspk), '.g', 'MarkerSize', 20);
plot(3,cond_3_nspk, '.k'); plot(3,nanmean(cond_3_nspk), '.g', 'MarkerSize', 20);
plot(4,cond_4_nspk, '.k'); plot(4,nanmean(cond_4_nspk), '.g', 'MarkerSize', 20);
plot(5,cond_5_nspk, '.k'); plot(5,nanmean(cond_5_nspk), '.g', 'MarkerSize', 20);
plot(6,cond_6_nspk, '.k'); plot(6,nanmean(cond_6_nspk), '.g', 'MarkerSize', 20);
plot(7,cond_7_nspk, '.k'); plot(7,nanmean(cond_7_nspk), '.g', 'MarkerSize', 20);
plot(8,cond_8_nspk, '.k'); plot(8,nanmean(cond_8_nspk), '.g', 'MarkerSize', 20);
hline(0,'--k');
set(gca, 'ylim', [-20 20], 'TickDir', 'out'); box off
title('instr min1')

% min2
for cellNum = 1:length(instr)
  
    cond_1_abs_min2(cellNum,:) = abs(instr(cellNum).units_t_hist.rate_pst(1,:) - instr(cellNum).units_t_hist.rate_pst_min2(1,:));
    c_1_norm
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
 plot(ts,nanmean(cond_1_abs)'); plot(ts,nanmean(cond_2_abs)'); plot(ts,nanmean(cond_3_abs)'); plot(ts,nanmean(cond_4_abs)'); 
 plot(ts,nanmean(cond_5_abs)'); plot(ts,nanmean(cond_6_abs)'); plot(ts,nanmean(cond_7_abs)'); plot(ts,nanmean(cond_8_abs)'); 
 set(gca,'xlim', [0 0.350]); vline(0,'--k'); title('instr min2')

% sacc nspk
for cellNum = 1:length(instr)
   cond_1_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_1) - nanmean(instr(cellNum).units_t_hist.nspk_1_min2); 
   cond_2_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_2) - nanmean(instr(cellNum).units_t_hist.nspk_2_min2); 
   cond_3_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_3) - nanmean(instr(cellNum).units_t_hist.nspk_3_min2); 
   cond_4_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_4) - nanmean(instr(cellNum).units_t_hist.nspk_4_min2); 
   cond_5_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_5) - nanmean(instr(cellNum).units_t_hist.nspk_5_min2); 
   cond_6_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_6) - nanmean(instr(cellNum).units_t_hist.nspk_6_min2); 
   cond_7_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_7) - nanmean(instr(cellNum).units_t_hist.nspk_7_min2); 
   cond_8_nspk(cellNum) = nanmean(instr(cellNum).units_t_hist.nspk_8) - nanmean(instr(cellNum).units_t_hist.nspk_8_min2);
end

figure; hold on;
plot(1,cond_1_nspk, '.k'); plot(1,nanmean(cond_1_nspk), '.g', 'MarkerSize', 20);
plot(2,cond_2_nspk, '.k'); plot(2,nanmean(cond_2_nspk), '.g', 'MarkerSize', 20);
plot(3,cond_3_nspk, '.k'); plot(3,nanmean(cond_3_nspk), '.g', 'MarkerSize', 20);
plot(4,cond_4_nspk, '.k'); plot(4,nanmean(cond_4_nspk), '.g', 'MarkerSize', 20);
plot(5,cond_5_nspk, '.k'); plot(5,nanmean(cond_5_nspk), '.g', 'MarkerSize', 20);
plot(6,cond_6_nspk, '.k'); plot(6,nanmean(cond_6_nspk), '.g', 'MarkerSize', 20);
plot(7,cond_7_nspk, '.k'); plot(7,nanmean(cond_7_nspk), '.g', 'MarkerSize', 20);
plot(8,cond_8_nspk, '.k'); plot(8,nanmean(cond_8_nspk), '.g', 'MarkerSize', 20);
hline(0,'--k'); title('instr min2')
set(gca, 'ylim', [-40 40], 'TickDir', 'out'); box off

%% sacc psth
% sacc(19)=[]; % leave out for now, deal with it later. 

% min1
ts=sacc(1).units_t_hist.ts; 
for cellNum = 1:length(sacc)
  
    cond_1_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(1,:) - sacc(cellNum).units_t_hist.rate_pst_min1(1,:));
    cond_2_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(2,:) - sacc(cellNum).units_t_hist.rate_pst_min1(2,:));
    cond_3_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(3,:) - sacc(cellNum).units_t_hist.rate_pst_min1(3,:));
    cond_4_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(4,:) - sacc(cellNum).units_t_hist.rate_pst_min1(4,:));
    cond_5_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(5,:) - sacc(cellNum).units_t_hist.rate_pst_min1(5,:));
    cond_6_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(6,:) - sacc(cellNum).units_t_hist.rate_pst_min1(6,:));
    cond_7_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(7,:) - sacc(cellNum).units_t_hist.rate_pst_min1(7,:));
    cond_8_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(8,:) - sacc(cellNum).units_t_hist.rate_pst_min1(8,:));
end
% raw plot
 figure; hold on; 
 plot(ts,nanmean(cond_1_abs)'); plot(ts,nanmean(cond_2_abs)'); plot(ts,nanmean(cond_3_abs)'); plot(ts,nanmean(cond_4_abs)'); 
 plot(ts,nanmean(cond_5_abs)'); plot(ts,nanmean(cond_6_abs)'); plot(ts,nanmean(cond_7_abs)'); plot(ts,nanmean(cond_8_abs)'); 
 set(gca,'xlim', [-0.150 0.150]); vline(0,'--k')

% sacc nspk
for cellNum = 1:length(sacc)
   cond_1_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_1) - nanmean(sacc(cellNum).units_t_hist.nspk_1_min1); 
   cond_2_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_2) - nanmean(sacc(cellNum).units_t_hist.nspk_2_min1); 
   cond_3_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_3) - nanmean(sacc(cellNum).units_t_hist.nspk_3_min1); 
   cond_4_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_4) - nanmean(sacc(cellNum).units_t_hist.nspk_4_min1); 
   cond_5_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_5) - nanmean(sacc(cellNum).units_t_hist.nspk_5_min1); 
   cond_6_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_6) - nanmean(sacc(cellNum).units_t_hist.nspk_6_min1); 
   cond_7_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_7) - nanmean(sacc(cellNum).units_t_hist.nspk_7_min1); 
   cond_8_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_8) - nanmean(sacc(cellNum).units_t_hist.nspk_8_min1); 
end

figure; hold on;
plot(1,cond_1_nspk, '.k'); plot(1,nanmean(cond_1_nspk), '.g', 'MarkerSize', 20);
plot(2,cond_2_nspk, '.k'); plot(2,nanmean(cond_2_nspk), '.g', 'MarkerSize', 20);
plot(3,cond_3_nspk, '.k'); plot(3,nanmean(cond_3_nspk), '.g', 'MarkerSize', 20);
plot(4,cond_4_nspk, '.k'); plot(4,nanmean(cond_4_nspk), '.g', 'MarkerSize', 20);
plot(5,cond_5_nspk, '.k'); plot(5,nanmean(cond_5_nspk), '.g', 'MarkerSize', 20);
plot(6,cond_6_nspk, '.k'); plot(6,nanmean(cond_6_nspk), '.g', 'MarkerSize', 20);
plot(7,cond_7_nspk, '.k'); plot(7,nanmean(cond_7_nspk), '.g', 'MarkerSize', 20);
plot(8,cond_8_nspk, '.k'); plot(8,nanmean(cond_8_nspk), '.g', 'MarkerSize', 20);
hline(0,'--k');
set(gca, 'ylim', [-40 40], 'TickDir', 'out'); box off
title('sacc min1')

% min2
for cellNum = 1:length(sacc)
  
    cond_1_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(1,:) - sacc(cellNum).units_t_hist.rate_pst_min2(1,:));
    cond_2_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(2,:) - sacc(cellNum).units_t_hist.rate_pst_min2(2,:));
    cond_3_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(3,:) - sacc(cellNum).units_t_hist.rate_pst_min2(3,:));
    cond_4_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(4,:) - sacc(cellNum).units_t_hist.rate_pst_min2(4,:));
    cond_5_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(5,:) - sacc(cellNum).units_t_hist.rate_pst_min2(5,:));
    cond_6_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(6,:) - sacc(cellNum).units_t_hist.rate_pst_min2(6,:));
    cond_7_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(7,:) - sacc(cellNum).units_t_hist.rate_pst_min2(7,:));
    cond_8_abs(cellNum,:) = abs(sacc(cellNum).units_t_hist.rate_pst(8,:) - sacc(cellNum).units_t_hist.rate_pst_min2(8,:));
end
% raw plot
 figure; hold on; 
 plot(ts,nanmean(cond_1_abs)'); plot(ts,nanmean(cond_2_abs)'); plot(ts,nanmean(cond_3_abs)'); plot(ts,nanmean(cond_4_abs)'); 
 plot(ts,nanmean(cond_5_abs)'); plot(ts,nanmean(cond_6_abs)'); plot(ts,nanmean(cond_7_abs)'); plot(ts,nanmean(cond_8_abs)'); 
 set(gca,'xlim', [-0.150 0.150]); vline(0,'--k'); title('sacc min2')

% sacc nspk
for cellNum = 1:length(sacc)
   cond_1_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_1) - nanmean(sacc(cellNum).units_t_hist.nspk_1_min2); 
   cond_2_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_2) - nanmean(sacc(cellNum).units_t_hist.nspk_2_min2); 
   cond_3_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_3) - nanmean(sacc(cellNum).units_t_hist.nspk_3_min2); 
   cond_4_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_4) - nanmean(sacc(cellNum).units_t_hist.nspk_4_min2); 
   cond_5_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_5) - nanmean(sacc(cellNum).units_t_hist.nspk_5_min2); 
   cond_6_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_6) - nanmean(sacc(cellNum).units_t_hist.nspk_6_min2); 
   cond_7_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_7) - nanmean(sacc(cellNum).units_t_hist.nspk_7_min2); 
   cond_8_nspk(cellNum) = nanmean(sacc(cellNum).units_t_hist.nspk_8) - nanmean(sacc(cellNum).units_t_hist.nspk_8_min2); 
end

figure; hold on;
plot(1,cond_1_nspk, '.k'); plot(1,nanmean(cond_1_nspk), '.g', 'MarkerSize', 20);
plot(2,cond_2_nspk, '.k'); plot(2,nanmean(cond_2_nspk), '.g', 'MarkerSize', 20);
plot(3,cond_3_nspk, '.k'); plot(3,nanmean(cond_3_nspk), '.g', 'MarkerSize', 20);
plot(4,cond_4_nspk, '.k'); plot(4,nanmean(cond_4_nspk), '.g', 'MarkerSize', 20);
plot(5,cond_5_nspk, '.k'); plot(5,nanmean(cond_5_nspk), '.g', 'MarkerSize', 20);
plot(6,cond_6_nspk, '.k'); plot(6,nanmean(cond_6_nspk), '.g', 'MarkerSize', 20);
plot(7,cond_7_nspk, '.k'); plot(7,nanmean(cond_7_nspk), '.g', 'MarkerSize', 20);
plot(8,cond_8_nspk, '.k'); plot(8,nanmean(cond_8_nspk), '.g', 'MarkerSize', 20);
hline(0,'--k'); title('sacc min2')
set(gca, 'ylim', [-40 40], 'TickDir', 'out'); box off




end

