function trial_hist= trial_history(units,cellNum)
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
trial_hist          = zeros(size(units(cellNum).trial.behav,2),1) ; % with that don't have have three consesequtive correct trials are 0

% starting from trial 3 cause trial 1 and 2 don't have a history
for i = 3:size(units(cellNum).trial.behav,2)
    
    
    these_trials = units(cellNum).trial.behav(i-2:i);
    if sum([these_trials.correctResponse]) == 6 % all 3 trials have correct responses
        
        code(:,1) = ismember(units(cellNum).trial.behav(i-2).conditionCode,proLegend(1,:));
        code(:,2) = ismember(units(cellNum).trial.behav(i-1 ).conditionCode,proLegend(1,:));
        code(:,3) = ismember(units(cellNum).trial.behav( i).conditionCode,proLegend(1,:));
        
        
        trial_hist(i,:) = find(ismember(pro_anti_encoder,code,'rows'));
        
        
    end
end

% extract neural responses for those trials
default_prs_pro_anti; % load prs 

for ii = 1:length(trial_hist)
    if trial_hist(ii) ~= 0
        instr.units_t_hist.code(ii) = trial_hist(ii);
        analyse_sacc_align=0; id = 'SS'
        % instr
        cond = trial_hist(ii);
        
        if cond == 1
            [instr.units_t_hist(cellNum).rate_pst_1, instr.units_t_hist(cellNum).instr.ts_pst_1] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
            instr.units_t_hist(cellNum).rate_pst_1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_1,prs.binwidth,prs.tsmooth); 
            instr.units_t_hist(cellNum).ntrls_1 = length(find(trial_hist==cond));
            % extract n-1 and n-2 trials
            [instr.units_t_hist(cellNum).rate_pst_2_min1, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_2_min1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_2_min1,prs.binwidth,prs.tsmooth); 
            
            [instr.units_t_hist(cellNum).rate_pst_2_min2, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_2_min2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_2_min2,prs.binwidth,prs.tsmooth); 
        elseif cond == 2
            [instr.units_t_hist(cellNum).rate_pst_2, instr.units_t_hist(cellNum).instr.ts_pst_2] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==2)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_2,prs.binwidth,prs.tsmooth); 
            instr.units_t_hist(cellNum).ntrls_2 = length(find(trial_hist==cond));
            % extract n-1 and n-2 trials
            [instr.units_t_hist(cellNum).rate_pst_2_min1, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_2_min1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_2_min1,prs.binwidth,prs.tsmooth); 
            
            [instr.units_t_hist(cellNum).rate_pst_2_min2, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_2_min2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_2_min2,prs.binwidth,prs.tsmooth);
        elseif cond == 3
            [instr.units_t_hist(cellNum).rate_pst_3, instr.units_t_hist(cellNum).instr.ts_pst_3] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==3)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_3 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_4,prs.binwidth,prs.tsmooth); 
            instr.units_t_hist(cellNum).ntrls_3 = length(find(trial_hist==cond));
            % extract n-1 and n-2 trials
            [instr.units_t_hist(cellNum).rate_pst_3_min1, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_3_min1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_3_min1,prs.binwidth,prs.tsmooth); 
            
            [instr.units_t_hist(cellNum).rate_pst_3_min2, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_3_min2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_3_min2,prs.binwidth,prs.tsmooth); 
        elseif cond == 4
            [instr.units_t_hist(cellNum).rate_pst_4, instr.units_t_hist(cellNum).instr.ts_pst_4] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==4)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_4 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_4,prs.binwidth,prs.tsmooth); 
            instr.units_t_hist(cellNum).ntrls_4 = length(find(trial_hist==cond));
            % extract n-1 and n-2 trials
            [instr.units_t_hist(cellNum).rate_pst_4_min1, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_4_min1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_4_min1,prs.binwidth,prs.tsmooth); 
            
            [instr.units_t_hist(cellNum).rate_pst_4_min2, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_4_min2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_4_min2,prs.binwidth,prs.tsmooth); 
            
        elseif cond == 5
            [instr.units_t_hist(cellNum).rate_pst_5, instr.units_t_hist(cellNum).instr.ts_pst_5] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_5 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_5,prs.binwidth,prs.tsmooth); 
            instr.units_t_hist(cellNum).ntrls_5 = length(find(trial_hist==cond));
            % extract n-1 and n-2 trials
            [instr.units_t_hist(cellNum).rate_pst_5_min1, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_5_min1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_5_min1,prs.binwidth,prs.tsmooth); 
            
            [instr.units_t_hist(cellNum).rate_pst_5_min2, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_5_min2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_5_min2,prs.binwidth,prs.tsmooth); 
        elseif cond == 6
            [instr.units_t_hist(cellNum).rate_pst_6, instr.units_t_hist(cellNum).instr.ts_pst_6] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_6 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_6,prs.binwidth,prs.tsmooth); 
            instr.units_t_hist(cellNum).ntrls_6 = length(find(trial_hist==cond));
            % extract n-1 and n-2 trials
            [instr.units_t_hist(cellNum).rate_pst_6_min1, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_6_min1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_6_min1,prs.binwidth,prs.tsmooth); 
            
            [instr.units_t_hist(cellNum).rate_pst_6_min2, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_6_min2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_6_min2,prs.binwidth,prs.tsmooth); 
        elseif cond == 7
            [instr.units_t_hist(cellNum).rate_pst_7, instr.units_t_hist(cellNum).instr.ts_pst_7] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_7 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_7,prs.binwidth,prs.tsmooth); 
            instr.units_t_hist(cellNum).ntrls_7 = length(find(trial_hist==cond));
            % extract n-1 and n-2 trials
            [instr.units_t_hist(cellNum).rate_pst_7_min1, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_7_min1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_7_min1,prs.binwidth,prs.tsmooth); 
            
            [instr.units_t_hist(cellNum).rate_pst_7_min2, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_7_min2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_7_min2,prs.binwidth,prs.tsmooth); 
        else
            [instr.units_t_hist(cellNum).rate_pst_8, instr.units_t_hist(cellNum).instr.ts_pst_8] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_8 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_8,prs.binwidth,prs.tsmooth); 
            instr.units_t_hist(cellNum).ntrls_8 = length(find(trial_hist==cond));
            % extract n-1 and n-2 trials
            [instr.units_t_hist(cellNum).rate_pst_8_min1, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-1),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_8_min1 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_8_min1,prs.binwidth,prs.tsmooth); 
            
            [instr.units_t_hist(cellNum).rate_pst_8_min2, ~] = Spiketimes2Rate(units(cellNum).trial.neural(find(trial_hist==cond)-2),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
            instr.units_t_hist(cellNum).rate_pst_8_min2 = smooth_pst(instr.units_t_hist(cellNum).rate_pst_8_min2,prs.binwidth,prs.tsmooth); 
        end
        
    end
end

%% sacc

analyse_sacc_align=1;
[units(cellNum).anti.neural.sacc.rate_pst,units(cellNum).anti.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to saccade onset
units(cellNum).anti.neural.sacc.rate_pst = smooth_pst(units(cellNum).anti.neural.sacc.rate_pst,prs.binwidth,prs.tsmooth);

analyse_sacc_align=0;
end

