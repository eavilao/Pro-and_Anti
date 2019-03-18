function units = extractWholeNeuronR(wholeNeuronResults)
% extractWholeNeuronResults
% This function extracts data from wholeNeuronResults from Pete's code and
% performs stats on simple spikes.

% In: wholeNeuronResults.m
% Out:  units containing spike times (Simple and Complex spikes). Every
% struct inside units is a neuron.
% aligned to trial onset and eye stats.

% Associated functions: SpikeTimes2Rate.m, SpikeTimes2RateTrial.m,
%                       SpikeTimes2CountTrial.m
%                       plot_ProAnti.m, eyeKinematics_ProAnti.m

% Windows
% Baseline = -0.3:-0.1 before saccade onset
% Instruction = 0:0.3 aligned to trial onset
% Saccade = -0.1:0.2 aligned to saccade

%% Create structure with relevant data (units) - extract main info from wholeNeuronResults
default_prs_pro_anti; % load parameters list
run_error_trials =0 ; % if run_error_trials = 1 - extract only error trials same way as normal, min trial number from default prs

% first find the ones that are not empty and
for i = 1:length(wholeNeuronResults)
    if  ~isempty(wholeNeuronResults(i).allStableTrials);
        cell_indx(i) = numel(wholeNeuronResults(i).selectedTrials.corProTrials)>=prs.min_trial & numel(wholeNeuronResults(i).selectedTrials.corAntiTrials)>=prs.min_trial;
    end
end
cell_indx = find(cell_indx);


for cellNum = 1:length(cell_indx);
    trials_with_spk=zeros(1,length(wholeNeuronResults(cell_indx(cellNum)).allStableTrials));
    
    for trialNum = 1:length(wholeNeuronResults(cell_indx(cellNum)).allStableTrials);
        
        if ~isempty(wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).alignedSpikes)
            trials_with_spk(trialNum)=1;
            % monkey
            units(cellNum).monk = wholeNeuronResults(cell_indx(cellNum)).monkey;
            % area
            units(cellNum).area = wholeNeuronResults(cell_indx(cellNum)).area;
            % coord
            units(cellNum).coord.depth = wholeNeuronResults(cell_indx(cellNum)).depth;
            units(cellNum).coord.loc = wholeNeuronResults(cell_indx(cellNum)).gridLoc;
            %behavioral
            units(cellNum).trial.behav(trialNum).saccadeOnset =  wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).saccadeTime; % saccade onset from trial start (?)
            units(cellNum).trial.behav(trialNum).goCueTime = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).goCueTime; %go cue from trial start
            units(cellNum).trial.behav(trialNum).saccAmplitude = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).saccadeAmplitude;
            units(cellNum).trial.behav(trialNum).saccDuration = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).saccadeDuration;
            units(cellNum).trial.behav(trialNum).saccPeakVel = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).saccadePeakVel;
            units(cellNum).trial.behav(trialNum).reactionTime = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).reactionTime;
            units(cellNum).trial.behav(trialNum).reward = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).relativeRewardTime;
            units(cellNum).trial.behav(trialNum).conditionCode = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).conditionCode;
            units(cellNum).trial.behav(trialNum).correctResponse = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).correctResponse;
            units(cellNum).trial.behav(trialNum).eyePos_x = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).eyePositionX;
            units(cellNum).trial.behav(trialNum).eyePos_y = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).eyePositionY;
            units(cellNum).trial.behav(trialNum).eyeVel_x = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).eyeVelocityX;
            units(cellNum).trial.behav(trialNum).eyeVel_y = wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).eyeVelocityY;
            units(cellNum).trial.behav(trialNum).eye_ts = [-3999:1:4000]/1000;
            units(cellNum).trial.behav(trialNum).eye_ts_sacc = units(cellNum).trial.behav(trialNum).eye_ts(units(cellNum).trial.behav(trialNum).eye_ts >= prs.eye_win_sacc(1) & ...
                units(cellNum).trial.behav(trialNum).eye_ts <= prs.eye_win_sacc(2));
            units(cellNum).trial.behav(trialNum).eyePos_x_sacc = units(cellNum).trial.behav(trialNum).eyePos_x(units(cellNum).trial.behav(trialNum).eye_ts >= prs.eye_win_sacc(1) & ...
                units(cellNum).trial.behav(trialNum).eye_ts <= prs.eye_win_sacc(2));
            units(cellNum).trial.behav(trialNum).eyePos_y_sacc = units(cellNum).trial.behav(trialNum).eyePos_y(units(cellNum).trial.behav(trialNum).eye_ts >= prs.eye_win_sacc(1) & ...
                units(cellNum).trial.behav(trialNum).eye_ts <= prs.eye_win_sacc(2));
            units(cellNum).trial.behav(trialNum).eyeVel_x_sacc = units(cellNum).trial.behav(trialNum).eyeVel_x(units(cellNum).trial.behav(trialNum).eye_ts >= prs.eye_win_sacc(1) & ...
                units(cellNum).trial.behav(trialNum).eye_ts <= prs.eye_win_sacc(2));
            units(cellNum).trial.behav(trialNum).eyeVel_y_sacc = units(cellNum).trial.behav(trialNum).eyeVel_y(units(cellNum).trial.behav(trialNum).eye_ts >= prs.eye_win_sacc(1) & ...
                units(cellNum).trial.behav(trialNum).eye_ts <= prs.eye_win_sacc(2));
            
            %neural data
            spks =  wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).alignedSpikes{1}; % contains spike times for SS aligned to trial onset
            if ~isempty(units(cellNum).trial.behav(trialNum).reward) & ~isnan(units(cellNum).trial.behav(trialNum).reward)
                
                indx_trial_spks = spks>prs.tspk(1) & spks < units(cellNum).trial.behav(trialNum).reward+prs.tspk(2); % just pick 100 ms before trial starts to reward +200 ms
                units(cellNum).trial.neural(trialNum).tspk_SS = spks(indx_trial_spks);
                units(cellNum).trial.neural(trialNum).tspk_SS_align_sacc =  units(cellNum).trial.neural(trialNum).tspk_SS-units(cellNum).trial.behav(trialNum).saccadeOnset; % contains spike times for CS aligned to sacc onset
                units(cellNum).trial.neural(trialNum).ts_eye_sacc = units(cellNum).trial.behav(trialNum).eye_ts-units(cellNum).trial.behav(trialNum).saccadeOnset;
            elseif run_error_trials
                
                indx_trial_spks = spks>prs.tspk(1) & spks < units(cellNum).trial.behav(trialNum).goCueTime+0.5+prs.tspk(2); % just pick 100 ms before trial starts to reward +200 ms
                units(cellNum).trial.neural(trialNum).tspk_SS = spks(indx_trial_spks);
                units(cellNum).trial.neural(trialNum).tspk_SS_align_sacc =  units(cellNum).trial.neural(trialNum).tspk_SS-units(cellNum).trial.behav(trialNum).saccadeOnset; % contains spike times for CS aligned to sacc onset
            else
                %keyboard
                units(cellNum).trial.neural(trialNum).tspk_SS = [];
                units(cellNum).trial.neural(trialNum).tspk_SS_align_sacc = [];
            end
            units(cellNum).id = 'SS';
            % select condition type Pro or Antisaccade
            if ismember(wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).conditionCode, prs.proConditions);
                units(cellNum).trial.behav(trialNum).condition = 'Prosaccade';
            else
                units(cellNum).trial.behav(trialNum).condition = 'Antisaccade';
            end
        end
        
    end
    % gather incorrect pro and anti trials
    response =   [wholeNeuronResults(cell_indx(cellNum)).allStableTrials.correctResponse];
    conditionCode = [wholeNeuronResults(cell_indx(cellNum)).allStableTrials.conditionCode];
    
    anti = ismember(conditionCode,  prs.antiConditions);
    pro = ismember(conditionCode, prs.proConditions );
    
    incorrect = ismember(response, 1);
    
    incorrectProTrials = find((incorrect+pro+trials_with_spk)==3)  ;
    incorrectAntiTrials = find((incorrect+anti+trials_with_spk)==3)  ;
    
    % Indexes to select Pro and Anti trials
    units(cellNum).pro.indx_correctProTrials = wholeNeuronResults(cell_indx(cellNum)).selectedTrials.corProTrials;
    units(cellNum).anti.indx_correctAntiTrials = wholeNeuronResults(cell_indx(cellNum)).selectedTrials.corAntiTrials;
    
    units(cellNum).pro.indx_incorrectProTrials = incorrectProTrials;
    units(cellNum).anti.indx_incorrectAntiTrials = incorrectAntiTrials;
    
    
    % pretend error trials are normal trials to not have to rewrite
    % everything
    if run_error_trials
        if numel(incorrectProTrials)<prs.min_trial |  numel(incorrectAntiTrials)<prs.min_trial
            throw(cellNum)=1;
        end
        
        units(cellNum).pro.indx_correctProTrials = incorrectProTrials;
        units(cellNum).anti.indx_correctAntiTrials = incorrectAntiTrials;
        
        
    end
end
clear cellNum trialNum %wholeNeuronResults

%% Spike times to rate for pro and antisaccades and behav
analyse_sacc_align = 0;

% throw out units with no error trials in pro or anti f
if run_error_trials
    units(find(throw))=[];
end

for cellNum = 1:length(units)
    id=units(cellNum).id;
    % Pro trials
    correctProTrials = units(cellNum).pro.indx_correctProTrials; % index pro trials
    % behav
    for trialNum = 1:length(units(cellNum).pro.indx_correctProTrials)
        units(cellNum).pro.behav.trial(trialNum).saccadeOnset = units(cellNum).trial.behav(correctProTrials(trialNum)).saccadeOnset;
        units(cellNum).pro.behav.trial(trialNum).goCueTime = units(cellNum).trial.behav(correctProTrials(trialNum)).goCueTime;
        units(cellNum).pro.behav.trial(trialNum).saccAmplitude = units(cellNum).trial.behav(correctProTrials(trialNum)).saccAmplitude;
        units(cellNum).pro.behav.trial(trialNum).saccDuration = units(cellNum).trial.behav(correctProTrials(trialNum)).saccDuration;
        units(cellNum).pro.behav.trial(trialNum).saccPeakVel = units(cellNum).trial.behav(correctProTrials(trialNum)).saccPeakVel;
        units(cellNum).pro.behav.trial(trialNum).reactionTime = units(cellNum).trial.behav(correctProTrials(trialNum)).reactionTime;
        units(cellNum).pro.behav.trial(trialNum).reward = units(cellNum).trial.behav(correctProTrials(trialNum)).reward;
        units(cellNum).pro.behav.trial(trialNum).conditionCode = units(cellNum).trial.behav(correctProTrials(trialNum)).conditionCode;
        units(cellNum).pro.behav.trial(trialNum).eye_ts_sacc = units(cellNum).trial.behav(correctProTrials(trialNum)).eye_ts_sacc;
        units(cellNum).pro.behav.trial(trialNum).eyePos_x_sacc = units(cellNum).trial.behav(correctProTrials(trialNum)).eyePos_x_sacc;
        units(cellNum).pro.behav.trial(trialNum).eyePos_y_sacc = units(cellNum).trial.behav(correctProTrials(trialNum)).eyePos_y_sacc;
        units(cellNum).pro.behav.trial(trialNum).eyeVel_x_sacc = units(cellNum).trial.behav(correctProTrials(trialNum)).eyeVel_x_sacc;
        units(cellNum).pro.behav.trial(trialNum).eyeVel_y_sacc = units(cellNum).trial.behav(correctProTrials(trialNum)).eyeVel_y_sacc;
        
    end
    
    
    
    % neural
    % instr
    units(cellNum).pro.neural.trial = units(cellNum).trial.neural(units(cellNum).pro.indx_correctProTrials);
    %instr
    try
        [units(cellNum).pro.neural.instr.rate_pst,units(cellNum).pro.neural.instr.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
    catch
        keyboard
    end
    
    units(cellNum).pro.neural.instr.rate_pst = smooth_pst(units(cellNum).pro.neural.instr.rate_pst,prs.binwidth,prs.tsmooth);
    % sacc aligned
    analyse_sacc_align=1;
    [units(cellNum).pro.neural.sacc.rate_pst,units(cellNum).pro.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to saccade onset
    units(cellNum).pro.neural.sacc.rate_pst = smooth_pst(units(cellNum).pro.neural.sacc.rate_pst,prs.binwidth,prs.tsmooth);
    
    %% sanity check -- for every cell, ranomize half of the trials and compute psth
    numtrials = []; pick_trials = []; 
    numtrials = length(units(cellNum).pro.neural.trial); % number of trials
    pick_trials = sort(randsample(numtrials, round(numtrials/2)));  % pick half of the trials randomly
    
    [units(cellNum).pro.neural.sacc.rate_pst_rand,~] = Spiketimes2Rate(units(cellNum).pro.neural.trial((pick_trials)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to saccade onset
    units(cellNum).pro.neural.sacc.rate_pst_rand = smooth_pst(units(cellNum).pro.neural.sacc.rate_pst_rand,prs.binwidth,prs.tsmooth);
    analyse_sacc_align=0;
    
    % Anti trials
    correctAntiTrials = units(cellNum).anti.indx_correctAntiTrials; % index pro trials
    % behav
    for trialNum = 1:length(units(cellNum).anti.indx_correctAntiTrials)
        units(cellNum).anti.behav.trial(trialNum).saccadeOnset = units(cellNum).trial.behav(correctAntiTrials(trialNum)).saccadeOnset;
        units(cellNum).anti.behav.trial(trialNum).goCueTime = units(cellNum).trial.behav(correctAntiTrials(trialNum)).goCueTime;
        units(cellNum).anti.behav.trial(trialNum).saccAmplitude = units(cellNum).trial.behav(correctAntiTrials(trialNum)).saccAmplitude;
        units(cellNum).anti.behav.trial(trialNum).saccDuration = units(cellNum).trial.behav(correctAntiTrials(trialNum)).saccDuration;
        units(cellNum).anti.behav.trial(trialNum).saccPeakVel = units(cellNum).trial.behav(correctAntiTrials(trialNum)).saccPeakVel;
        units(cellNum).anti.behav.trial(trialNum).reactionTime = units(cellNum).trial.behav(correctAntiTrials(trialNum)).reactionTime;
        units(cellNum).anti.behav.trial(trialNum).reward = units(cellNum).trial.behav(correctAntiTrials(trialNum)).reward;
        units(cellNum).anti.behav.trial(trialNum).conditionCode = units(cellNum).trial.behav(correctAntiTrials(trialNum)).conditionCode;
        units(cellNum).anti.behav.trial(trialNum).eye_ts_sacc = units(cellNum).trial.behav(correctAntiTrials(trialNum)).eye_ts_sacc;
        units(cellNum).anti.behav.trial(trialNum).eyePos_x_sacc = units(cellNum).trial.behav(correctAntiTrials(trialNum)).eyePos_x_sacc;
        units(cellNum).anti.behav.trial(trialNum).eyePos_y_sacc = units(cellNum).trial.behav(correctAntiTrials(trialNum)).eyePos_y_sacc;
        units(cellNum).anti.behav.trial(trialNum).eyeVel_x_sacc = units(cellNum).trial.behav(correctAntiTrials(trialNum)).eyeVel_x_sacc;
        units(cellNum).anti.behav.trial(trialNum).eyeVel_y_sacc = units(cellNum).trial.behav(correctAntiTrials(trialNum)).eyeVel_y_sacc;
    end
    %neural
    units(cellNum).anti.neural.trial = units(cellNum).trial.neural(units(cellNum).anti.indx_correctAntiTrials);
    % intr
    [units(cellNum).anti.neural.instr.rate_pst,units(cellNum).anti.neural.instr.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id); % aligned to trial onset
    units(cellNum).anti.neural.instr.rate_pst = smooth_pst(units(cellNum).anti.neural.instr.rate_pst,prs.binwidth,prs.tsmooth);
    % sacc aligned
    analyse_sacc_align=1;
    [units(cellNum).anti.neural.sacc.rate_pst,units(cellNum).anti.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to saccade onset
    units(cellNum).anti.neural.sacc.rate_pst = smooth_pst(units(cellNum).anti.neural.sacc.rate_pst,prs.binwidth,prs.tsmooth);
    
    %% sanity check -- for every cell, pick random half of the trials and compute psth
    numtrials = []; pick_trials = []; 
    numtrials = length(units(cellNum).anti.neural.trial); % number of trials
    pick_trials = sort(randsample(numtrials, round(numtrials/2)));  % pick half of the trials randomly
    
    [units(cellNum).anti.neural.sacc.rate_pst_rand,~] = Spiketimes2Rate(units(cellNum).anti.neural.trial((pick_trials)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id); % aligned to saccade onset
    units(cellNum).anti.neural.sacc.rate_pst_rand = smooth_pst(units(cellNum).anti.neural.sacc.rate_pst_rand,prs.binwidth,prs.tsmooth);

    analyse_sacc_align=0;
end

%% For every trial - Spike time count and spike times to rate for pro and antisaccades

for cellNum = 1:length(units)
    id=units(cellNum).id;
    %pro trials
    correctProTrials = units(cellNum).pro.indx_correctProTrials; % index pro trials
    % get spk counts
    for trialNum = 1:length(correctProTrials)
        
        [units(cellNum).pro.neural.instr.spkCount(trialNum,:),units(cellNum).pro.neural.instr.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).pro.neural.trial(trialNum),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id); % spk counts
        analyse_sacc_align=1;
        [units(cellNum).pro.neural.sacc.spkCount(trialNum,:),units(cellNum).pro.neural.sacc.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).pro.neural.trial(trialNum),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
        analyse_sacc_align=0;
    end
    % Take out probability of spk in pro
    for i=1:length(units(cellNum).pro.neural.instr.spkCount(1,:))
        % instr aligned
        units(cellNum).pro.neural.instr.pbDist(1,i)= sum(units(cellNum).pro.neural.instr.spkCount(:,i))/length(correctProTrials); % compute probability of spk instr onset
        units(cellNum).pro.neural.instr.pbDist_sem = sqrt(units(cellNum).pro.neural.instr.pbDist(1,i)*(1-units(cellNum).pro.neural.instr.pbDist(1,i))/length(correctProTrials));
        % sacc aligned
        units(cellNum).pro.neural.sacc.pbDist(1,i)= sum(units(cellNum).pro.neural.sacc.spkCount(:,i))/length(correctProTrials); % compute probability of spk sacc aligned
        units(cellNum).pro.neural.sacc.pbDist_sem = sqrt(units(cellNum).pro.neural.sacc.pbDist(1,i)*(1-units(cellNum).pro.neural.sacc.pbDist(1,i))/length(correctProTrials));
    end
    
    
    % anti trials
    correctAntiTrials = units(cellNum).anti.indx_correctAntiTrials;
    for trialNum = 1:length(correctAntiTrials)
        [units(cellNum).anti.neural.instr.spkCount(trialNum,:),units(cellNum).anti.neural.instr.trial(trialNum).ts] = Spiketimes2CountTrial(units(cellNum).anti.neural.trial(trialNum),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id); % spk counts
        analyse_sacc_align=1;
        [units(cellNum).anti.neural.sacc.spkCount(trialNum,:),units(cellNum).anti.neural.sacc.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).anti.neural.trial(trialNum),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
        analyse_sacc_align=0;
        
    end
    
    % Take out probability of spk in anti
    for i=1:length(units(cellNum).anti.neural.instr.spkCount(1,:))
        % instr aligned
        units(cellNum).anti.neural.instr.pbDist(1,i)= sum(units(cellNum).anti.neural.instr.spkCount(:,i))/length(correctAntiTrials); % compute probability of spk instr onset
        units(cellNum).anti.neural.instr.pbDist_sem = sqrt(units(cellNum).pro.neural.instr.pbDist(1,i)*(1-units(cellNum).anti.neural.instr.pbDist(1,i))/length(correctAntiTrials));
        % sacc aligned
        units(cellNum).anti.neural.sacc.pbDist(1,i)= sum(units(cellNum).anti.neural.sacc.spkCount(:,i))/length(correctAntiTrials); % compute probability of spk sacc aligned
        units(cellNum).anti.neural.sacc.pbDist_sem = sqrt(units(cellNum).anti.neural.sacc.pbDist(1,i)*(1-units(cellNum).anti.neural.sacc.pbDist(1,i))/length(correctAntiTrials));
    end
    
    % nspk pro - all bins and windows per trial
    units(cellNum).pro.neural.instr.tspk_all = []; units(cellNum).pro.neural.sacc.tspk_all = []; 
    for trialNum = 1:length(correctProTrials)
        
        timepoints_instr = -1:prs.binwidth:1; timepoints_sacc = -1:prs.binwidth:1;
        %         timepoints_instr = -0.1:binwidth:1; timepoints_sacc = -0.8:binwidth:0.3;
        
        % align to instr
        [~,units(cellNum).pro.neural.instr.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum).tspk_SS,prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
        
        % win
        tspk_instr= units(cellNum).pro.neural.trial(trialNum).tspk_SS>prs.instruction_win(1) &units(cellNum).pro.neural.trial(trialNum).tspk_SS<prs.instruction_win(2);
        t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
        units(cellNum).pro.neural.instr.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS(tspk_instr);
        [units(cellNum).pro.neural.instr.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.instr.tspk(trialNum,:)),timepoints_instr,prs.binwidth,analyse_sacc_align,id);
        
        % align to sacc
        analyse_sacc_align = 1;
        [~,units(cellNum).pro.neural.sacc.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
        
        % win
        tspk_sacc= units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc>prs.saccade_win(1) & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc<prs.saccade_win(2);
        t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
        units(cellNum).pro.neural.sacc.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc(tspk_sacc);
        [units(cellNum).pro.neural.sacc.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.sacc.tspk(trialNum,:)),timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
        analyse_sacc_align = 0;
        
        % concatenate tspk 
        units(cellNum).pro.neural.instr.tspk_all = [units(cellNum).pro.neural.instr.tspk_all ; units(cellNum).pro.neural.instr.tspk{trialNum,:}];
        units(cellNum).pro.neural.sacc.tspk_all = [units(cellNum).pro.neural.sacc.tspk_all ; units(cellNum).pro.neural.sacc.tspk{trialNum,:}];
        
        % baseline
        tspk_base = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc>prs.baseline_win(1) & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc<prs.baseline_win(2);
        t_base = prs.timepoints_sacc >prs.baseline_win(1) & prs.timepoints_sacc<prs.baseline_win(2); timepoints_base = prs.timepoints_sacc(t_base);
        units(cellNum).pro.neural.base.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc(tspk_base);
        [units(cellNum).pro.neural.base.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.base.tspk(trialNum,:)),timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
        analyse_sacc_align = 0;
    end
    
    % nspk anti - windows per trial
    units(cellNum).anti.neural.instr.tspk_all = []; units(cellNum).anti.neural.sacc.tspk_all = []; 
    for trialNum = 1:length(correctAntiTrials)
        % align to instr
        [~,units(cellNum).anti.neural.instr.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum).tspk_SS,prs.timepoints_instr,prs.binwidth,analyse_sacc_align,id);
        
        % win
        tspk_instr= units(cellNum).anti.neural.trial(trialNum).tspk_SS>prs.instruction_win(1) &units(cellNum).anti.neural.trial(trialNum).tspk_SS<prs.instruction_win(2);
        t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
        units(cellNum).anti.neural.instr.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS(tspk_instr);
        [units(cellNum).anti.neural.instr.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.instr.tspk(trialNum,:)),timepoints_instr,prs.binwidth,analyse_sacc_align,id);
        
        % align to sacc
        analyse_sacc_align = 1;
        [~,units(cellNum).anti.neural.sacc.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
        
        % win
        tspk_sacc= units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc>prs.saccade_win(1) &units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc<prs.saccade_win(2);
        t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc < prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
        units(cellNum).anti.neural.sacc.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc(tspk_sacc);
        [units(cellNum).anti.neural.sacc.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.sacc.tspk(trialNum,:)),timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
        analyse_sacc_align = 0;
        
          % concatenate tspk 
        units(cellNum).anti.neural.instr.tspk_all = [units(cellNum).anti.neural.instr.tspk_all ; units(cellNum).anti.neural.instr.tspk{trialNum,:}];
        units(cellNum).anti.neural.sacc.tspk_all = [units(cellNum).anti.neural.sacc.tspk_all ; units(cellNum).anti.neural.sacc.tspk{trialNum,:}];
        
        % baseline
        tspk_base = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc>prs.baseline_win(1) & units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc<prs.baseline_win(2);
        t_base = prs.timepoints_sacc >prs.baseline_win(1) & prs.timepoints_sacc<prs.baseline_win(2); timepoints_base = prs.timepoints_sacc(t_base);
        units(cellNum).anti.neural.base.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc(tspk_base);
        [units(cellNum).anti.neural.base.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.base.tspk(trialNum,:)),timepoints_sacc,prs.binwidth,analyse_sacc_align,id);
        analyse_sacc_align = 0;
        
    end
    
    
    %% probability of spk in pro and anti  - bigger window (5 bins = 50 ms)
    % spike probability density for 50 ms
    
    win_size = prs.win_size;
    % instr pro
    for trialNum = 1:length(correctProTrials)
        thisTrial = units(cellNum).pro.neural.instr.nspkCount(trialNum,:);
        indx_beg = 1;indx_win_count = win_size;
        for win_num = 1:length(thisTrial(1,:))-win_size
            units(cellNum).pro.neural.instr.spkCount_win(trialNum,win_num) = sum(thisTrial(indx_beg:indx_win_count))/win_size;
            units(cellNum).pro.neural.instr.ts_spkCount_win(win_num,:) = units(cellNum).pro.neural.instr.ts_pst(indx_beg:indx_win_count);
            indx_beg = indx_beg+1; indx_win_count = indx_win_count+1;
        end
    end
    
    % instr anti
    for trialNum = 1:length(correctAntiTrials)
        thisTrial = units(cellNum).anti.neural.instr.nspkCount(trialNum,:);
        indx_beg = 1;indx_win_count = win_size;
        for win_num = 1:length(thisTrial(1,:))-win_size
            units(cellNum).anti.neural.instr.spkCount_win(trialNum,win_num) = sum(thisTrial(indx_beg:indx_win_count))/win_size;
            units(cellNum).anti.neural.instr.ts_spkCount_win(win_num,:) = units(cellNum).anti.neural.instr.ts_pst(indx_beg:indx_win_count);
            indx_beg = indx_beg+1; indx_win_count = indx_win_count+1;
        end
    end
    
    % sacc pro
    for trialNum = 1:length(correctProTrials)
        thisTrial = units(cellNum).pro.neural.sacc.nspkCount(trialNum,:);
        indx_beg = 1;indx_win_count = win_size;
        for win_num = 1:length(thisTrial(1,:))-win_size
            units(cellNum).pro.neural.sacc.spkCount_win(trialNum,win_num) = sum(thisTrial(indx_beg:indx_win_count))/win_size;
            units(cellNum).pro.neural.sacc.ts_spkCount_win(win_num,:) = units(cellNum).pro.neural.sacc.ts_pst(indx_beg:indx_win_count);
            indx_beg = indx_beg+1; indx_win_count = indx_win_count+1;
        end
    end
    
    % sacc anti
    for trialNum = 1:length(correctAntiTrials)
        thisTrial = units(cellNum).anti.neural.sacc.nspkCount(trialNum,:);
        indx_beg = 1;indx_win_count = win_size;
        for win_num = 1:length(thisTrial(1,:))-win_size
            units(cellNum).anti.neural.sacc.spkCount_win(trialNum,win_num) = sum(thisTrial(indx_beg:indx_win_count))/win_size;
            units(cellNum).anti.neural.sacc.ts_spkCount_win(win_num,:) = units(cellNum).anti.neural.sacc.ts_pst(indx_beg:indx_win_count);
            indx_beg = indx_beg+1; indx_win_count = indx_win_count+1;
        end
    end
    
end

%% STATS %%
%%%%%%%%%%%
%%%%%%%%%%%
%% Eye kinematics
eyeKin = eyeKinematics_ProAnti(units); % extract

proAmp = vertcat(eyeKin(1,:).proAmp);
antiAmp = vertcat(eyeKin(1,:).antiAmp);
proDur = vertcat(eyeKin(1,:).proDur);
antiDur = vertcat(eyeKin(1,:).antiDur);
proPV = vertcat(eyeKin(1,:).proPV);
antiPV = vertcat(eyeKin(1,:).antiPV);
proRT = vertcat(eyeKin(1,:).proRT)*1000;
antiRT = vertcat(eyeKin(1,:).antiRT)*1000;

% stats per cell
for cellNum = 1:length(units)
    %pro
    units(cellNum).pro.stats.eye.amp.mu = nanmean(eyeKin(cellNum).proAmp); units(cellNum).pro.stats.eye.amp.sig = nanstd(eyeKin(cellNum).proAmp);
    units(cellNum).pro.stats.eye.dur.mu = nanmean(eyeKin(cellNum).proDur); units(cellNum).pro.stats.eye.dur.sig = nanstd(eyeKin(cellNum).proDur);
    units(cellNum).pro.stats.eye.pv.mu = nanmean(eyeKin(cellNum).proPV); units(cellNum).pro.stats.eye.pv.sig = nanstd(eyeKin(cellNum).proPV);
    units(cellNum).pro.stats.eye.rt.mu = nanmean(eyeKin(cellNum).proRT); units(cellNum).pro.stats.eye.rt.sig = nanstd(eyeKin(cellNum).proRT);
    %anti
    units(cellNum).anti.stats.eye.amp.mu = nanmean(eyeKin(cellNum).antiAmp); units(cellNum).anti.stats.eye.amp.sig = nanstd(eyeKin(cellNum).antiAmp);
    units(cellNum).anti.stats.eye.dur.mu = nanmean(eyeKin(cellNum).antiDur); units(cellNum).anti.stats.eye.dur.sig = nanstd(eyeKin(cellNum).antiDur);
    units(cellNum).anti.stats.eye.pv.mu = nanmean(eyeKin(cellNum).antiPV); units(cellNum).anti.stats.eye.pv.sig = nanstd(eyeKin(cellNum).antiPV);
    units(cellNum).anti.stats.eye.rt.mu = nanmean(eyeKin(cellNum).antiRT); units(cellNum).anti.stats.eye.rt.sig = nanstd(eyeKin(cellNum).antiRT);
end

%% Neural %%%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%

%% stats per neuron - compare windows against baseline activity for pro and anti (concatenate all trials and make comparison)
%% Extract spks in windows and compute mean and sem
for cellNum = 1:length(units)
    ntrls_pro = length(units(cellNum).pro.neural.trial);
    ntrls_anti = length(units(cellNum).anti.neural.trial);
    % windows
    instr_win = units(cellNum).pro.neural.instr.ts_pst > prs.instruction_win(1) & units(cellNum).pro.neural.instr.ts_pst <= prs.instruction_win(2);
    sacc_win = units(cellNum).pro.neural.sacc.ts_pst > prs.saccade_win(1) & units(cellNum).pro.neural.sacc.ts_pst <= prs.saccade_win(2);
    base_win = units(cellNum).pro.neural.sacc.ts_pst > prs.baseline_win(1) & units(cellNum).pro.neural.sacc.ts_pst <= prs.baseline_win(2);
    t_instr= units(cellNum).pro.neural.instr.ts_pst(instr_win); %time
    t_sacc = units(cellNum).pro.neural.sacc.ts_pst(sacc_win);
    sacc_on = units(cellNum).pro.neural.sacc.ts_pst > 0 & units(cellNum).pro.neural.sacc.ts_pst <= prs.saccade_win(2);
    sacc_pre = units(cellNum).pro.neural.sacc.ts_pst > prs.saccade_win(1) & units(cellNum).pro.neural.sacc.ts_pst <= -0.050;
    
    %get spks pro
    instr_spks_pro = units(cellNum).pro.neural.instr.rate_pst(instr_win);
    sacc_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(sacc_win);
    base_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(base_win);
    sacc_on_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(sacc_on);
    sacc_pre_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(sacc_pre);
    
    %pst, mean and sem
    %instr
    units(cellNum).pro.neural.instr.ts_pst_win = units(cellNum).pro.neural.instr.ts_pst(instr_win);
    units(cellNum).pro.neural.instr.rate_pst_win = units(cellNum).pro.neural.instr.rate_pst(instr_win);
    units(cellNum).pro.neural.instr.rate_mu = mean(instr_spks_pro);
    units(cellNum).pro.neural.instr.rate_sig = std(instr_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.instr.rate_std = std(instr_spks_pro);
    
    [units(cellNum).pro.neural.instr.peak_resp, indx_max] = max(units(cellNum).pro.neural.instr.rate_pst_win);
    units(cellNum).pro.neural.instr.peak_resp_time = t_instr(indx_max);
    %sacc
    units(cellNum).pro.neural.sacc.ts_pst_win = units(cellNum).pro.neural.sacc.ts_pst(sacc_win);
    units(cellNum).pro.neural.sacc.rate_pst_win = units(cellNum).pro.neural.sacc.rate_pst(sacc_win);
    units(cellNum).pro.neural.sacc.rate_mu = mean(sacc_spks_pro);
    units(cellNum).pro.neural.sacc.rate_sig = std(sacc_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.sacc.rate_std = std(sacc_spks_pro);
    
    [units(cellNum).pro.neural.sacc.peak_resp, indx_max] = max(units(cellNum).pro.neural.sacc.rate_pst_win);
    units(cellNum).pro.neural.sacc.peak_resp_time = t_sacc(indx_max);
    %base
    units(cellNum).pro.neural.base.ts_pst_win = units(cellNum).pro.neural.sacc.ts_pst(base_win);
    units(cellNum).pro.neural.base.rate_pst_win = units(cellNum).pro.neural.sacc.rate_pst(base_win);
    units(cellNum).pro.neural.base.rate_mu = mean(base_spks_pro);
    units(cellNum).pro.neural.base.rate_sig = std(base_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.base.rate_std = std(base_spks_pro);
    
    % sacc after 0 > or < than baseline
    if mean(sacc_on_spks_pro) > mean(sacc_pre_spks_pro)
    units(cellNum).pro.neural.exc = 1; units(cellNum).pro.neural.sup = 0; 
    else
         units(cellNum).pro.neural.exc = 0; units(cellNum).pro.neural.sup = 1; 
    end
    
    %get spks anti
    instr_spks_anti = units(cellNum).anti.neural.instr.rate_pst(instr_win);
    sacc_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(sacc_win);
    base_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(base_win);
    sacc_on_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(sacc_on);
    sacc_pre_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(sacc_pre);
    
    %pst, mean and sem
    %instr
    units(cellNum).anti.neural.instr.ts_pst_win = units(cellNum).anti.neural.instr.ts_pst(instr_win);
    units(cellNum).anti.neural.instr.rate_pst_win = units(cellNum).anti.neural.instr.rate_pst(instr_win);
    units(cellNum).anti.neural.instr.rate_mu = mean(instr_spks_anti);
    units(cellNum).anti.neural.instr.rate_sig = std(instr_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.instr.rate_std = std(instr_spks_anti);
    
    [units(cellNum).anti.neural.instr.peak_resp, indx_max] = max(units(cellNum).anti.neural.instr.rate_pst_win);
    units(cellNum).anti.neural.instr.peak_resp_time = t_instr(indx_max);
    
    % sacc
    units(cellNum).anti.neural.sacc.ts_pst_win = units(cellNum).anti.neural.sacc.ts_pst(sacc_win);
    units(cellNum).anti.neural.sacc.rate_pst_win = units(cellNum).anti.neural.sacc.rate_pst(sacc_win);
    units(cellNum).anti.neural.sacc.rate_mu = mean(sacc_spks_anti);
    units(cellNum).anti.neural.sacc.rate_sig = std(sacc_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.sacc.rate_std = std(sacc_spks_anti);
    
    [units(cellNum).anti.neural.sacc.peak_resp, indx_max] = max(units(cellNum).anti.neural.sacc.rate_pst_win);
    units(cellNum).anti.neural.sacc.peak_resp_time = t_sacc(indx_max);
    
    % base
    units(cellNum).anti.neural.base.ts_pst_win = units(cellNum).anti.neural.sacc.ts_pst(base_win);
    units(cellNum).anti.neural.base.rate_pst_win = units(cellNum).anti.neural.sacc.rate_pst(base_win);
    units(cellNum).anti.neural.base.rate_mu = mean(base_spks_anti);
    units(cellNum).anti.neural.base.rate_sig = std(base_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.base.rate_std = std(base_spks_anti);
    
    
     % sacc after 0 > or < than baseline
    if mean(sacc_on_spks_anti) > mean(sacc_pre_spks_anti)
    units(cellNum).anti.neural.exc = 1; units(cellNum).anti.neural.sup = 0; 
    else
         units(cellNum).anti.neural.exc = 0; units(cellNum).anti.neural.sup = 1; 
    end
    
    %% compare windows against baseline activity for pro and anti - nspk
    % pro
    [units(cellNum).stats.pro.pval.instrVSbase_nspk, units(cellNum).stats.pro.flags.instrVSbase_nspk] = signrank(units(cellNum).pro.neural.instr.nspk,units(cellNum).pro.neural.base.nspk);
    [units(cellNum).stats.pro.pval.saccVSbase_nspk, units(cellNum).stats.pro.flags.saccVSbase_nspk] = signrank(units(cellNum).pro.neural.sacc.nspk,units(cellNum).pro.neural.base.nspk);
    
    % anti
    [units(cellNum).stats.anti.pval.instrVSbase_nspk, units(cellNum).stats.anti.flags.instrVSbase_nspk] = signrank(units(cellNum).anti.neural.instr.nspk,units(cellNum).anti.neural.base.nspk);
    [units(cellNum).stats.anti.pval.saccVSbase_nspk, units(cellNum).stats.anti.flags.saccVSbase_nspk] = signrank(units(cellNum).anti.neural.sacc.nspk,units(cellNum).anti.neural.base.nspk);
    
    %     %% compare windows against baseline activity for pro and anti - rate
    %     <--- stats should be done on nspks
    %     % pro
    %     [units(cellNum).stats.pro.flags.instrVSbase, units(cellNum).stats.pro.pval.instrVSbase] = ttest(units(cellNum).pro.neural.instr.rate_pst_win,units(cellNum).pro.neural.base.rate_pst_win);
    %     [units(cellNum).stats.pro.flags.saccVSbase, units(cellNum).stats.pro.pval.saccVSbase] = ttest(units(cellNum).pro.neural.sacc.rate_pst_win,units(cellNum).pro.neural.base.rate_pst_win);
    %
    %     % anti
    %     [units(cellNum).stats.anti.flags.instrVSbase, units(cellNum).stats.anti.pval.instrVSbase] = ttest(units(cellNum).anti.neural.instr.rate_pst_win,units(cellNum).anti.neural.base.rate_pst_win);
    %     [units(cellNum).stats.anti.flags.saccVSbase, units(cellNum).stats.anti.pval.saccVSbase] = ttest(units(cellNum).anti.neural.sacc.rate_pst_win,units(cellNum).anti.neural.base.rate_pst_win);
    %
    
    %% compare pro vs anti - nspk
    %paired t-test
    [units(cellNum).stats.instr.flags.proVsAnti_instr, units(cellNum).stats.instr.pval.proVsAnti_instr] = ttest2(units(cellNum).pro.neural.instr.nspk,units(cellNum).anti.neural.instr.nspk);
    [units(cellNum).stats.sacc.flags.proVsAnti_sacc, units(cellNum).stats.sacc.pval.proVsAnti_sacc] = ttest2(units(cellNum).pro.neural.sacc.nspk,units(cellNum).anti.neural.sacc.nspk);
    
    %kstest2
    [units(cellNum).stats.instr.flags.proVsAnti_instr_ks_nspk, units(cellNum).stats.instr.pval.proVsAnti_instr_ks_nspk] = kstest2(units(cellNum).pro.neural.instr.nspk , units(cellNum).anti.neural.instr.nspk);
    [units(cellNum).stats.sacc.flags.proVsAnti_sacc_ks_nspk, units(cellNum).stats.sacc.pval.proVsAnti_sacc_ks_nspk] = kstest2(units(cellNum).pro.neural.sacc.nspk , units(cellNum).anti.neural.sacc.nspk);
    
    [units(cellNum).stats.instr.flags.proVsAnti_instr_ks_t_spk, units(cellNum).stats.instr.pval.proVsAnti_instr_ks_t_spk] = kstest2(units(cellNum).pro.neural.instr.tspk_all , units(cellNum).anti.neural.instr.tspk_all);
    [units(cellNum).stats.sacc.flags.proVsAnti_sacc_ks_t_spk, units(cellNum).stats.sacc.pval.proVsAnti_sacc_ks_t_spk] = kstest2(units(cellNum).pro.neural.sacc.tspk_all , units(cellNum).anti.neural.sacc.tspk_all);
    
    
    %% Compute change in FR from meam and baseline
    % pro
    units(cellNum).pro.neural.instr.delta_rate = units(cellNum).pro.neural.instr.rate_pst_win - units(cellNum).pro.neural.instr.rate_mu;
    units(cellNum).pro.neural.sacc.delta_rate = units(cellNum).pro.neural.sacc.rate_pst_win - units(cellNum).pro.neural.sacc.rate_mu;
    
    units(cellNum).pro.neural.sacc.delta_rate_base = units(cellNum).pro.neural.sacc.rate_pst_win - units(cellNum).pro.neural.base.rate_mu;
    units(cellNum).pro.neural.instr.delta_rate_base = units(cellNum).pro.neural.instr.rate_pst_win - units(cellNum).pro.neural.base.rate_mu;
    
    % anti
    units(cellNum).anti.neural.instr.delta_rate = units(cellNum).anti.neural.instr.rate_pst_win - units(cellNum).anti.neural.instr.rate_mu;
    units(cellNum).anti.neural.sacc.delta_rate = units(cellNum).anti.neural.sacc.rate_pst_win - units(cellNum).anti.neural.sacc.rate_mu;
    
    units(cellNum).anti.neural.sacc.delta_rate_base = units(cellNum).anti.neural.sacc.rate_pst_win - units(cellNum).anti.neural.base.rate_mu;
    units(cellNum).anti.neural.instr.delta_rate_base = units(cellNum).anti.neural.instr.rate_pst_win - units(cellNum).anti.neural.base.rate_mu;
    
    
    %% compute DDI
    % pro instr
    trialNum_pro = size(units(cellNum).pro.neural.sacc.nspk,1);
    SSE_pro_instr = sum((units(cellNum).pro.neural.instr.nspk-mean(units(cellNum).pro.neural.instr.nspk)).^2);
    units(cellNum).stats.pro.instr.DDI = (max(units(cellNum).pro.neural.sacc.nspk)-min(units(cellNum).pro.neural.sacc.nspk))/((max(units(cellNum).pro.neural.sacc.nspk)-min(units(cellNum).pro.neural.sacc.nspk))+2*sqrt(SSE_pro_instr/(trialNum_pro-2)));
    % pro sacc
    SSE_pro = sum((units(cellNum).pro.neural.sacc.nspk-mean(units(cellNum).pro.neural.sacc.nspk)).^2);
    units(cellNum).stats.pro.sacc.DDI = (max(units(cellNum).pro.neural.sacc.nspk)-min(units(cellNum).pro.neural.sacc.nspk))/((max(units(cellNum).pro.neural.sacc.nspk)-min(units(cellNum).pro.neural.sacc.nspk))+2*sqrt(SSE_pro/(trialNum_pro-2)));
    
    % anti
    trialNum_anti = size(units(cellNum).anti.neural.sacc.nspk,1);
    SSE_anti_instr = sum((units(cellNum).anti.neural.instr.nspk-mean(units(cellNum).anti.neural.instr.nspk)).^2);
    units(cellNum).stats.anti.instr.DDI = (max(units(cellNum).anti.neural.sacc.nspk)-min(units(cellNum).anti.neural.sacc.nspk))/((max(units(cellNum).anti.neural.sacc.nspk)-min(units(cellNum).anti.neural.sacc.nspk))+2*sqrt(SSE_anti_instr/(trialNum_anti-2)));
    % pro sacc
    SSE_anti = sum((units(cellNum).anti.neural.sacc.nspk-mean(units(cellNum).anti.neural.sacc.nspk)).^2);
    units(cellNum).stats.anti.sacc.DDI = (max(units(cellNum).anti.neural.sacc.nspk)-min(units(cellNum).anti.neural.sacc.nspk))/((max(units(cellNum).anti.neural.sacc.nspk)-min(units(cellNum).anti.neural.sacc.nspk))+2*sqrt(SSE_anti/(trialNum_anti-2)));
    
    
    %% normalized FR - z-scored aligned to sacc
    %% pro --  aligned to saccade
    units(cellNum).pro.neural.sacc.norm.rate_pst = (units(cellNum).pro.neural.sacc.rate_pst - mean(units(cellNum).pro.neural.sacc.rate_pst))/...
        std(units(cellNum).pro.neural.sacc.rate_pst);
    %% anti --  aligned to saccade
    units(cellNum).anti.neural.sacc.norm.rate_pst = (units(cellNum).anti.neural.sacc.rate_pst - mean(units(cellNum).anti.neural.sacc.rate_pst))/...
        std(units(cellNum).anti.neural.sacc.rate_pst);
    %% normalize anti based on prosaccades (as described by Maarten, March 2019)
    units(cellNum).anti.neural.sacc.norm.rate_pst_pro = (units(cellNum).anti.neural.sacc.rate_pst - mean(units(cellNum).pro.neural.sacc.rate_pst))/...
        std(units(cellNum).pro.neural.sacc.rate_pst);
    % grab normalized baseline vals
    norm_pro_base = mean(units(cellNum).pro.neural.sacc.norm.rate_pst(base_win));
    % compute change in firing rate for pro and anti with normalized 
    units(cellNum).pro.neural.sacc.norm.delta_rate = units(cellNum).pro.neural.sacc.norm.rate_pst(sacc_win) - norm_pro_base;
    units(cellNum).anti.neural.sacc.norm.delta_rate = units(cellNum).anti.neural.sacc.norm.rate_pst_pro(sacc_win) - norm_pro_base; 
    
    %% compute change index in norm -- CHECK! 
    units(cellNum).stats.instr.change_indx = (units(cellNum).pro.neural.instr.rate_mu  - units(cellNum).anti.neural.instr.rate_mu)/(units(cellNum).pro.neural.instr.rate_mu  + units(cellNum).anti.neural.instr.rate_mu);
    units(cellNum).stats.sacc.change_indx = (units(cellNum).pro.neural.sacc.rate_mu  - units(cellNum).anti.neural.sacc.rate_mu)/(units(cellNum).pro.neural.sacc.rate_mu  + units(cellNum).anti.neural.sacc.rate_mu);
    
    %% statistical test to compare if pro == anti - aligned to instr using spk count
    % H0:pro=anti  versus HA:pro?anti
    correctProTrials = units(cellNum).pro.indx_correctProTrials; % index pro trials
    correctAntiTrials = units(cellNum).anti.indx_correctAntiTrials;
    ntrls_pro = length(correctProTrials);
    ntrls_anti = length(correctAntiTrials);
    
    % instr
    for i=1:length(units(cellNum).pro.neural.instr.spkCount(1,:))
        p=(ntrls_pro*units(cellNum).pro.neural.instr.pbDist(i)+ntrls_anti*units(cellNum).anti.neural.instr.pbDist(i))/(ntrls_pro+ntrls_anti);
        units(cellNum).stats.instr.pval.pbDist_testStat(i) = (units(cellNum).pro.neural.instr.pbDist(1,i)-units(cellNum).anti.neural.instr.pbDist(1,i))/(sqrt(p*(1-p)*(1/ntrls_pro + 1/ntrls_anti)));
        if units(cellNum).stats.instr.pval.pbDist_testStat(i) > prs.signif_criteria | units(cellNum).stats.instr.pval.pbDist_testStat(i) < -prs.signif_criteria
            units(cellNum).stats.instr.flags.pbDist(i) = 1;
        else
            units(cellNum).stats.instr.flags.pbDist(i) = 0;
        end
    end
    
    % sacc
    for i=1:length(units(cellNum).anti.neural.sacc.spkCount(1,:))
        p=(ntrls_pro*units(cellNum).pro.neural.sacc.pbDist(i)+ntrls_anti*units(cellNum).anti.neural.sacc.pbDist(i))/(ntrls_pro+ntrls_anti);
        units(cellNum).stats.sacc.pval.pbDist_testStat(i) = (units(cellNum).pro.neural.sacc.pbDist(1,i)-units(cellNum).anti.neural.sacc.pbDist(1,i))/(sqrt(p*(1-p)*(1/ntrls_pro + 1/ntrls_anti)));
        if units(cellNum).stats.sacc.pval.pbDist_testStat(i) > prs.signif_criteria | units(cellNum).stats.sacc.pval.pbDist_testStat(i) < -prs.signif_criteria
            units(cellNum).stats.sacc.flags.pbDist(i) = 1;
        else
            units(cellNum).stats.sacc.flags.pbDist(i) = 0;
        end
    end
    
    %% statistical test to compare if pro == anti - aligned to instr using spk count with bigger window
    
    for win_num = 1:length(units(cellNum).pro.neural.instr.spkCount_win(1,:));
        % instr pro vs anti
        [units(cellNum).stats.instr.flags.spk_count_bigWin(win_num),units(cellNum).stats.instr.pval.spk_count_bigWin(win_num)] = ttest2(units(cellNum).pro.neural.instr.spkCount_win(:,win_num),...
            units(cellNum).anti.neural.instr.spkCount_win(:,win_num));
        
        % sacc pro vs anti
        [units(cellNum).stats.sacc.flags.spk_count_bigWin(win_num),units(cellNum).stats.sacc.pval.spk_count_bigWin(win_num)] = ttest2(units(cellNum).pro.neural.sacc.spkCount_win(:,win_num),...
            units(cellNum).anti.neural.sacc.spkCount_win(:,win_num));
    end
    
     %% Sliding window to test if neuron is exc or supp after saccade onset
    % detect saccade related using sliding window (5 bins) 200 ms before
    % sacc onset to 200 ms after sacc onset and store time at > or < than 2*std
    
%         indx_win = prs.slide_win_size; win_size_prev = prs.slide_win_size_prev;
%         t = units(cellNum).pro.neural.sacc.ts_pst;
%         indx_beg = find(t>0,1); indx_end = find(t>0.150,1);
%     %
%     %     % pro
%         rate_pst = units(cellNum).pro.neural.sacc.rate_pst;
%         exc_pro=0; sup_pro=0;exc_anti=0; sup_anti=0;
%     
%         for indx = indx_beg:indx_end
%             r_mu = mean(rate_pst(round(indx-win_size/2):round(indx+win_size/2)));
%             r_thresh1 = mean(rate_pst(indx-win_size_prev:indx))+...
%                 2*std(rate_pst(indx-win_size_prev:indx));
%             r_thresh2 = mean(rate_pst(indx-win_size_prev:indx))-...
%                 2*std(rate_pst(indx-win_size_prev:indx));
%             if (r_thresh1>1 && r_mu>r_thresh1), exc_pro=exc_pro+1;else exc_pro=0;end
%             if (r_thresh2>1 && r_mu<r_thresh2), sup_pro=sup_pro+1;else sup_pro=0;end
%             if exc_pro==3, events.rise.t_on = t(indx); units(cellNum).pro.events.rise.type_on = 'exc';
%                 break;
%             elseif sup_pro==3, units(cellNum).pro.events.rise.t_on = t(indx); units(cellNum).pro.events.rise.type_on = 'sup';
%                 break;
%             end
%         end
    %
    %     % anti
    %     rate_pst = units(cellNum).anti.neural.sacc.rate_pst;
    %     exc_anti=0; sup_anti=0;
    %
    %     for indx = indx_beg:indx_end
    %         r_mu = mean(rate_pst(round(indx-win_size/2):round(indx+win_size/2)));
    %         r_thresh1 = mean(rate_pst(indx-win_size_prev:indx))+...
    %             2*std(rate_pst(indx-win_size_prev:indx));
    %         r_thresh2 = mean(rate_pst(indx-win_size_prev:indx))-...
    %             2*std(rate_pst(indx-win_size_prev:indx));
    %         if (r_thresh1>1 && r_mu>r_thresh1), exc_anti=exc_anti+1;else exc_anti=0;end
    %         if (r_thresh2>1 && r_mu<r_thresh2), sup_anti=sup_anti+1;else sup_anti=0;end
    %         if exc_anti==3, events.rise.t_on = t(indx); events.rise.type_on = 'exc';
    %             break;
    %         elseif sup_anti==3, units(cellNum).pro.events.rise.t_on = t(indx); units(cellNum).pro.events.rise.type_on = 'sup';
    %             break;
    %         end
    %     end
    
    %% Sliding window to test diff pro vs anti
    % detect saccade related using sliding window (5 bins) 200 ms before
    % sacc onset to 200 ms after sacc onset and store time at > or < than 2*std
    
    %     indx_win = prs.slide_win_size; win_size_prev = prs.slide_win_size_prev;
    %     t = units(cellNum).pro.neural.sacc.ts_pst;
    %     indx_beg = find(t>-0.2,1); indx_end = find(t>0.2,1);
    %
    %     % pro
    %     rate_pst = units(cellNum).pro.neural.sacc.rate_pst;
    %     exc_pro=0; sup_pro=0;exc_anti=0; sup_anti=0;
    %
    %     for indx = indx_beg:indx_end
    %         r_mu = mean(rate_pst(round(indx-win_size/2):round(indx+win_size/2)));
    %         r_thresh1 = mean(rate_pst(indx-win_size_prev:indx))+...
    %             2*std(rate_pst(indx-win_size_prev:indx));
    %         r_thresh2 = mean(rate_pst(indx-win_size_prev:indx))-...
    %             2*std(rate_pst(indx-win_size_prev:indx));
    %         if (r_thresh1>1 && r_mu>r_thresh1), exc_pro=exc_pro+1;else exc_pro=0;end
    %         if (r_thresh2>1 && r_mu<r_thresh2), sup_pro=sup_pro+1;else sup_pro=0;end
    %         if exc_pro==3, events.rise.t_on = t(indx); events.rise.type_on = 'exc';
    %             break;
    %         elseif sup_pro==3, units(cellNum).pro.events.rise.t_on = t(indx); units(cellNum).pro.events.rise.type_on = 'sup';
    %             break;
    %         end
    %     end
    %
    %     % anti
    %     rate_pst = units(cellNum).anti.neural.sacc.rate_pst;
    %     exc_anti=0; sup_anti=0;
    %
    %     for indx = indx_beg:indx_end
    %         r_mu = mean(rate_pst(round(indx-win_size/2):round(indx+win_size/2)));
    %         r_thresh1 = mean(rate_pst(indx-win_size_prev:indx))+...
    %             2*std(rate_pst(indx-win_size_prev:indx));
    %         r_thresh2 = mean(rate_pst(indx-win_size_prev:indx))-...
    %             2*std(rate_pst(indx-win_size_prev:indx));
    %         if (r_thresh1>1 && r_mu>r_thresh1), exc_anti=exc_anti+1;else exc_anti=0;end
    %         if (r_thresh2>1 && r_mu<r_thresh2), sup_anti=sup_anti+1;else sup_anti=0;end
    %         if exc_anti==3, events.rise.t_on = t(indx); events.rise.type_on = 'exc';
    %             break;
    %         elseif sup_anti==3, units(cellNum).pro.events.rise.t_on = t(indx); units(cellNum).pro.events.rise.type_on = 'sup';
    %             break;
    %         end
    %     end
    
end

%% coefficient of variation and CV2
% CV_isi = mean_isi/std_isi
% CV2 =  2x|ISIn+1-ISIn|/(ISIn+1 + ISIn) Holt et al. 199

% CV= std dev/mean


for cellNum = 1:length(units)
    correctProTrials = units(cellNum).pro.indx_correctProTrials; % index pro trials
    correctAntiTrials = units(cellNum).anti.indx_correctAntiTrials;
    % pro
    for trialNum = 1:length(correctProTrials)
        isi_tspk = diff(units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc); std_isi = std(isi_tspk); mean_isi = mean(isi_tspk);
        % compute CV isi
        units(cellNum).pro.neural.trial(trialNum).cv_isi= mean_isi/std_isi;
        %compute CV2
        for i=2:length(isi_tspk)-1
            cv2_trialNum(i) = 2*((abs(isi_tspk(i+1)-isi_tspk(i))))/(isi_tspk(i+1)+isi_tspk(i));
            units(cellNum).pro.neural.trial(trialNum).cv2 = mean(cv2_trialNum);
        end
        
    end
    units(cellNum).pro.meanCV_isi =  mean([units(cellNum).pro.neural.trial.cv_isi]);
    
    %anti
    for trialNum = 1:length(correctAntiTrials)
        isi_tspk = diff(units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc); std_isi = std(isi_tspk); mean_isi = mean(isi_tspk);
        % compute CV isi
        units(cellNum).anti.neural.trial(trialNum).cv_isi= mean_isi/std_isi;
        %compute CV2
        for i=2:length(isi_tspk)-1
            cv2_trialNum(i) = 2*((abs(isi_tspk(i+1)-isi_tspk(i))))/(isi_tspk(i+1)+isi_tspk(i));
            units(cellNum).anti.neural.trial(trialNum).cv2 = mean(cv2_trialNum);
        end
        
    end
    
    units(cellNum).anti.meanCV_isi =  mean([units(cellNum).anti.neural.trial.cv_isi]);
    
end

%% Modulation ratio
for cellNum = 1:length(units)
    units(cellNum).stats.mod_ratio = units(cellNum).pro.neural.sacc.rate_pst_win / units(cellNum).anti.neural.sacc.rate_pst_win;
end


end


