function [units,pop] = extractWholeNeuronR(wholeNeuronResults)
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

% first find the ones that are not empty and have at least 10 trials. 
for i = 1:length(wholeNeuronResults)
    if  ~isempty(wholeNeuronResults(i).allStableTrials)
        cell_indx(i) = numel(wholeNeuronResults(i).selectedTrials.corProTrials)>=prs.min_trial & numel(wholeNeuronResults(i).selectedTrials.corAntiTrials)>=prs.min_trial;
    end
end
cell_indx = find(cell_indx);


for cellNum = 1:length(cell_indx)
    trials_with_spk=zeros(1,length(wholeNeuronResults(cell_indx(cellNum)).allStableTrials));
    
    for trialNum = 1:length(wholeNeuronResults(cell_indx(cellNum)).allStableTrials)
        
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
                units(cellNum).trial.neural(trialNum).tspk_SS_align_instrDir = units(cellNum).trial.neural(trialNum).tspk_SS - (units(cellNum).trial.behav(trialNum).goCueTime-0.1); 
                units(cellNum).trial.neural(trialNum).ts_eye_sacc = units(cellNum).trial.behav(trialNum).eye_ts-units(cellNum).trial.behav(trialNum).saccadeOnset;
                
                
            elseif run_error_trials
                
                indx_trial_spks = spks>prs.tspk(1) & spks < units(cellNum).trial.behav(trialNum).goCueTime+0.5+prs.tspk(2); % just pick 100 ms before trial starts to reward +200 ms
                units(cellNum).trial.neural(trialNum).tspk_SS = spks(indx_trial_spks);
                units(cellNum).trial.neural(trialNum).tspk_SS_align_sacc =  units(cellNum).trial.neural(trialNum).tspk_SS-units(cellNum).trial.behav(trialNum).saccadeOnset; % contains spike times for CS aligned to sacc onset
                units(cellNum).trial.neural(trialNum).tspk_SS_align_instrDir = units(cellNum).trial.neural(trialNum).tspk_SS - (units(cellNum).trial.behav(trialNum).goCueTime-0.1);
            else
                %keyboard
                units(cellNum).trial.neural(trialNum).tspk_SS = [];
                units(cellNum).trial.neural(trialNum).tspk_SS_align_sacc = [];
                units(cellNum).trial.neural(trialNum).tspk_SS_align_instrDir = [];
            end
            units(cellNum).id = 'SS';
            % select condition type Pro or Antisaccade
            if ismember(wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).conditionCode, prs.proConditions)
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
analyse_sacc_align = 0; analyse_instrDir_align = 0; analyse_instr_back = 0;

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
        [units(cellNum).pro.neural.instr.rate_pst,units(cellNum).pro.neural.instr.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,prs.timepoints_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to trial onset
    catch
        keyboard
    end
    
    units(cellNum).pro.neural.instr.rate_pst = smooth_pst(units(cellNum).pro.neural.instr.rate_pst,prs.binwidth,prs.tsmooth);
    
    % sacc aligned
    analyse_sacc_align=1;
    [units(cellNum).pro.neural.sacc.rate_pst,units(cellNum).pro.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to saccade onset
    units(cellNum).pro.neural.sacc.rate_pst = smooth_pst(units(cellNum).pro.neural.sacc.rate_pst,prs.binwidth,prs.tsmooth);
    analyse_sacc_align=0;
    
    % instr dir aligned
    analyse_instrDir_align = 1; 
    [units(cellNum).pro.neural.instrDir.rate_pst,units(cellNum).pro.neural.instrDir.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,prs.timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to saccade onset
    units(cellNum).pro.neural.instrDir.rate_pst = smooth_pst(units(cellNum).pro.neural.instrDir.rate_pst,prs.binwidth,prs.tsmooth);
    analyse_instrDir_align = 0; 
    %% sanity check -- for every cell, randomize half of the trials and compute psth (5 times)
    analyse_sacc_align=1;
    for j = 1:100
    numtrials = []; pick_trials = []; 
    numtrials = length(units(cellNum).pro.neural.trial); % number of trials
    pick_trials = sort(randsample(numtrials, round(numtrials/2)));  % pick half of the trials randomly
    other_trials = 1:numtrials; other_half = other_trials(~ismember(other_trials,pick_trials)); 
    % units(cellNum).pro.neural.sacc.rate_pst_computed_half_indx(j,:) = pick_trials; 
    % units(cellNum).pro.neural.sacc.rate_pst_other_half_indx(j,:) = other_half;
    
    [units(cellNum).pro.neural.sacc.rate_pst_rand(j,:),~] = Spiketimes2Rate(units(cellNum).pro.neural.trial((pick_trials)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to saccade onset
    units(cellNum).pro.neural.sacc.rate_pst_rand(j,:) = smooth_pst(units(cellNum).pro.neural.sacc.rate_pst_rand(j,:),prs.binwidth,prs.tsmooth);
    
    [units(cellNum).pro.neural.sacc.rate_pst_rand_other_half(j,:),~] = Spiketimes2Rate(units(cellNum).pro.neural.trial((other_half)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to saccade onset
    units(cellNum).pro.neural.sacc.rate_pst_rand_other_half(j,:) = smooth_pst(units(cellNum).pro.neural.sacc.rate_pst_rand(j,:),prs.binwidth,prs.tsmooth); 
    end
    analyse_sacc_align=0;
%%  
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
    [units(cellNum).anti.neural.instr.rate_pst,units(cellNum).anti.neural.instr.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,prs.timepoints_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to trial onset
    units(cellNum).anti.neural.instr.rate_pst = smooth_pst(units(cellNum).anti.neural.instr.rate_pst,prs.binwidth,prs.tsmooth);
    % sacc aligned
    analyse_sacc_align=1;
    [units(cellNum).anti.neural.sacc.rate_pst,units(cellNum).anti.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to saccade onset
    units(cellNum).anti.neural.sacc.rate_pst = smooth_pst(units(cellNum).anti.neural.sacc.rate_pst,prs.binwidth,prs.tsmooth);
    
    % instr dir aligned
    analyse_instrDir_align = 1; 
    [units(cellNum).anti.neural.instrDir.rate_pst,units(cellNum).anti.neural.instrDir.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,prs.timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to saccade onset
    units(cellNum).anti.neural.instrDir.rate_pst = smooth_pst(units(cellNum).anti.neural.instrDir.rate_pst,prs.binwidth,prs.tsmooth);
    analyse_instrDir_align = 0; 
    
    %% sanity check -- for every cell, pick random half of the trials and compute psth
    for j = 1:100;
        numtrials = []; pick_trials = [];
        numtrials = length(units(cellNum).anti.neural.trial); % number of trials
        pick_trials = sort(randsample(numtrials, round(numtrials/2)));  % pick half of the trials randomly
        other_trials = 1:numtrials; other_half = other_trials(~ismember(other_trials,pick_trials));
        % units(cellNum).anti.neural.sacc.rate_pst_computed_half_indx(j,:) = pick_trials;
        % units(cellNum).anti.neural.sacc.rate_pst_other_half_indx(j,:) = other_half;
        
        [units(cellNum).anti.neural.sacc.rate_pst_rand(j,:),~] = Spiketimes2Rate(units(cellNum).anti.neural.trial((pick_trials)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to saccade onset
        units(cellNum).anti.neural.sacc.rate_pst_rand(j,:) = smooth_pst(units(cellNum).anti.neural.sacc.rate_pst_rand(j,:),prs.binwidth,prs.tsmooth);
        
        [units(cellNum).anti.neural.sacc.rate_pst_rand_other_half(j,:),~] = Spiketimes2Rate(units(cellNum).anti.neural.trial((other_half)),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % aligned to saccade onset
        units(cellNum).anti.neural.sacc.rate_pst_rand_other_half(j,:) = smooth_pst(units(cellNum).anti.neural.sacc.rate_pst_rand(j,:),prs.binwidth,prs.tsmooth);
    end
    analyse_sacc_align=0;
end

%% For every trial - Spike time count and spike times to rate for pro and antisaccades

for cellNum = 1:length(units)
    id=units(cellNum).id;
    %pro trials
    correctProTrials = units(cellNum).pro.indx_correctProTrials; % index pro trials
    % get spk counts
    for trialNum = 1:length(correctProTrials)
        
        [units(cellNum).pro.neural.instr.spkCount(trialNum,:),units(cellNum).pro.neural.instr.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).pro.neural.trial(trialNum),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % spk counts
        analyse_sacc_align=1;
        [units(cellNum).pro.neural.sacc.spkCount(trialNum,:),units(cellNum).pro.neural.sacc.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).pro.neural.trial(trialNum),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_sacc_align=0;
        analyse_instrDir_align=1;
        [units(cellNum).pro.neural.instrDir.spkCount(trialNum,:),units(cellNum).pro.neural.instrDir.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).pro.neural.trial(trialNum),prs.timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        [units(cellNum).pro.neural.instr_back.spkCount(trialNum,:),units(cellNum).pro.neural.instr_back.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).pro.neural.trial(trialNum), prs.timepoints_instrDir, prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); 
        analyse_instrDir_align=0;
    end
    % Take out probability of spk in pro
    for i=1:length(units(cellNum).pro.neural.instr.spkCount(1,:))
        % instr aligned
        units(cellNum).pro.neural.instr.pbDist(1,i)= sum(units(cellNum).pro.neural.instr.spkCount(:,i))/length(correctProTrials); % compute probability of spk instr onset
        units(cellNum).pro.neural.instr.pbDist_sem = sqrt(units(cellNum).pro.neural.instr.pbDist(1,i)*(1-units(cellNum).pro.neural.instr.pbDist(1,i))/length(correctProTrials));
    end
    
    for i=1:length(units(cellNum).pro.neural.instr_back.spkCount(1,:))
        % instr aligned
        units(cellNum).pro.neural.instr_back.pbDist(1,i)= sum(units(cellNum).pro.neural.instr_back.spkCount(:,i))/length(correctProTrials); % compute probability of spk instr onset
        units(cellNum).pro.neural.instr_back.pbDist_sem = sqrt(units(cellNum).pro.neural.instr_back.pbDist(1,i)*(1-units(cellNum).pro.neural.instr_back.pbDist(1,i))/length(correctProTrials));
    end
    
    for i=1:length(units(cellNum).pro.neural.sacc.spkCount(1,:))
        % sacc aligned
        units(cellNum).pro.neural.sacc.pbDist(1,i)= sum(units(cellNum).pro.neural.sacc.spkCount(:,i))/length(correctProTrials); % compute probability of spk sacc aligned
        units(cellNum).pro.neural.sacc.pbDist_sem = sqrt(units(cellNum).pro.neural.sacc.pbDist(1,i)*(1-units(cellNum).pro.neural.sacc.pbDist(1,i))/length(correctProTrials));
    end
%     for i = 1:length(units(cellNum).pro.neural.instrDir.spkCount(1,:))
%         units(cellNum).pro.neural.instrDir.pbDist(1,i)= sum(units(cellNum).pro.neural.instrDir.spkCount(:,i))/length(correctProTrials); % compute probability of spk instr onset
%         units(cellNum).pro.neural.instrDir.pbDist_sem = sqrt(units(cellNum).pro.neural.instrDir.pbDist(1,i)*(1-units(cellNum).pro.neural.instrDir.pbDist(1,i))/length(correctProTrials));
%     
%         units(cellNum).pro.neural.instr_back.pbDist(1,i)= sum(units(cellNum).pro.neural.instr_back.spkCount(:,i))/length(correctProTrials); % compute probability of spk instr onset
%         units(cellNum).pro.neural.instr_back.pbDist_sem = sqrt(units(cellNum).pro.neural.instr_back.pbDist(1,i)*(1-units(cellNum).pro.neural.instr_back.pbDist(1,i))/length(correctProTrials));
%     
%     end
    
    % anti trials
    correctAntiTrials = units(cellNum).anti.indx_correctAntiTrials;
    for trialNum = 1:length(correctAntiTrials)
        [units(cellNum).anti.neural.instr.spkCount(trialNum,:),units(cellNum).anti.neural.instr.trial(trialNum).ts] = Spiketimes2CountTrial(units(cellNum).anti.neural.trial(trialNum),prs.timepoints_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); % spk counts
        analyse_sacc_align=1;
        [units(cellNum).anti.neural.sacc.spkCount(trialNum,:),units(cellNum).anti.neural.sacc.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).anti.neural.trial(trialNum),prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_sacc_align=0;
        analyse_instrDir_align=1;
        [units(cellNum).anti.neural.instrDir.spkCount(trialNum,:),units(cellNum).anti.neural.instrDir.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).anti.neural.trial(trialNum),prs.timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        [units(cellNum).anti.neural.instr_back.spkCount(trialNum,:),units(cellNum).anti.neural.instr_back.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).anti.neural.trial(trialNum), prs.timepoints_instrDir, prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id); 
        analyse_instrDir_align=0;
        
    end
    
    % Take out probability of spk in anti
    for i=1:length(units(cellNum).anti.neural.instr.spkCount(1,:))
        % instr aligned
        units(cellNum).anti.neural.instr.pbDist(1,i)= sum(units(cellNum).anti.neural.instr.spkCount(:,i))/length(correctAntiTrials); % compute probability of spk instr onset
        units(cellNum).anti.neural.instr.pbDist_sem = sqrt(units(cellNum).pro.neural.instr.pbDist(1,i)*(1-units(cellNum).anti.neural.instr.pbDist(1,i))/length(correctAntiTrials));
    end
    
     for i=1:length(units(cellNum).anti.neural.instr_back.spkCount(1,:))
        % instr aligned
        units(cellNum).anti.neural.instr_back.pbDist(1,i)= sum(units(cellNum).anti.neural.instr_back.spkCount(:,i))/length(correctAntiTrials); % compute probability of spk instr onset
        units(cellNum).anti.neural.instr_back.pbDist_sem = sqrt(units(cellNum).anti.neural.instr_back.pbDist(1,i)*(1-units(cellNum).anti.neural.instr_back.pbDist(1,i))/length(correctAntiTrials));
    end
        
     for i=1:length(units(cellNum).anti.neural.sacc.spkCount(1,:))
        % sacc aligned
        units(cellNum).anti.neural.sacc.pbDist(1,i)= sum(units(cellNum).anti.neural.sacc.spkCount(:,i))/length(correctAntiTrials); % compute probability of spk sacc aligned
        units(cellNum).anti.neural.sacc.pbDist_sem = sqrt(units(cellNum).anti.neural.sacc.pbDist(1,i)*(1-units(cellNum).anti.neural.sacc.pbDist(1,i))/length(correctAntiTrials));
    end
    
    % nspk pro - all bins and windows per trial
    units(cellNum).pro.neural.instr.tspk_all = []; units(cellNum).pro.neural.sacc.tspk_all = []; units(cellNum).pro.neural.instrDir.tspk_all = []; 
    for trialNum = 1:length(correctProTrials)
        
        % timepoints_instr = -1:prs.binwidth:1; timepoints_sacc = -1:prs.binwidth:1;
        
        %% align to instr
        [~,units(cellNum).pro.neural.instr.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum).tspk_SS,prs.timepoints_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % win
        tspk_instr= units(cellNum).pro.neural.trial(trialNum).tspk_SS>prs.instruction_win(1) &units(cellNum).pro.neural.trial(trialNum).tspk_SS<prs.instruction_win(2);
        t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
        units(cellNum).pro.neural.instr.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS(tspk_instr);
        [units(cellNum).pro.neural.instr.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.instr.tspk(trialNum,:)),timepoints_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        %% align to sacc
        analyse_sacc_align = 1;
        [~,units(cellNum).pro.neural.sacc.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % win
        tspk_sacc= units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc>prs.saccade_win(1) & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc<prs.saccade_win(2);
        t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc<prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
        units(cellNum).pro.neural.sacc.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc(tspk_sacc);
        [units(cellNum).pro.neural.sacc.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.sacc.tspk(trialNum,:)),timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        % presacc (for comparing exc and sup
        tspk_presacc = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc>prs.presaccade_win(1) &units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc<prs.presaccade_win(2);
        t_presacc = prs.timepoints_sacc > prs.presaccade_win(1) & prs.timepoints_sacc < prs.presaccade_win(2); timepoints_presacc = prs.timepoints_sacc(t_presacc);
        units(cellNum).pro.neural.presacc.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc(tspk_presacc);
        [units(cellNum).pro.neural.presacc.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.presacc.tspk(trialNum,:)),timepoints_presacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_sacc_align = 0;
        
        %% align to instr dir 
        analyse_instrDir_align = 1;
        [~,units(cellNum).pro.neural.instr_dir.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir,prs.timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % win
        tspk_instrDir = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir>prs.instr_dir(1) & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir<prs.instr_dir(2);
        t_instrDir = prs.timepoints_instrDir > prs.instr_dir(1) & prs.timepoints_instrDir<prs.instr_dir(2); timepoints_instrDir = prs.timepoints_instrDir(t_instrDir);
        units(cellNum).pro.neural.instrDir.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir(tspk_instrDir);
        [units(cellNum).pro.neural.instrDir.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.instrDir.tspk(trialNum,:)),timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_instrDir_align = 0;
        
         %% align to instr back 
        analyse_instrDir_align = 1;
        [~,units(cellNum).pro.neural.instr_back.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir,prs.timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % win
        tspk_instr_back = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir>prs.instr_back(1) & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir < prs.instr_back(2);
        t_instr_back = prs.timepoints_instrDir > prs.instr_back(1) & prs.timepoints_instrDir<prs.instr_back(2); timepoints_instr_back = prs.timepoints_instrDir(t_instr_back);
        units(cellNum).pro.neural.instr_back.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir(tspk_instr_back);
        [units(cellNum).pro.neural.instr_back.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.instr_back.tspk(trialNum,:)),timepoints_instr_back,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_instrDir_align = 0;
        
        %% align to go cue
        analyse_instrDir_align = 1; 
        tspk_goCue = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir>prs.go_cue(1) & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir<prs.go_cue(2);
        t_goCue = prs.timepoints_instrDir > prs.go_cue(1) & prs.timepoints_instrDir<prs.go_cue(2); timepoints_goCue = prs.timepoints_instrDir(t_goCue);
        units(cellNum).pro.neural.goCue.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_instrDir(tspk_goCue);
        [units(cellNum).pro.neural.goCue.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.goCue.tspk(trialNum,:)),timepoints_goCue,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_instrDir_align=0;
        
        % concatenate tspk 
        units(cellNum).pro.neural.instr.tspk_all = [units(cellNum).pro.neural.instr.tspk_all ; units(cellNum).pro.neural.instr.tspk{trialNum,:}];
        units(cellNum).pro.neural.sacc.tspk_all = [units(cellNum).pro.neural.sacc.tspk_all ; units(cellNum).pro.neural.sacc.tspk{trialNum,:}];
        units(cellNum).pro.neural.instrDir.tspk_all = [units(cellNum).pro.neural.instrDir.tspk_all ; units(cellNum).pro.neural.instrDir.tspk{trialNum,:}];
        
        %% baseline
        tspk_base = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc>prs.baseline_win(1) & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc<prs.baseline_win(2);
        t_base = prs.timepoints_sacc >prs.baseline_win(1) & prs.timepoints_sacc<prs.baseline_win(2); timepoints_base = prs.timepoints_sacc(t_base);
        units(cellNum).pro.neural.base.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc(tspk_base);
        [units(cellNum).pro.neural.base.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.base.tspk(trialNum,:)),timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_sacc_align = 0;
        
        % baseline instr
        tspk_base_instr = units(cellNum).pro.neural.trial(trialNum).tspk_SS>prs.baseline_instr(1) & units(cellNum).pro.neural.trial(trialNum).tspk_SS<prs.baseline_instr(2);
        t_base_instr = prs.timepoints_instr > prs.baseline_instr(1) & prs.timepoints_instr<prs.baseline_instr(2); timepoints_base_instr = prs.timepoints_instr(t_base_instr);
        units(cellNum).pro.neural.base.tspk_instr{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS(tspk_base_instr);
        [units(cellNum).pro.neural.base.nspk_instr(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.base.tspk_instr(trialNum,:)),timepoints_base_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_sacc_align = 0;
        
    end
    
    % nspk anti - windows per trial
    units(cellNum).anti.neural.instr.tspk_all = []; units(cellNum).anti.neural.sacc.tspk_all = []; units(cellNum).anti.neural.instrDir.tspk_all = []; 
    for trialNum = 1:length(correctAntiTrials)
        %% align to instr
        [~,units(cellNum).anti.neural.instr.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum).tspk_SS,prs.timepoints_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % win
        tspk_instr= units(cellNum).anti.neural.trial(trialNum).tspk_SS>prs.instruction_win(1) &units(cellNum).anti.neural.trial(trialNum).tspk_SS<prs.instruction_win(2);
        t_instr = prs.timepoints_instr > prs.instruction_win(1) & prs.timepoints_instr<prs.instruction_win(2); timepoints_instr = prs.timepoints_instr(t_instr);
        units(cellNum).anti.neural.instr.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS(tspk_instr);
        [units(cellNum).anti.neural.instr.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.instr.tspk(trialNum,:)),timepoints_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
               
        %% align to sacc
        analyse_sacc_align = 1;
        [~,units(cellNum).anti.neural.sacc.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc,prs.timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % win
        tspk_sacc = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc>prs.saccade_win(1) &units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc<prs.saccade_win(2);
        t_sacc = prs.timepoints_sacc > prs.saccade_win(1) & prs.timepoints_sacc < prs.saccade_win(2); timepoints_sacc = prs.timepoints_sacc(t_sacc);
        units(cellNum).anti.neural.sacc.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc(tspk_sacc);
        [units(cellNum).anti.neural.sacc.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.sacc.tspk(trialNum,:)),timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % presacc (for comparing exc and sup
        tspk_presacc = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc>prs.presaccade_win(1) &units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc<prs.presaccade_win(2);
        t_presacc = prs.timepoints_sacc > prs.presaccade_win(1) & prs.timepoints_sacc < prs.presaccade_win(2); timepoints_presacc = prs.timepoints_sacc(t_presacc);
        units(cellNum).anti.neural.presacc.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc(tspk_presacc);
        [units(cellNum).anti.neural.presacc.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.presacc.tspk(trialNum,:)),timepoints_presacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_sacc_align = 0;
        
        %% align to instr dir 
        analyse_instrDir_align = 1;
        [~,units(cellNum).anti.neural.instrDir.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir,prs.timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % win
        tspk_instrDir = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir>prs.instr_dir(1) & units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir<prs.instr_dir(2);
        t_instrDir = prs.timepoints_instrDir > prs.instr_dir(1) & prs.timepoints_instrDir<prs.instr_dir(2); timepoints_instrDir = prs.timepoints_instrDir(t_instrDir);
        units(cellNum).anti.neural.instrDir.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir(tspk_instrDir);
        [units(cellNum).anti.neural.instrDir.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.instrDir.tspk(trialNum,:)),timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_instrDir_align = 0;
        
        %% align to instr back 
        analyse_instrDir_align = 1;
        [~,units(cellNum).anti.neural.instr_back.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir,prs.timepoints_instrDir,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        
        % win
        tspk_instr_back = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir>prs.instr_back(1) & units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir<prs.instr_back(2);
        t_instr_back = prs.timepoints_instrDir > prs.instr_back(1) & prs.timepoints_instrDir<prs.instr_back(2); timepoints_instr_back = prs.timepoints_instrDir(t_instr_back);
        units(cellNum).anti.neural.instr_back.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir(tspk_instr_back);
        [units(cellNum).anti.neural.instr_back.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.instr_back.tspk(trialNum,:)),timepoints_instr_back,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_instrDir_align = 0;
        
        %% align to go cue
        analyse_instrDir_align = 1; 
        tspk_goCue = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir>prs.go_cue(1) & units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir<prs.go_cue(2);
        t_goCue = prs.timepoints_instrDir > prs.go_cue(1) & prs.timepoints_instrDir<prs.go_cue(2); timepoints_goCue = prs.timepoints_instrDir(t_goCue);
        units(cellNum).anti.neural.goCue.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_instrDir(tspk_goCue);
        [units(cellNum).anti.neural.goCue.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.goCue.tspk(trialNum,:)),timepoints_goCue,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_instrDir_align=0;
        
        
        % concatenate tspk 
        units(cellNum).anti.neural.instr.tspk_all = [units(cellNum).anti.neural.instr.tspk_all ; units(cellNum).anti.neural.instr.tspk{trialNum,:}];
        units(cellNum).anti.neural.sacc.tspk_all = [units(cellNum).anti.neural.sacc.tspk_all ; units(cellNum).anti.neural.sacc.tspk{trialNum,:}];
        units(cellNum).anti.neural.instrDir.tspk_all = [units(cellNum).anti.neural.instrDir.tspk_all ; units(cellNum).anti.neural.instrDir.tspk{trialNum,:}]; 
        
        % baseline
        tspk_base = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc>prs.baseline_win(1) & units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc<prs.baseline_win(2);
        t_base = prs.timepoints_sacc >prs.baseline_win(1) & prs.timepoints_sacc<prs.baseline_win(2); timepoints_base = prs.timepoints_sacc(t_base);
        units(cellNum).anti.neural.base.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc(tspk_base);
        [units(cellNum).anti.neural.base.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.base.tspk(trialNum,:)),timepoints_sacc,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
        analyse_sacc_align = 0;
        
        % baseline instr
        tspk_base_instr = units(cellNum).anti.neural.trial(trialNum).tspk_SS>prs.baseline_instr(1) & units(cellNum).anti.neural.trial(trialNum).tspk_SS<prs.baseline_instr(2);
        t_base_instr = prs.timepoints_instr >prs.baseline_instr(1) & prs.timepoints_instr<prs.baseline_instr(2); timepoints_base_instr = prs.timepoints_instr(t_base_instr);
        units(cellNum).anti.neural.base.tspk_instr{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS(tspk_base_instr);
        [units(cellNum).anti.neural.base.nspk_instr(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.base.tspk_instr(trialNum,:)),timepoints_base_instr,prs.binwidth,analyse_sacc_align,analyse_instrDir_align,id);
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

%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%% Stats Eye %%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
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

%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%% Neural %%%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%

%% stats per neuron - compare windows against baseline activity for pro and anti (concatenate all trials and make comparison)
%% Extract spks in windows and compute mean and sem
for cellNum = 1:length(units)
    ntrls_pro = length(units(cellNum).pro.neural.trial);
    ntrls_anti = length(units(cellNum).anti.neural.trial);
    % windows
    instr_win = units(cellNum).pro.neural.instr.ts_pst > prs.instruction_win(1) & units(cellNum).pro.neural.instr.ts_pst <= prs.instruction_win(2);
    sacc_win = units(cellNum).pro.neural.sacc.ts_pst >= prs.saccade_win(1) & units(cellNum).pro.neural.sacc.ts_pst <= prs.saccade_win(2);
    instr_back_win = units(cellNum).pro.neural.instrDir.ts_pst > prs.instr_back(1) & units(cellNum).pro.neural.instrDir.ts_pst < prs.instr_back(2);
    instrDir_win = units(cellNum).pro.neural.instrDir.ts_pst > prs.instr_dir(1) & units(cellNum).pro.neural.instrDir.ts_pst <= prs.instr_dir(2); 
    goCue_win = units(cellNum).pro.neural.instrDir.ts_pst > prs.go_cue(1) & units(cellNum).pro.neural.instrDir.ts_pst <= prs.go_cue(2); 
    base_win = units(cellNum).pro.neural.sacc.ts_pst > prs.baseline_win(1) & units(cellNum).pro.neural.sacc.ts_pst <= prs.baseline_win(2);
    base_instr = units(cellNum).pro.neural.instr.ts_pst > prs.baseline_instr(1) & units(cellNum).pro.neural.instr.ts_pst <= prs.baseline_instr(2);
    t_instr = units(cellNum).pro.neural.instr.ts_pst(instr_win); %time
    t_instr_back = units(cellNum).pro.neural.instrDir.ts_pst(instr_back_win); 
    t_sacc = units(cellNum).pro.neural.sacc.ts_pst(sacc_win);
    t_instrDir = units(cellNum).pro.neural.instrDir.ts_pst(instrDir_win); 
    t_goCue = units(cellNum).pro.neural.instrDir.ts_pst(goCue_win);
    sacc_on = units(cellNum).pro.neural.sacc.ts_pst > 0 & units(cellNum).pro.neural.sacc.ts_pst <= prs.saccade_win(2);
    % sacc_pre = units(cellNum).pro.neural.sacc.ts_pst > prs.saccade_win(1) & units(cellNum).pro.neural.sacc.ts_pst <= -0.050;
    sacc_pre = units(cellNum).pro.neural.sacc.ts_pst >= -0.150 & units(cellNum).pro.neural.sacc.ts_pst <= 0;
    
    %get spks pro
    instr_spks_pro = units(cellNum).pro.neural.instr.rate_pst(instr_win);
    instr_back_spks_pro = units(cellNum).pro.neural.instrDir.rate_pst(instr_back_win);
    sacc_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(sacc_win);
    instrDir_spks_pro = units(cellNum).pro.neural.instrDir.rate_pst(instrDir_win);
    goCue_spks_pro = units(cellNum).pro.neural.instrDir.rate_pst(goCue_win); 
    sacc_on_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(sacc_on);
    sacc_pre_spks_pro = units(cellNum).pro.neural.instr.rate_pst(sacc_pre);
    base_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(base_win);
    base_spks_instr = units(cellNum).pro.neural.instr.rate_pst(base_instr);
    
    %pst, mean and sem
    %instr
    units(cellNum).pro.neural.instr.ts_pst_win = units(cellNum).pro.neural.instr.ts_pst(instr_win);
    units(cellNum).pro.neural.instr.rate_pst_win = units(cellNum).pro.neural.instr.rate_pst(instr_win);
    units(cellNum).pro.neural.instr.rate_mu = mean(instr_spks_pro);
    units(cellNum).pro.neural.instr.rate_sig = std(instr_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.instr.rate_std = std(instr_spks_pro);
    
    %instr back
    units(cellNum).pro.neural.instr_back.ts_pst_win = units(cellNum).pro.neural.instrDir.ts_pst(instr_back_win);
    units(cellNum).pro.neural.instr_back.rate_pst_win = units(cellNum).pro.neural.instrDir.rate_pst(instr_back_win);
    units(cellNum).pro.neural.instr_back.rate_mu = mean(instr_back_spks_pro);
    units(cellNum).pro.neural.instr_back.rate_sig = std(instr_back_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.instr_back.rate_std = std(instr_back_spks_pro);
    
    %sacc
    units(cellNum).pro.neural.sacc.ts_pst_win = units(cellNum).pro.neural.sacc.ts_pst(sacc_win);
    units(cellNum).pro.neural.sacc.rate_pst_win = units(cellNum).pro.neural.sacc.rate_pst(sacc_win);
    units(cellNum).pro.neural.sacc.rate_mu = mean(sacc_spks_pro);
    units(cellNum).pro.neural.sacc.rate_sig = std(sacc_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.sacc.rate_std = std(sacc_spks_pro);
    
    % instr Dir
    units(cellNum).pro.neural.instrDir.ts_pst_win = units(cellNum).pro.neural.instrDir.ts_pst(instrDir_win);
    units(cellNum).pro.neural.instrDir.rate_pst_win = units(cellNum).pro.neural.instrDir.rate_pst(instrDir_win);
    units(cellNum).pro.neural.instrDir.rate_mu = mean(instrDir_spks_pro);
    units(cellNum).pro.neural.instrDir.rate_sig = std(instrDir_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.instrDir.rate_std = std(instrDir_spks_pro);
    
    % go Cue
    units(cellNum).pro.neural.goCue.ts_pst_win = units(cellNum).pro.neural.instrDir.ts_pst(goCue_win); 
    units(cellNum).pro.neural.goCue.rate_pst_win = units(cellNum).pro.neural.instrDir.rate_pst(goCue_win);
    units(cellNum).pro.neural.goCue.rate_mu = mean(goCue_spks_pro);
    units(cellNum).pro.neural.goCue.rate_sig = std(goCue_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.goCue.rate_std = std(goCue_spks_pro);
    
    %base
    units(cellNum).pro.neural.base.ts_pst_win = units(cellNum).pro.neural.sacc.ts_pst(base_win);
    units(cellNum).pro.neural.base.rate_pst_win = units(cellNum).pro.neural.sacc.rate_pst(base_win);
    units(cellNum).pro.neural.base.rate_mu = mean(base_spks_pro);
    units(cellNum).pro.neural.base.rate_sig = std(base_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.base.rate_std = std(base_spks_pro);
    
    % base instr pro
    units(cellNum).pro.neural.base.ts_pst_win_instr = units(cellNum).pro.neural.instr.ts_pst(base_instr);
    units(cellNum).pro.neural.base.rate_instr_pst_win = units(cellNum).pro.neural.instr.rate_pst(base_instr);
    units(cellNum).pro.neural.base.rate_instr_mu = mean(base_spks_instr);
    units(cellNum).pro.neural.base.rate_instr_sig = std(base_spks_instr)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.base.rate_instr_std = std(base_spks_instr);
    
    % sacc after 0 > or < than baseline
    if mean(units(cellNum).pro.neural.sacc.nspk) > mean(units(cellNum).pro.neural.presacc.nspk)
        units(cellNum).pro.neural.exc = 1; units(cellNum).pro.neural.sup = 0;
    else
        units(cellNum).pro.neural.exc = 0; units(cellNum).pro.neural.sup = 1;
    end
    
     % instr > or < than baseline
     if mean(units(cellNum).pro.neural.instr.nspk) > mean(units(cellNum).pro.neural.base.nspk)
         units(cellNum).pro.neural.exc_instr = 1; units(cellNum).pro.neural.sup_instr = 0;
     else
         units(cellNum).pro.neural.exc_instr = 0; units(cellNum).pro.neural.sup_instr = 1;
     end
    
     % instr > or < than baseline
     if mean(units(cellNum).pro.neural.instr_back.nspk) > mean(units(cellNum).pro.neural.base.nspk)
         units(cellNum).pro.neural.exc_instr_back = 1; units(cellNum).pro.neural.sup_instr_back = 0;
     else
         units(cellNum).pro.neural.exc_instr_back = 0; units(cellNum).pro.neural.sup_instr_back = 1;
     end
     
     %% Detect latencies pro
     
%      % instr pro
%      if units(cellNum).pro.neural.exc_instr
%          [units(cellNum).pro.neural.instr.peak_resp, indx_max] = max(units(cellNum).pro.neural.instr.rate_pst_win);
%          units(cellNum).pro.neural.instr.peak_resp_time = t_instr(indx_max);
%          
%          [units(cellNum).pro.neural.instr.min_resp, indx_min] = min(units(cellNum).pro.neural.instr.rate_pst_win);
%          units(cellNum).pro.neural.instr.min_resp_time = t_instr(indx_min);
%      else
%          [units(cellNum).pro.neural.instr.peak_resp, indx_max] = max(flip(units(cellNum).pro.neural.instr.rate_pst_win));
%          units(cellNum).pro.neural.instr.peak_resp_time = t_instr(indx_max);
%          
%          [units(cellNum).pro.neural.instr.min_resp, indx_min] = min(units(cellNum).pro.neural.instr.rate_pst_win);
%          units(cellNum).pro.neural.instr.min_resp_time = t_instr(indx_min);
%      end
     
     % instr back pro
     if units(cellNum).pro.neural.exc_instr_back
         [units(cellNum).pro.neural.instr_back.peak_resp, indx_max] = max(units(cellNum).pro.neural.instr_back.rate_pst_win);
         units(cellNum).pro.neural.instr_back.peak_resp_time = t_instr_back(indx_max);
         
         [units(cellNum).pro.neural.instr_back.min_resp, indx_min] = min(units(cellNum).pro.neural.instr_back.rate_pst_win);
         units(cellNum).pro.neural.instr_back.min_resp_time = t_instr_back(indx_min);
     else
         [units(cellNum).pro.neural.instr_back.peak_resp, indx_max] = max(flip(units(cellNum).pro.neural.instr_back.rate_pst_win));
         units(cellNum).pro.neural.instr_back.peak_resp_time = t_instr_back(indx_max);
         
         [units(cellNum).pro.neural.instr_back.min_resp, indx_min] = min(flip(units(cellNum).pro.neural.instr_back.rate_pst_win));
         units(cellNum).pro.neural.instr_back.min_resp_time = t_instr_back(indx_min);
     end
     
     % sacc pro
     if units(cellNum).pro.neural.exc
         [units(cellNum).pro.neural.sacc.peak_resp, indx_max] = max(units(cellNum).pro.neural.sacc.rate_pst_win);
         units(cellNum).pro.neural.sacc.peak_resp_time = t_sacc(indx_max);
         
         [units(cellNum).pro.neural.sacc.min_resp, indx_min] = min(units(cellNum).pro.neural.sacc.rate_pst_win);
         units(cellNum).pro.neural.sacc.min_resp_time = t_sacc(indx_min);
         
     else
         [units(cellNum).pro.neural.sacc.peak_resp, indx_max] = max(flip(units(cellNum).pro.neural.sacc.rate_pst_win));
         units(cellNum).pro.neural.sacc.peak_resp_time = t_sacc(indx_max);
         
         [units(cellNum).pro.neural.sacc.min_resp, indx_min] = min(flip(units(cellNum).pro.neural.sacc.rate_pst_win));
         units(cellNum).pro.neural.sacc.min_resp_time = t_sacc(indx_min);
     end
     
     %%
     %get spks anti
    instr_spks_anti = units(cellNum).anti.neural.instr.rate_pst(instr_win);
    sacc_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(sacc_win);
    instr_back_spks_anti = units(cellNum).anti.neural.instrDir.rate_pst(instr_back_win);
    instrDir_spks_anti = units(cellNum).anti.neural.instrDir.rate_pst(instrDir_win);
    goCue_spks_anti = units(cellNum).anti.neural.instrDir.rate_pst(goCue_win); 
    base_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(base_win);
    sacc_on_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(sacc_on);
    sacc_pre_spks_anti = units(cellNum).anti.neural.instr.rate_pst(sacc_pre);
    base_spks_instr = units(cellNum).anti.neural.instr.rate_pst(base_instr);
    
    %pst, mean and sem
    %instr
    units(cellNum).anti.neural.instr.ts_pst_win = units(cellNum).anti.neural.instr.ts_pst(instr_win);
    units(cellNum).anti.neural.instr.rate_pst_win = units(cellNum).anti.neural.instr.rate_pst(instr_win);
    units(cellNum).anti.neural.instr.rate_mu = mean(instr_spks_anti);
    units(cellNum).anti.neural.instr.rate_sig = std(instr_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.instr.rate_std = std(instr_spks_anti);
    
    % sacc
    units(cellNum).anti.neural.sacc.ts_pst_win = units(cellNum).anti.neural.sacc.ts_pst(sacc_win);
    units(cellNum).anti.neural.sacc.rate_pst_win = units(cellNum).anti.neural.sacc.rate_pst(sacc_win);
    units(cellNum).anti.neural.sacc.rate_mu = mean(sacc_spks_anti);
    units(cellNum).anti.neural.sacc.rate_sig = std(sacc_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.sacc.rate_std = std(sacc_spks_anti);
    
        %instr back
    units(cellNum).anti.neural.instr_back.ts_pst_win = units(cellNum).anti.neural.instrDir.ts_pst(instr_back_win);
    units(cellNum).anti.neural.instr_back.rate_pst_win = units(cellNum).anti.neural.instrDir.rate_pst(instr_back_win);
    units(cellNum).anti.neural.instr_back.rate_mu = mean(instr_back_spks_anti);
    units(cellNum).anti.neural.instr_back.rate_sig = std(instr_back_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.instr_back.rate_std = std(instr_back_spks_anti);
    
     % instr Dir
    units(cellNum).anti.neural.instrDir.ts_pst_win = units(cellNum).anti.neural.instrDir.ts_pst(instrDir_win);
    units(cellNum).anti.neural.instrDir.rate_pst_win = units(cellNum).anti.neural.instrDir.rate_pst(instrDir_win);
    units(cellNum).anti.neural.instrDir.rate_mu = mean(instrDir_spks_anti);
    units(cellNum).anti.neural.instrDir.rate_sig = std(instrDir_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.instrDir.rate_std = std(instrDir_spks_anti);
    
    % go Cue
    units(cellNum).anti.neural.goCue.ts_pst_win = units(cellNum).anti.neural.instrDir.ts_pst(goCue_win); 
    units(cellNum).anti.neural.goCue.rate_pst_win = units(cellNum).anti.neural.instrDir.rate_pst(goCue_win);
    units(cellNum).anti.neural.goCue.rate_mu = mean(goCue_spks_anti);
    units(cellNum).anti.neural.goCue.rate_sig = std(goCue_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.goCue.rate_std = std(goCue_spks_anti);
    
    % base
    units(cellNum).anti.neural.base.ts_pst_win = units(cellNum).anti.neural.sacc.ts_pst(base_win);
    units(cellNum).anti.neural.base.rate_pst_win = units(cellNum).anti.neural.sacc.rate_pst(base_win);
    units(cellNum).anti.neural.base.rate_mu = mean(base_spks_anti);
    units(cellNum).anti.neural.base.rate_sig = std(base_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.base.rate_std = std(base_spks_anti);
    
    % base instr
    units(cellNum).anti.neural.base.ts_pst_win_instr = units(cellNum).anti.neural.instr.ts_pst(base_instr);
    units(cellNum).anti.neural.base.rate_instr_pst_win = units(cellNum).anti.neural.instr.rate_pst(base_instr);
    units(cellNum).anti.neural.base.rate_instr_mu = mean(base_spks_instr);
    units(cellNum).anti.neural.base.rate_instr_sig = std(base_spks_instr)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.base.rate_instr_std = std(base_spks_instr);
   
    % sacc after 0 > or < than baseline
    if  mean(units(cellNum).anti.neural.sacc.nspk) > mean(units(cellNum).anti.neural.presacc.nspk)
        units(cellNum).anti.neural.exc = 1; units(cellNum).anti.neural.sup = 0;
    else
        units(cellNum).anti.neural.exc = 0; units(cellNum).anti.neural.sup = 1;
    end
   
    % instr > or < than baseline
     if mean(units(cellNum).anti.neural.instr.nspk) > mean(units(cellNum).anti.neural.base.nspk)
         units(cellNum).anti.neural.exc_instr = 1; units(cellNum).anti.neural.sup_instr = 0;
     else
         units(cellNum).anti.neural.exc_instr = 0; units(cellNum).anti.neural.sup_instr = 1;
     end
    
   % instr > or < than baseline
     if mean(units(cellNum).anti.neural.instr_back.nspk) > mean(units(cellNum).anti.neural.base.nspk)
         units(cellNum).anti.neural.exc_instr_back = 1; units(cellNum).anti.neural.sup_instr_back = 0;
     else
         units(cellNum).anti.neural.exc_instr_back = 0; units(cellNum).anti.neural.sup_instr_back = 1;
     end
    
    %% Detect latencies
    
    % instr anti
%     if units(cellNum).anti.neural.exc_instr
%         [units(cellNum).anti.neural.instr.peak_resp, indx_max] = max(units(cellNum).anti.neural.instr.rate_pst_win);
%         units(cellNum).anti.neural.instr.peak_resp_time = t_instr(indx_max);
%         
%         [units(cellNum).anti.neural.instr.min_resp, indx_min] = min(units(cellNum).anti.neural.instr.rate_pst_win);
%         units(cellNum).anti.neural.instr.min_resp_time = t_instr(indx_min);
%     else
%         [units(cellNum).anti.neural.instr.peak_resp, indx_max] = max(flip(units(cellNum).anti.neural.instr.rate_pst_win));
%         units(cellNum).anti.neural.instr.peak_resp_time = t_instr(indx_max);
%         
%         [units(cellNum).anti.neural.instr.min_resp, indx_min] = min(units(cellNum).anti.neural.instr.rate_pst_win);
%         units(cellNum).anti.neural.instr.min_resp_time = t_instr(indx_min);
%     end

% instr back pro
     if units(cellNum).anti.neural.exc_instr_back
         [units(cellNum).anti.neural.instr_back.peak_resp, indx_max] = max(units(cellNum).anti.neural.instr_back.rate_pst_win);
         units(cellNum).anti.neural.instr_back.peak_resp_time = t_instr_back(indx_max);
         
         [units(cellNum).anti.neural.instr_back.min_resp, indx_min] = min(units(cellNum).anti.neural.instr_back.rate_pst_win);
         units(cellNum).anti.neural.instr_back.min_resp_time = t_instr_back(indx_min);
     else
         [units(cellNum).anti.neural.instr_back.peak_resp, indx_max] = max(flip(units(cellNum).anti.neural.instr_back.rate_pst_win));
         units(cellNum).anti.neural.instr_back.peak_resp_time = t_instr_back(indx_max);
         
         [units(cellNum).anti.neural.instr_back.min_resp, indx_min] = min(flip(units(cellNum).anti.neural.instr_back.rate_pst_win));
         units(cellNum).anti.neural.instr_back.min_resp_time = t_instr_back(indx_min);
     end
    
    
    % sacc anti
    if units(cellNum).anti.neural.exc
        [units(cellNum).anti.neural.sacc.peak_resp, indx_max] = max(units(cellNum).anti.neural.sacc.rate_pst_win);
        units(cellNum).anti.neural.sacc.peak_resp_time = t_sacc(indx_max);
        
        [units(cellNum).anti.neural.sacc.min_resp, indx_min] = min(units(cellNum).anti.neural.sacc.rate_pst_win);
        units(cellNum).anti.neural.sacc.min_resp_time = t_sacc(indx_min);
    else
        [units(cellNum).anti.neural.sacc.peak_resp, indx_max] = max(flip(units(cellNum).anti.neural.sacc.rate_pst_win));
        units(cellNum).anti.neural.sacc.peak_resp_time = t_sacc(indx_max);
        
        [units(cellNum).anti.neural.sacc.min_resp, indx_min] = min(flip(units(cellNum).anti.neural.sacc.rate_pst_win));
        units(cellNum).anti.neural.sacc.min_resp_time = t_sacc(indx_min);
    end
    
    %% compare windows against baseline activity for pro and anti - nspk
    % pro
    [units(cellNum).stats.pro.pval.instrVSbase_nspk, units(cellNum).stats.pro.flags.instrVSbase_nspk] = signrank(units(cellNum).pro.neural.instr.nspk,units(cellNum).pro.neural.base.nspk_instr);
    [units(cellNum).stats.pro.pval.instr_backVSbase_nspk, units(cellNum).stats.pro.flags.instr_backVSbase_nspk] = signrank(units(cellNum).pro.neural.instr_back.nspk,units(cellNum).pro.neural.base.nspk_instr);
    %[units(cellNum).stats.pro.pval.saccVSbase_nspk, units(cellNum).stats.pro.flags.saccVSbase_nspk] = signrank(units(cellNum).pro.neural.sacc.nspk,units(cellNum).pro.neural.base.nspk);
    [units(cellNum).stats.pro.pval.saccVSbase_nspk, units(cellNum).stats.pro.flags.saccVSbase_nspk] = signrank(units(cellNum).pro.neural.sacc.nspk,units(cellNum).pro.neural.base.nspk_instr);
    [units(cellNum).stats.pro.pval.saccVSinstrDir_nspk, units(cellNum).stats.pro.flags.saccVSinstrDir_nspk] = signrank(units(cellNum).pro.neural.sacc.nspk,units(cellNum).pro.neural.instrDir.nspk);
    [units(cellNum).stats.pro.pval.instrDirVSbase_nspk, units(cellNum).stats.pro.flags.instrDirVSbase_nspk] = signrank(units(cellNum).pro.neural.instrDir.nspk,units(cellNum).pro.neural.base.nspk_instr);
    [units(cellNum).stats.pro.pval.goCueVSinstrDir_nspk, units(cellNum).stats.pro.flags.goCueVSinstrDir_nspk] = signrank(units(cellNum).pro.neural.goCue.nspk,units(cellNum).pro.neural.instrDir.nspk);
    
    % stat for comparing exc and sup
    [units(cellNum).pro.neural.categ.flag, units(cellNum).pro.neural.categ.pval] = ttest2(units(cellNum).pro.neural.presacc.nspk, units(cellNum).pro.neural.sacc.nspk); 
    
    % anti
    [units(cellNum).stats.anti.pval.instrVSbase_nspk, units(cellNum).stats.anti.flags.instrVSbase_nspk] = signrank(units(cellNum).anti.neural.instr.nspk,units(cellNum).anti.neural.base.nspk_instr);
    [units(cellNum).stats.anti.pval.instr_backVSbase_nspk, units(cellNum).stats.anti.flags.instr_backVSbase_nspk] = signrank(units(cellNum).anti.neural.instr_back.nspk,units(cellNum).anti.neural.base.nspk_instr);
    % [units(cellNum).stats.anti.pval.saccVSbase_nspk, units(cellNum).stats.anti.flags.saccVSbase_nspk] = signrank(units(cellNum).anti.neural.sacc.nspk,units(cellNum).anti.neural.base.nspk);
    [units(cellNum).stats.anti.pval.saccVSbase_nspk, units(cellNum).stats.anti.flags.saccVSbase_nspk] = signrank(units(cellNum).anti.neural.sacc.nspk,units(cellNum).anti.neural.base.nspk_instr);
    [units(cellNum).stats.anti.pval.saccVSinstrDir_nspk, units(cellNum).stats.anti.flags.saccVSinstrDir_nspk] = signrank(units(cellNum).anti.neural.sacc.nspk,units(cellNum).anti.neural.instrDir.nspk);
    [units(cellNum).stats.anti.pval.instrDirVSinstr_nspk, units(cellNum).stats.anti.flags.instrDirVSinstr_nspk] = signrank(units(cellNum).anti.neural.instrDir.nspk,units(cellNum).anti.neural.instr.nspk);
    [units(cellNum).stats.anti.pval.goCueVSinstrDir_nspk, units(cellNum).stats.anti.flags.goCueVSinstrDir_nspk] = signrank(units(cellNum).anti.neural.goCue.nspk,units(cellNum).anti.neural.instrDir.nspk);
    % stat for comparing exc and supp
    [units(cellNum).anti.neural.categ.flag, units(cellNum).anti.neural.categ.pval] = ttest2(units(cellNum).anti.neural.presacc.nspk, units(cellNum).anti.neural.sacc.nspk); 
    
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
    [units(cellNum).stats.instr_back.flags.proVsAnti_instr, units(cellNum).stats.instr_back.pval.proVsAnti_instr] = ttest2(units(cellNum).pro.neural.instr_back.nspk,units(cellNum).anti.neural.instr_back.nspk);
    [units(cellNum).stats.sacc.flags.proVsAnti_sacc, units(cellNum).stats.sacc.pval.proVsAnti_sacc] = ttest2(units(cellNum).pro.neural.sacc.nspk,units(cellNum).anti.neural.sacc.nspk);
    
    %kstest2
    [units(cellNum).stats.instr.flags.proVsAnti_instr_ks_nspk, units(cellNum).stats.instr.pval.proVsAnti_instr_ks_nspk] = kstest2(units(cellNum).pro.neural.instr.nspk , units(cellNum).anti.neural.instr.nspk);
    [units(cellNum).stats.instr_back.flags.proVsAnti_instr_ks_nspk, units(cellNum).stats.instr_back.pval.proVsAnti_instr_ks_nspk] = kstest2(units(cellNum).pro.neural.instr_back.nspk , units(cellNum).anti.neural.instr_back.nspk);
    [units(cellNum).stats.sacc.flags.proVsAnti_sacc_ks_nspk, units(cellNum).stats.sacc.pval.proVsAnti_sacc_ks_nspk] = kstest2(units(cellNum).pro.neural.sacc.nspk , units(cellNum).anti.neural.sacc.nspk);
    
    [units(cellNum).stats.instr.flags.proVsAnti_instr_ks_t_spk, units(cellNum).stats.instr.pval.proVsAnti_instr_ks_t_spk] = kstest2(units(cellNum).pro.neural.instr.tspk_all , units(cellNum).anti.neural.instr.tspk_all);
    [units(cellNum).stats.sacc.flags.proVsAnti_sacc_ks_t_spk, units(cellNum).stats.sacc.pval.proVsAnti_sacc_ks_t_spk] = kstest2(units(cellNum).pro.neural.sacc.tspk_all , units(cellNum).anti.neural.sacc.tspk_all);
    
    
    %% Compute change in FR from mean and baseline
    % pro
    units(cellNum).pro.neural.instr.delta_rate = units(cellNum).pro.neural.instr.rate_pst_win - units(cellNum).pro.neural.base.rate_instr_mu; % mean
    units(cellNum).pro.neural.sacc.delta_rate = units(cellNum).pro.neural.sacc.rate_pst_win - units(cellNum).pro.neural.sacc.rate_mu; % mean
    
    units(cellNum).pro.neural.sacc.delta_rate_base = units(cellNum).pro.neural.sacc.rate_pst_win - units(cellNum).pro.neural.base.rate_mu;  %baseline
    units(cellNum).pro.neural.instr.delta_rate_base = units(cellNum).pro.neural.instr.rate_pst_win - units(cellNum).pro.neural.base.rate_instr_mu; % baseline from ITI
    
    % anti
    units(cellNum).anti.neural.instr.delta_rate = units(cellNum).anti.neural.instr.rate_pst_win - units(cellNum).anti.neural.base.rate_instr_mu; % mean
    units(cellNum).anti.neural.sacc.delta_rate = units(cellNum).anti.neural.sacc.rate_pst_win - units(cellNum).anti.neural.sacc.rate_mu;
    
    units(cellNum).anti.neural.sacc.delta_rate_base = units(cellNum).anti.neural.sacc.rate_pst_win - units(cellNum).anti.neural.base.rate_mu; % baseline
    
    units(cellNum).anti.neural.instr.delta_rate_base = units(cellNum).anti.neural.instr.rate_pst_win - units(cellNum).anti.neural.base.rate_instr_mu; % baseline from ITI
    
    
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
    
    
    %% normalized FR - z-scored aligned to sacc     NEEDS TO CHANGE - CALCULATE ONLY IN TIME WINDOW
    % pro --  aligned to saccade
    units(cellNum).pro.neural.sacc.norm.rate_pst = (units(cellNum).pro.neural.sacc.rate_pst - mean(units(cellNum).pro.neural.sacc.rate_pst))/...
        std(units(cellNum).pro.neural.sacc.rate_pst);
    % anti --  aligned to saccade
    units(cellNum).anti.neural.sacc.norm.rate_pst = (units(cellNum).anti.neural.sacc.rate_pst - mean(units(cellNum).anti.neural.sacc.rate_pst))/...
        std(units(cellNum).anti.neural.sacc.rate_pst);
    % normalize anti based on prosaccades (as described by Maarten, March 2019)
    units(cellNum).anti.neural.sacc.norm.rate_pst_pro = (units(cellNum).anti.neural.sacc.rate_pst - mean(units(cellNum).pro.neural.sacc.rate_pst))/...
        std(units(cellNum).pro.neural.sacc.rate_pst);
    % grab normalized baseline vals
    norm_pro_base = mean(units(cellNum).pro.neural.sacc.norm.rate_pst(base_win));
    % compute change in firing rate for pro and anti with normalized 
    units(cellNum).pro.neural.sacc.norm.delta_rate = units(cellNum).pro.neural.sacc.norm.rate_pst(sacc_win) - norm_pro_base;
    units(cellNum).anti.neural.sacc.norm.delta_rate = units(cellNum).anti.neural.sacc.norm.rate_pst_pro(sacc_win) - norm_pro_base;
    
    %% Normalized FR only in time window using pro as norm base (will make a difference from above as the (baseline or reference) mean is calculated only in the time window
     % pro --  aligned to saccade
    units(cellNum).pro.neural.sacc.norm.rate_pst_win = (units(cellNum).pro.neural.sacc.rate_pst_win - (mean(units(cellNum).pro.neural.sacc.rate_pst_win)))/...
        (std(units(cellNum).pro.neural.sacc.rate_pst_win));
    units(cellNum).anti.neural.sacc.norm.rate_pst_pro_win = (units(cellNum).anti.neural.sacc.rate_pst_win - (mean(units(cellNum).pro.neural.sacc.rate_pst_win)))/...
        (std(units(cellNum).pro.neural.sacc.rate_pst_win));
    
    %% normalized FR - z-scored aligned to instr onset NEEDS TO CHANGE - CALCULATE ONLY IN TIME WINDOW (?)
    % pro --  aligned to instr
    units(cellNum).pro.neural.instr.norm.rate_pst = (units(cellNum).pro.neural.instr.rate_pst - mean(units(cellNum).pro.neural.instr.rate_pst))/...
        std(units(cellNum).pro.neural.instr.rate_pst);
    % anti --  aligned to instr
    units(cellNum).anti.neural.instr.norm.rate_pst = (units(cellNum).anti.neural.instr.rate_pst - mean(units(cellNum).anti.neural.instr.rate_pst))/...
        std(units(cellNum).anti.neural.instr.rate_pst);
    % normalize anti based on prosaccades (as described by Maarten, March 2019)
    units(cellNum).anti.neural.instr.norm.rate_pst_pro = (units(cellNum).anti.neural.instr.rate_pst - mean(units(cellNum).pro.neural.instr.rate_pst))/...
        std(units(cellNum).pro.neural.instr.rate_pst);
    % grab normalized baseline vals
    norm_pro_base = mean(units(cellNum).pro.neural.instr.norm.rate_pst(base_instr));
    % compute change in firing rate for pro and anti with normalized 
    units(cellNum).pro.neural.instr.norm.delta_rate = units(cellNum).pro.neural.instr.norm.rate_pst(instr_win) - norm_pro_base;
    units(cellNum).anti.neural.instr.norm.delta_rate = units(cellNum).anti.neural.instr.norm.rate_pst_pro(instr_win) - norm_pro_base; 
    
    %% compute change index in norm -- CHECK! 
%     units(cellNum).stats.instr_back.change_indx = (units(cellNum).pro.neural.instr_back.rate_mu  - units(cellNum).anti.neural.instr_back.rate_mu)/(units(cellNum).pro.neural.instr_back.rate_mu  + units(cellNum).anti.neural.instr_back.rate_mu);
%     units(cellNum).stats.sacc.change_indx = (units(cellNum).pro.neural.sacc.rate_mu  - units(cellNum).anti.neural.sacc.rate_mu)/(units(cellNum).pro.neural.sacc.rate_mu  + units(cellNum).anti.neural.sacc.rate_mu);
    
    %% compute modulation depth (or index) for each cell for pro and anti (max-min)/(max+min) sacc-base/sacc+base
     % instr
     units(cellNum).stats.pro.instr_back.mod_depth = (units(cellNum).pro.neural.instr_back.peak_resp - units(cellNum).pro.neural.base.rate_instr_mu)/ (units(cellNum).pro.neural.instr_back.peak_resp + units(cellNum).pro.neural.base.rate_instr_mu);
     units(cellNum).stats.anti.instr_back.mod_depth = (units(cellNum).anti.neural.instr_back.peak_resp - units(cellNum).anti.neural.base.rate_instr_mu)/ (units(cellNum).anti.neural.instr_back.peak_resp + units(cellNum).anti.neural.base.rate_instr_mu);
     
     % sacc
     units(cellNum).stats.pro.sacc.mod_depth = (units(cellNum).pro.neural.sacc.peak_resp - mean(units(cellNum).pro.neural.presacc.nspk))/(units(cellNum).pro.neural.sacc.peak_resp + mean(units(cellNum).pro.neural.presacc.nspk));
     units(cellNum).stats.anti.sacc.mod_depth = (units(cellNum).anti.neural.sacc.peak_resp - mean(units(cellNum).anti.neural.presacc.nspk))/(units(cellNum).anti.neural.sacc.peak_resp + mean(units(cellNum).anti.neural.presacc.nspk));
     
     %% Just another mod index for pro and anti separately (R_max - R_min)/R_max
     % instr
     units(cellNum).stats.pro.instr.selectivity_indx = (units(cellNum).pro.neural.instr_back.peak_resp - units(cellNum).pro.neural.base.rate_instr_mu)/units(cellNum).pro.neural.instr_back.peak_resp;
     units(cellNum).stats.anti.instr.selectivity_indx = (units(cellNum).anti.neural.instr_back.peak_resp - units(cellNum).anti.neural.base.rate_instr_mu)/units(cellNum).anti.neural.instr_back.peak_resp;
     
     % sacc
     units(cellNum).stats.pro.sacc.selectivity_indx = (units(cellNum).pro.neural.sacc.peak_resp - mean(units(cellNum).pro.neural.presacc.nspk))/units(cellNum).pro.neural.sacc.peak_resp;
     units(cellNum).stats.anti.sacc.selectivity_indx = (units(cellNum).anti.neural.sacc.peak_resp - mean(units(cellNum).pro.neural.presacc.nspk))/units(cellNum).anti.neural.sacc.peak_resp;

        %% Compute modulation (another) mod index as describred in Wypych et al 2012. 
%     % instr
%     units(cellNum).stats.pro.instr_back.mod_indx_mu = (units(cellNum).pro.neural.instr_back.nspk)/(units(cellNum).pro.neural.instr_back.rate_mu - units(cellNum).pro.neural.base.rate_instr_mu); 
%     units(cellNum).stats.anti.instr_back.mod_indx_mu = (units(cellNum).anti.neural.instr_back.nspk)/(units(cellNum).anti.neural.instr_back.rate_mu - units(cellNum).anti.neural.base.rate_instr_mu); 
%     
%     % sacc
%     units(cellNum).stats.pro.sacc.mod_indx_mu = (units(cellNum).pro.neural.sacc.nspk)/(units(cellNum).pro.neural.sacc.rate_mu - units(cellNum).pro.neural.base.rate_mu);
%     units(cellNum).stats.anti.sacc.mod_indx_mu = (units(cellNum).anti.neural.sacc.nspk)/(units(cellNum).anti.neural.sacc.rate_mu - units(cellNum).anti.neural.base.rate_mu);
    
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
    
     % instr back
    for i=1:length(units(cellNum).pro.neural.instr_back.spkCount(1,:))
        p=(ntrls_pro*units(cellNum).pro.neural.instr_back.pbDist(i)+ntrls_anti*units(cellNum).anti.neural.instr_back.pbDist(i))/(ntrls_pro+ntrls_anti);
        units(cellNum).stats.instr_back.pval.pbDist_testStat(i) = (units(cellNum).pro.neural.instr_back.pbDist(1,i)-units(cellNum).anti.neural.instr_back.pbDist(1,i))/(sqrt(p*(1-p)*(1/ntrls_pro + 1/ntrls_anti)));
        if units(cellNum).stats.instr_back.pval.pbDist_testStat(i) > prs.signif_criteria | units(cellNum).stats.instr_back.pval.pbDist_testStat(i) < -prs.signif_criteria
            units(cellNum).stats.instr_back.flags.pbDist(i) = 1;
        else
            units(cellNum).stats.instr_back.flags.pbDist(i) = 0;
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
    end
    
    for win_num = 1:length(units(cellNum).pro.neural.sacc.spkCount_win(1,:));
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
%         % anti
%         rate_pst = units(cellNum).anti.neural.sacc.rate_pst;
%         exc_anti=0; sup_anti=0;
%     
%         for indx = indx_beg:indx_end
%             r_mu = mean(rate_pst(round(indx-win_size/2):round(indx+win_size/2)));
%             r_thresh1 = mean(rate_pst(indx-win_size_prev:indx))+...
%                 2*std(rate_pst(indx-win_size_prev:indx));
%             r_thresh2 = mean(rate_pst(indx-win_size_prev:indx))-...
%                 2*std(rate_pst(indx-win_size_prev:indx));
%             if (r_thresh1>1 && r_mu>r_thresh1), exc_anti=exc_anti+1;else exc_anti=0;end
%             if (r_thresh2>1 && r_mu<r_thresh2), sup_anti=sup_anti+1;else sup_anti=0;end
%             if exc_anti==3, events.rise.t_on = t(indx); events.rise.type_on = 'exc';
%                 break;
%             elseif sup_anti==3, units(cellNum).pro.events.rise.t_on = t(indx); units(cellNum).pro.events.rise.type_on = 'sup';
%                 break;
%             end
%         end
    
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


%% Compute comparison of FR for all neurons using nspk instr (this doesnt make sense, it needs to be divided by area)
spks_pro = []; spks_anti = []; 
for cellNum = 1:length(units)
    spks_pro = [spks_pro ; units(cellNum).pro.neural.instr.nspk];
    spks_anti = [spks_anti ; units(cellNum).anti.neural.instr.nspk];
end

pop.pro.instr.all_spk = spks_pro; pop.anti.instr.all_spk = spks_anti; 
% mean 
pop.pro.instr.nspk_mu = mean(spks_pro); pop.anti.instr.nspk_mu = mean(spks_anti); % mean
pop.pro.instr.nspk_std = std(spks_pro); pop.anti.instr.nspk_std = std(spks_anti); % std dev


% compute stats using kstest 
[pop.stats.instr.pro_anti_nspk_pVal, ~] = ranksum(spks_pro, spks_anti); 

%% Compute comparison of FR for all neurons using nspk sacc
spks_pro = []; spks_anti = []; 
for cellNum = 1:length(units)
    spks_pro = [spks_pro ; units(cellNum).pro.neural.sacc.nspk];
    spks_anti = [spks_anti ; units(cellNum).anti.neural.sacc.nspk];
end

pop.pro.sacc.all_spk = spks_pro; pop.anti.sacc.all_spk = spks_anti; 

% mean 
pop.pro.sacc.nspk_mu = mean(spks_pro); pop.anti.sacc.nspk_mu = mean(spks_anti); % mean
pop.pro.sacc.nspk_std = std(spks_pro); pop.anti.sacc.nspk_std = std(spks_anti); % std dev

% compute stats using kstest 
[pop.stats.sacc.pro_anti_nspk_pVal, ~] = ranksum(spks_pro, spks_anti);

%% regress activity to saccade kinematics (per area)

% % gather data vermis
% recArea = 'vermis';
% for cellNum = 1:length(units)
%     indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
% end
% indx_area = find(indx_area);
% 
% r_all_pro = []; kin_all_pro = []; r_all_anti = []; kin_all_anti = [];      
% 
% for i = 1:length(indx_area)
%     r_pro = []; amp_pro = []; dur_pro = []; pv_pro = []; rt_pro = [];  r_anti = []; amp_anti = []; dur_anti = []; pv_anti = []; rt_anti = [];
%     
%     
%     % run for pro and anti separately (unequal nr of trials)
%     for j = 1:length(units(indx_area(i)).pro.behav.trial)
%         r_pro(j) = units(indx_area(i)).pro.neural.sacc.nspk(j);
%         amp_pro(j) = units(indx_area(i)).pro.behav.trial(j).saccAmplitude;
%         dur_pro(j) = units(indx_area(i)).pro.behav.trial(j).saccDuration;
%         pv_pro(j) = units(indx_area(i)).pro.behav.trial(j).saccPeakVel;
%         rt_pro(j) = units(indx_area(i)).pro.behav.trial(j).reactionTime;
%     end
%     for j = 1:length(units(indx_area(i)).anti.behav.trial)
%         r_anti(j) = units(indx_area(i)).anti.neural.sacc.nspk(j);
%         amp_anti(j) = units(indx_area(i)).anti.behav.trial(j).saccAmplitude;
%         dur_anti(j) = units(indx_area(i)).anti.behav.trial(j).saccDuration;
%         pv_anti(j) = units(indx_area(i)).anti.behav.trial(j).saccPeakVel;
%         rt_anti(j) = units(indx_area(i)).anti.behav.trial(j).reactionTime;
%     end
%     
%     r_all_pro = [ r_all_pro ; r_pro' ];
%     kin_all_pro = [ kin_all_pro ; amp_pro' dur_pro' pv_pro' rt_pro' ones(size(amp_pro,2),1)];
%     
%     r_all_anti = [ r_all_anti ; r_anti' ];
%     kin_all_anti = [ kin_all_anti ; amp_anti' dur_anti' pv_anti' rt_anti' ones(size(amp_anti,2),1)];
%     
% end
% 
%  % save all r and kin vermis
% pop.pro.vermis.r_all = r_all_pro; pop.pro.vermis.kin_all = kin_all_pro; pop.pro.vermis.labels = [{'amp'} , {'dur'}, {'pv'}, {'rt'}]; 
% pop.anti.vermis.r_all = r_all_anti; pop.anti.vermis.kin_all = kin_all_anti; pop.anti.vermis.labels = [{'amp'} , {'dur'}, {'pv'}, {'rt'}];
% 
% [pop.stats.sacc.pro.regress.vermis.coeff_pro, pop.stats.sacc.pro.regress.vermis.CI_pro, pop.stats.sacc.pro.regress.vermis.rsq_pro, pop.stats.sacc.pro.regress.vermis.reg_stats_pro] = regress(r_all_pro,kin_all_pro); 
% [pop.stats.sacc.anti.regress.vermis.coeff_anti, pop.stats.sacc.anti.regress.vermis.CI_anti, pop.stats.sacc.anti.regress.vermis.rsq_anti, pop.stats.sacc.anti.regress.vermis.reg_stats_anti] = regress(r_all_anti,kin_all_anti); 

 %% gather data lateral
% recArea = 'lateral';
% for cellNum = 1:length(units)
%     indx_area(cellNum) = strcmp(units(cellNum).area, recArea);
% end
% indx_area = find(indx_area);
% 
% r_all_pro = []; kin_all_pro = []; r_all_anti = []; kin_all_anti = [];      
% 
% for i = 1:length(indx_area)
%     r_pro = []; amp_pro = []; dur_pro = []; pv_pro = []; rt_pro = [];  r_anti = []; amp_anti = []; dur_anti = []; pv_anti = []; rt_anti = [];
%     
%     
%     % run for pro and anti separately (unequal nr of trials)
%     for j = 1:length(units(indx_area(i)).pro.behav.trial)
%         r_pro(j) = units(indx_area(i)).pro.neural.sacc.nspk(j);
%         amp_pro(j) = units(indx_area(i)).pro.behav.trial(j).saccAmplitude;
%         dur_pro(j) = units(indx_area(i)).pro.behav.trial(j).saccDuration;
%         pv_pro(j) = units(indx_area(i)).pro.behav.trial(j).saccPeakVel;
%         rt_pro(j) = units(indx_area(i)).pro.behav.trial(j).reactionTime;
%     end
%     for j = 1:length(units(indx_area(i)).anti.behav.trial)
%         r_anti(j) = units(indx_area(i)).anti.neural.sacc.nspk(j);
%         amp_anti(j) = units(indx_area(i)).anti.behav.trial(j).saccAmplitude;
%         dur_anti(j) = units(indx_area(i)).anti.behav.trial(j).saccDuration;
%         pv_anti(j) = units(indx_area(i)).anti.behav.trial(j).saccPeakVel;
%         rt_anti(j) = units(indx_area(i)).anti.behav.trial(j).reactionTime;
%     end
%     
%     r_all_pro = [ r_all_pro ; r_pro' ];
%     kin_all_pro = [ kin_all_pro ; amp_pro' dur_pro' pv_pro' rt_pro' ones(size(amp_pro,2),1)];
%     
%     r_all_anti = [ r_all_anti ; r_anti' ];
%     kin_all_anti = [ kin_all_anti ; amp_anti' dur_anti' pv_anti' rt_anti' ones(size(amp_anti,2),1)];
%     
% end
% 
%  % save all r and kin vermis
% pop.pro.lateral.r_all = r_all_pro; pop.pro.lateral.kin_all = kin_all_pro; pop.pro.lateral.labels = [{'amp'} , {'dur'}, {'pv'}, {'rt'}]; 
% pop.anti.lateral.r_all = r_all_anti; pop.anti.lateral.kin_all = kin_all_anti;  pop.anti.lateral.labels = [{'amp'} , {'dur'}, {'pv'}, {'rt'}]; 
% 
% [pop.stats.sacc.pro.regress.lateral.coeff_pro, pop.stats.sacc.pro.regress.lateral.CI_pro, pop.stats.sacc.pro.regress.lateral.rsq_pro, pop.stats.sacc.pro.regress.lateral.reg_stats_pro] = regress(r_all_pro,kin_all_pro); 
% [pop.stats.sacc.anti.regress.lateral.coeff_anti, pop.stats.sacc.anti.regress.lateral.CI_anti, pop.stats.sacc.anti.regress.lateral.rsq_anti, pop.stats.sacc.anti.regress.lateral.reg_stats_anti] = regress(r_all_anti,kin_all_anti); 


%% regress activity to saccade kinematics per cell
% Reaction time is extra, but not eye kinematic.
pop.stats.sacc.vermis.r_pro_all = []; pop.stats.sacc.vermis.r_anti_all = [] ; pop.stats.sacc.lateral.r_pro_all = []; pop.stats.sacc.lateral.r_anti_all = []; 
for i = 1:length(units)
    r_pro = []; amp_pro = []; dur_pro = []; pv_pro = []; rt_pro = [];  r_anti = []; amp_anti = []; dur_anti = []; pv_anti = []; rt_anti = [];
    r_all_pro = []; kin_all_pro = []; r_all_anti = []; kin_all_anti = [];
    
    % run for pro and anti separately (unequal nr of trials)
    for j = 1:length(units(i).pro.behav.trial)
        r_pro(:,j) = units(i).pro.neural.sacc.nspk(j);
        amp_pro(:,j) = units(i).pro.behav.trial(j).saccAmplitude;
        dur_pro(:,j) = units(i).pro.behav.trial(j).saccDuration;
        pv_pro(:,j) = units(i).pro.behav.trial(j).saccPeakVel;
        % rt_pro(:,j) = units(i).pro.behav.trial(j).reactionTime;
    end
    for j = 1:length(units(i).anti.behav.trial)
        r_anti(:,j) = units(i).anti.neural.sacc.nspk(j); 
        amp_anti(:,j) = units(i).anti.behav.trial(j).saccAmplitude;
        dur_anti(:,j) = units(i).anti.behav.trial(j).saccDuration;
        pv_anti(:,j) = units(i).anti.behav.trial(j).saccPeakVel;
        % rt_anti(:,j) = units(i).anti.behav.trial(j).reactionTime;
    end
    
    r_all_pro = [ r_all_pro ; r_pro' ];
    % kin_all_pro = [ kin_all_pro ; amp_pro' dur_pro' pv_pro' rt_pro' ones(size(amp_pro,2),1)];
    kin_all_pro = [ kin_all_pro ; amp_pro' dur_pro' pv_pro' ones(size(amp_pro,2),1)];
    
    
    r_all_anti = [ r_all_anti ; r_anti' ];
    % kin_all_anti = [ kin_all_anti ; amp_anti' dur_anti' pv_anti' rt_anti' ones(size(amp_anti,2),1)];
    kin_all_anti = [ kin_all_anti ; amp_anti' dur_anti' pv_anti' ones(size(amp_anti,2),1)];
    
    % save all r and eye kinmeatics
    %     pop.kin(i).pro.r_all_pro = r_all_pro; pop.kin(i).pro.kin_all_pro = kin_all_pro;
    %     pop.kin(i).anti.r_all_anti = r_all_anti;  pop.kin(i).anti.kin_all_anti = kin_all_anti;
    %
    %      [units(i).stats.sacc.regress.coeff_pro, units(i).stats.sacc.regress.CI_pro, units(i).stats.sacc.regress.rsq_pro, units(i).stats.sacc.regress.reg_stats_pro] = ...
    %         regress(r_all_pro,kin_all_pro);
    %
    %     [units(i).stats.sacc.regress.coeff_anti, units(i).stats.sacc.regress.CI_anti, units(i).stats.sacc.regress.rsq_anti, units(i).stats.sacc.regress.reg_stats_anti] = ...
    %         regress(r_all_anti,kin_all_anti);
    
    % save all r and eye kinmeatics
    pop.kin(i).pro.r_all_pro = r_all_pro; pop.kin(i).pro.kin_all_pro = kin_all_pro;
    pop.kin(i).anti.r_all_anti = r_all_anti;  pop.kin(i).anti.kin_all_anti = kin_all_anti;
    
    %      [units(i).stats.sacc.regress.coeff_pro, units(i).stats.sacc.regress.CI_pro, units(i).stats.sacc.regress.rsq_pro, units(i).stats.sacc.regress.reg_stats_pro] = ...
    %         regress(r_all_pro,kin_all_pro);
    %
    %     [units(i).stats.sacc.regress.coeff_anti, units(i).stats.sacc.regress.CI_anti, units(i).stats.sacc.regress.rsq_anti, units(i).stats.sacc.regress.reg_stats_anti] = ...
    %         regress(r_all_anti,kin_all_anti);
    
    %% amplitude
    [units(i).stats.sacc.regress.coeff_pro_amp, units(i).stats.sacc.regress.CI_pro_amp, units(i).stats.sacc.regress.rsq_pro_amp, units(i).stats.sacc.regress.reg_stats_pro_amp] = ...
        regress(r_all_pro,[amp_pro' ones(size(amp_pro,2),1)]);
    
    [units(i).stats.sacc.regress.coeff_anti_amp, units(i).stats.sacc.regress.CI_anti_amp, units(i).stats.sacc.regress.rsq_anti_amp, units(i).stats.sacc.regress.reg_stats_anti_amp] = ...
        regress(r_all_anti,[amp_anti' ones(size(amp_anti,2),1)]);
    
    % corrected: firing rate/max coefficient
    pop.stats.sacc.regress(i).coeff_pro_amp_corrected = mean(r_all_pro./units(i).stats.sacc.regress.coeff_pro_amp(1));
    pop.stats.sacc.regress(i).coeff_anti_amp_corrected = mean(r_all_anti./units(i).stats.sacc.regress.coeff_anti_amp(1));
    
    %% peak velocity
    
    [units(i).stats.sacc.regress.coeff_pro_peakVel, units(i).stats.sacc.regress.CI_pro_peakVel, units(i).stats.sacc.regress.rsq_pro_peakVel, units(i).stats.sacc.regress.reg_stats_pro_peakVel] = ...
        regress(r_all_pro,[pv_pro' ones(size(pv_pro,2),1)]);
    
    [units(i).stats.sacc.regress.coeff_anti_peakVel, units(i).stats.sacc.regress.CI_anti_peakVel, units(i).stats.sacc.regress.rsq_anti_peakVel, units(i).stats.sacc.regress.reg_stats_anti_peakVel] = ...
        regress(r_all_anti,[pv_anti' ones(size(pv_anti,2),1)]);
    
    % corrected: firing rate/max coefficient
    pop.stats.sacc.regress(i).coeff_pro_pv_corrected = mean(r_all_pro'./units(i).stats.sacc.regress.coeff_pro_peakVel(1));
    pop.stats.sacc.regress(i).coeff_anti_pv_corrected = mean(r_all_anti'./units(i).stats.sacc.regress.coeff_anti_peakVel(1));
    
    % save firing rate for all per area
    if strcmp(units(i).area, 'lateral')
        pop.stats.sacc.lateral.r_pro_all = [pop.stats.sacc.lateral.r_pro_all ; r_all_pro];
        pop.stats.sacc.lateral.r_anti_all = [pop.stats.sacc.lateral.r_anti_all ; r_all_anti];
    else
        pop.stats.sacc.vermis.r_pro_all = [pop.stats.sacc.vermis.r_pro_all ; r_all_pro];
        pop.stats.sacc.vermis.r_anti_all = [pop.stats.sacc.vermis.r_anti_all ; r_all_anti];
    end
    %% Compute Modulation sensitivity (De Zeeuw et al. 1995)
    % Magnitude sensitivity = sqrt (SS firing rate in pro)^2 + (Pro regress coefficient amplitude)^2 + (Pro regress coefficient duration)^2 + (Pro regress coefficient peak velocity)^2 )
    
    % pop.stats.sacc.pro.mag_sensitivity(i) = sqrt((units(i).stats.sacc.regress.coeff_pro(2)^2) + (units(i).stats.sacc.regress.coeff_pro(3)^2) + (units(i).stats.sacc.regress.coeff_pro(4)^2));
    % pop.stats.sacc.anti.mag_sensitivity(i) = sqrt((units(i).stats.sacc.regress.coeff_anti(2)^2) + (units(i).stats.sacc.regress.coeff_anti(3)^2) + (units(i).stats.sacc.regress.coeff_anti(4)^2));
    
    
    
end


%% Modulation ratio
for cellNum = 1:length(units)
    units(cellNum).stats.mod_ratio_instr = units(cellNum).pro.neural.instr_back.rate_pst_win ./ units(cellNum).anti.neural.instr_back.rate_pst_win;
    units(cellNum).stats.mod_ratio = units(cellNum).pro.neural.sacc.rate_pst_win ./ units(cellNum).anti.neural.sacc.rate_pst_win;
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

%%
pop.indx_sel.vermis = createTable_pro_anti_omv(units); 
pop.indx_sel.lateral = createTable_pro_anti_lat(units); 

%% Compute kstest on pro vs anti on all exc and sup neurons
% r_exc_pro = []; r_sup_pro = [];  r_exc_anti = []; r_sup_anti = [];
% 
% cnt_exc=1; cnt_sup=1;
% for cellNum = 1:length(units)
%     if strcmp(units(cellNum).area, 'vermis') && units(cellNum).pro.neural.exc==1
%         indx_exc(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
%     elseif strcmp(units(cellNum).area,  'vermis') && units(cellNum).pro.neural.sup==1
%         indx_sup(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
%     end
% end
% 
% for i = 1:length(indx_exc_pro), r_exc_pro = [r_exc_pro; units(indx_exc_pro(i)).pro.neural.sacc.nspk]; end
% for i = 1:length(indx_sup_pro), r_sup_pro = [r_sup_pro; units(indx_sup_pro(i)).pro.neural.sacc.nspk]; end
% for i = 1:length(indx_exc_anti), r_exc_anti = [r_exc_anti; units(indx_exc_anti(i)).anti.neural.sacc.nspk]; end
% for i = 1:length(indx_sup_anti), r_sup_anti = [r_sup_anti; units(indx_sup_anti(i)).anti.neural.sacc.nspk]; end
% 
% [pop.stats.sacc.vermis.pro_anti_exc.flag, pop.stats.sacc.vermis.pro_anti_exc.pVal] = kstest2(r_exc_pro, r_exc_anti); 
% [pop.stats.sacc.vermis.pro_anti_sup.flag, pop.stats.sacc.vermis.pro_anti_sup.pVal] = kstest2(r_sup_pro, r_sup_anti); 
% 
% % lateral
% cnt_exc=1; cnt_sup=1;
% for cellNum = 1:length(units)
%     if strcmp(units(cellNum).area, 'lateral') && units(cellNum).pro.neural.exc==1
%         indx_exc(cnt_exc) = cellNum; cnt_exc=cnt_exc+1;
%     elseif strcmp(units(cellNum).area,  'lateral') && units(cellNum).pro.neural.sup==1
%         indx_sup(cnt_sup) = cellNum; cnt_sup=cnt_sup+1;
%     end
% end
% 
% for i = 1:length(indx_exc_pro), r_exc_pro = [r_exc_pro; units(indx_exc_pro(i)).pro.neural.sacc.nspk]; end
% for i = 1:length(indx_sup_pro), r_sup_pro = [r_sup_pro; units(indx_sup_pro(i)).pro.neural.sacc.nspk]; end
% for i = 1:length(indx_exc_anti), r_exc_anti = [r_exc_anti; units(indx_exc_anti(i)).anti.neural.sacc.nspk]; end
% for i = 1:length(indx_sup_anti), r_sup_anti = [r_sup_anti; units(indx_sup_anti(i)).anti.neural.sacc.nspk]; end
% 
% [pop.stats.sacc.lateral.pro_anti_exc.flag, pop.stats.sacc.lateral.pro_anti_exc.pVal] = kstest2(r_exc_pro, r_exc_anti); 
% [pop.stats.sacc.lateral.pro_anti_sup.flag, pop.stats.sacc.lateral.pro_anti_sup.pVal] = kstest2(r_sup_pro, r_sup_anti); 

%% Compute kstest on pro vs anti on selected neurons
r_exc_pro = []; r_sup_pro = [];  r_exc_anti = []; r_sup_anti = [];

% vermis
indx_exc_pro = pop.indx_sel.vermis.sacc.all.pro.exc;
indx_sup_pro = pop.indx_sel.vermis.sacc.all.pro.sup;
indx_exc_anti = pop.indx_sel.vermis.sacc.all.anti.exc;
indx_sup_anti = pop.indx_sel.vermis.sacc.all.anti.sup;

for i = 1:length(indx_exc_pro), r_exc_pro = [r_exc_pro; units(indx_exc_pro(i)).pro.neural.sacc.nspk]; end
for i = 1:length(indx_sup_pro), r_sup_pro = [r_sup_pro; units(indx_sup_pro(i)).pro.neural.sacc.nspk]; end
for i = 1:length(indx_exc_anti), r_exc_anti = [r_exc_anti; units(indx_exc_anti(i)).anti.neural.sacc.nspk]; end
for i = 1:length(indx_sup_anti), r_sup_anti = [r_sup_anti; units(indx_sup_anti(i)).anti.neural.sacc.nspk]; end

[pop.stats.sacc.vermis.pro_anti_exc_sel.flag, pop.stats.sacc.vermis.pro_anti_exc_sel.pVal] = kstest2(r_exc_pro, r_exc_anti); 
[pop.stats.sacc.vermis.pro_anti_sup_sel.flag, pop.stats.sacc.vermis.pro_anti_sup_sel.pVal] = kstest2(r_sup_pro, r_sup_anti); 

% lateral
r_exc_pro = []; r_sup_pro = [];  r_exc_anti = []; r_sup_anti = [];
indx_exc_pro = pop.indx_sel.lateral.sacc.all.pro.exc;
indx_sup_pro = pop.indx_sel.lateral.sacc.all.pro.sup;
indx_exc_anti = pop.indx_sel.lateral.sacc.all.anti.exc;
indx_sup_anti = pop.indx_sel.lateral.sacc.all.anti.sup;

for i = 1:length(indx_exc_pro), r_exc_pro = [r_exc_pro; units(indx_exc_pro(i)).pro.neural.sacc.nspk]; end
for i = 1:length(indx_sup_pro), r_sup_pro = [r_sup_pro; units(indx_sup_pro(i)).pro.neural.sacc.nspk]; end
for i = 1:length(indx_exc_anti), r_exc_anti = [r_exc_anti; units(indx_exc_anti(i)).anti.neural.sacc.nspk]; end
for i = 1:length(indx_sup_anti), r_sup_anti = [r_sup_anti; units(indx_sup_anti(i)).anti.neural.sacc.nspk]; end

[pop.stats.sacc.lateral.pro_anti_exc_sel.flag, pop.stats.sacc.lateral.pro_anti_exc_sel.pVal] = kstest2(r_exc_pro, r_exc_anti); 
[pop.stats.sacc.lateral.pro_anti_sup_sel.flag, pop.stats.sacc.lateral.pro_anti_sup_sel.pVal] = kstest2(r_sup_pro, r_sup_anti); 


