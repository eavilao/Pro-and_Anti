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
%uiopen;
tsmooth = 0.050;

% first find the ones that are not empty
for i = 1:length(wholeNeuronResults)
    if  ~isempty(wholeNeuronResults(i).allStableTrials);
      cell_indx(i) = numel(wholeNeuronResults(i).selectedTrials.corProTrials)>=5 & numel(wholeNeuronResults(i).selectedTrials.corAntiTrials)>=5;  
    end
end
cell_indx = find(cell_indx);

for cellNum = 1:length(cell_indx);
    for trialNum = 1:length(wholeNeuronResults(cell_indx(cellNum)).allStableTrials);
        if ~isempty(wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).alignedSpikes)
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
            
            %neural data
            spks =  wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).alignedSpikes{1}; % contains spike times for SS aligned to trial onset
            if ~isempty(units(cellNum).trial.behav(trialNum).reward)
                indx_trial_spks = spks>-0.1 & spks < units(cellNum).trial.behav(trialNum).reward+0.2; % just pick 100 ms before trial starts to reward +200 ms
                units(cellNum).trial.neural(trialNum).tspk_SS = spks(indx_trial_spks);
                units(cellNum).trial.neural(trialNum).tspk_SS_align_sacc =  units(cellNum).trial.neural(trialNum).tspk_SS-units(cellNum).trial.behav(trialNum).saccadeOnset; % contains spike times for CS aligned to sacc onset
            else
                units(cellNum).trial.neural(trialNum).tspk_SS = [];
                units(cellNum).trial.neural(trialNum).tspk_SS_align_sacc = [];
            end
            units(cellNum).id = 'SS';
            % select condition type Pro or Antisaccade
            if ismember(wholeNeuronResults(cell_indx(cellNum)).allStableTrials(trialNum).conditionCode, [1 4 5 8 10 11 14 15]);
                units(cellNum).trial.behav(trialNum).condition = 'Prosaccade';
            else
                units(cellNum).trial.behav(trialNum).condition = 'Antisaccade';
            end
        end
        
        % Indexes to select Pro and Anti trials
        units(cellNum).pro.indx_correctProTrials = wholeNeuronResults(cell_indx(cellNum)).selectedTrials.corProTrials;
        units(cellNum).anti.indx_correctAntiTrials = wholeNeuronResults(cell_indx(cellNum)).selectedTrials.corAntiTrials;
        
    end
end

clear cellNum trialNum %wholeNeuronResults

%% Spike times to rate for pro and antisaccades and behav
binwidth = 0.01;
timepoints_instr = -0.1:binwidth:1;
timepoints_sacc = -0.8:binwidth:0.3;
analyse_sacc_align = 0;

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
    end
    
    % neural
    % instr
    units(cellNum).pro.neural.trial = units(cellNum).trial.neural(units(cellNum).pro.indx_correctProTrials);
    %instr
    [units(cellNum).pro.neural.instr.rate_pst,units(cellNum).pro.neural.instr.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,timepoints_instr,binwidth,analyse_sacc_align,id); % aligned to trial onset
    units(cellNum).pro.neural.instr.rate_pst = smooth_pst(units(cellNum).pro.neural.instr.rate_pst,binwidth,tsmooth);
    % sacc aligned
    analyse_sacc_align=1;
    [units(cellNum).pro.neural.sacc.rate_pst,units(cellNum).pro.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,timepoints_sacc,binwidth,analyse_sacc_align,id); % aligned to saccade onset
    units(cellNum).pro.neural.sacc.rate_pst = smooth_pst(units(cellNum).pro.neural.sacc.rate_pst,binwidth,tsmooth);
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
    end
    %neural
    units(cellNum).anti.neural.trial = units(cellNum).trial.neural(units(cellNum).anti.indx_correctAntiTrials);
    % intr
    [units(cellNum).anti.neural.instr.rate_pst,units(cellNum).anti.neural.instr.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,timepoints_instr,binwidth,analyse_sacc_align,id); % aligned to trial onset
    units(cellNum).anti.neural.instr.rate_pst = smooth_pst(units(cellNum).anti.neural.instr.rate_pst,binwidth,tsmooth);
    % sacc aligned
    analyse_sacc_align=1;
    [units(cellNum).anti.neural.sacc.rate_pst,units(cellNum).anti.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,timepoints_sacc,binwidth,analyse_sacc_align,id); % aligned to saccade onset
    units(cellNum).anti.neural.sacc.rate_pst = smooth_pst(units(cellNum).anti.neural.sacc.rate_pst,binwidth,tsmooth);
    analyse_sacc_align=0;
end

%% For every trial - Spike time count and spike times to rate for pro and antisaccades

for cellNum = 1:length(units)
    id=units(cellNum).id;
    %pro trials
    correctProTrials = units(cellNum).pro.indx_correctProTrials; % index pro trials
    % get spk counts
    for trialNum = 1:length(correctProTrials)
        timepoints_instr = -0.1:binwidth:1; timepoints_sacc = -0.8:binwidth:0.3;
        [units(cellNum).pro.neural.instr.spkCount(trialNum,:),units(cellNum).pro.neural.instr.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).pro.neural.trial(trialNum),timepoints_instr,binwidth,analyse_sacc_align,id); % spk counts
        analyse_sacc_align=1;
        [units(cellNum).pro.neural.sacc.spkCount(trialNum,:),units(cellNum).pro.neural.sacc.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).pro.neural.trial(trialNum),timepoints_sacc,binwidth,analyse_sacc_align,id);
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
        [units(cellNum).anti.neural.instr.spkCount(trialNum,:),units(cellNum).anti.neural.instr.trial(trialNum).ts] = Spiketimes2CountTrial(units(cellNum).anti.neural.trial(trialNum),timepoints_instr,binwidth,analyse_sacc_align,id); % spk counts
        analyse_sacc_align=1;
        [units(cellNum).anti.neural.sacc.spkCount(trialNum,:),units(cellNum).anti.neural.sacc.ts_spkCount(trialNum,:)] = Spiketimes2CountTrial(units(cellNum).anti.neural.trial(trialNum),timepoints_sacc,binwidth,analyse_sacc_align,id);
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
    for trialNum = 1:length(correctProTrials)
        timepoints_instr = -0.1:binwidth:1; timepoints_sacc = -0.8:binwidth:0.3;
        % align to instr
        [~,units(cellNum).pro.neural.instr.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum).tspk_SS,timepoints_instr,binwidth,analyse_sacc_align,id);
        
        % win
        tspk_instr= units(cellNum).pro.neural.trial(trialNum).tspk_SS>0 &units(cellNum).pro.neural.trial(trialNum).tspk_SS<0.301;
        t_instr = timepoints_instr > 0 & timepoints_instr<0.301; timepoints_instr = timepoints_instr(t_instr);
        units(cellNum).pro.neural.instr.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS(tspk_instr);
        [units(cellNum).pro.neural.instr.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.instr.tspk),timepoints_instr,binwidth,analyse_sacc_align,id);
        
        % align to sacc
        analyse_sacc_align = 1;
        [~,units(cellNum).pro.neural.sacc.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc,timepoints_sacc,binwidth,analyse_sacc_align,id);
        
        % win
        tspk_sacc= units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc>-0.1 & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc<0.201;
        t_sacc = timepoints_sacc > -0.1 & timepoints_sacc<0.2; timepoints_sacc = timepoints_sacc(t_sacc);
        units(cellNum).pro.neural.sacc.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc(tspk_sacc);
        [units(cellNum).pro.neural.sacc.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.sacc.tspk),timepoints_sacc,binwidth,analyse_sacc_align,id);
        analyse_sacc_align = 0;
        
        % baseline
        tspk_base = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc>-0.3 & units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc<-0.101;
        t_base = timepoints_sacc >-0.3 & timepoints_sacc<-0.1; timepoints_base = timepoints_sacc(t_base);
        units(cellNum).pro.neural.base.tspk{trialNum,:} = units(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc(tspk_sacc);
        [units(cellNum).pro.neural.base.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).pro.neural.base.tspk),timepoints_sacc,binwidth,analyse_sacc_align,id);
        analyse_sacc_align = 0;
    end
    
    % nspk anti - windows per trial
    for trialNum = 1:length(correctAntiTrials)
        timepoints_instr = -0.1:binwidth:1; timepoints_sacc = -0.8:binwidth:0.3;
        % align to instr
        [~,units(cellNum).anti.neural.instr.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum).tspk_SS,timepoints_instr,binwidth,analyse_sacc_align,id);
        
        % win
        tspk_instr= units(cellNum).anti.neural.trial(trialNum).tspk_SS>0 &units(cellNum).anti.neural.trial(trialNum).tspk_SS<0.301;
        t_instr = timepoints_instr > 0 & timepoints_instr<0.301; timepoints_instr = timepoints_instr(t_instr);
        units(cellNum).anti.neural.instr.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS(tspk_instr);
        [units(cellNum).anti.neural.instr.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.instr.tspk),timepoints_instr,binwidth,analyse_sacc_align,id);
        
        % align to sacc
        analyse_sacc_align = 1;
        [~,units(cellNum).anti.neural.sacc.nspkCount(trialNum,:),~] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc,timepoints_sacc,binwidth,analyse_sacc_align,id);
        
        % win
        tspk_sacc= units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc>-0.1 &units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc<0.2;
        t_sacc = timepoints_sacc > -0.101 & timepoints_sacc<0.2; timepoints_sacc = timepoints_sacc(t_sacc);
        units(cellNum).anti.neural.sacc.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc(tspk_sacc);
        [units(cellNum).anti.neural.sacc.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.sacc.tspk),timepoints_sacc,binwidth,analyse_sacc_align,id);
        analyse_sacc_align = 0;
        
        % baseline
        tspk_base = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc>-0.3 & units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc<-0.1;
        t_base = timepoints_sacc >-0.3 & timepoints_sacc<-0.1; timepoints_base = timepoints_sacc(t_base);
        units(cellNum).anti.neural.base.tspk{trialNum,:} = units(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc(tspk_sacc);
        [units(cellNum).anti.neural.base.nspk(trialNum,:),~,~] = Spiketimes2RateTrial(cell2mat(units(cellNum).anti.neural.base.tspk),timepoints_sacc,binwidth,analyse_sacc_align,id);
        analyse_sacc_align = 0;
        
    end
    
    %% probability of spk in pro and anti  - bigger window (5 bins = 50 ms)
    
    win_size = 5; % num of bins to take in window. 1 bin = 10 ms
    indx_win = win_size;
    % instr pro
    for trialNum = 1:length(correctProTrials)
        thisTrial = units(cellNum).pro.neural.instr.nspkCount(trialNum,:);
        indx_beg = 1;indx_win_count = win_size;
        for win_num = 1:length(thisTrial(1,:))-indx_win
            units(cellNum).pro.neural.instr.spkCount_win(trialNum,win_num) = sum(thisTrial(indx_beg:indx_win_count))/win_size;
            indx_beg = indx_beg+1; indx_win_count = indx_win_count+1;
        end
    end
    
    % instr anti
    for trialNum = 1:length(correctAntiTrials)
        thisTrial = units(cellNum).anti.neural.instr.nspkCount(trialNum,:);
        indx_beg = 1;indx_win_count = win_size;
        for win_num = 1:length(thisTrial(1,:))-indx_win
            units(cellNum).anti.neural.instr.spkCount_win(trialNum,win_num) = sum(thisTrial(indx_beg:indx_win_count))/win_size;
            indx_beg = indx_beg+1; indx_win_count = indx_win_count+1;
        end
    end
    
    % sacc pro
    for trialNum = 1:length(correctProTrials)
        thisTrial = units(cellNum).pro.neural.sacc.nspkCount(trialNum,:);
        indx_beg = 1;indx_win_count = win_size;
        for win_num = 1:length(thisTrial(1,:))-indx_win
            units(cellNum).pro.neural.sacc.spkCount_win(trialNum,win_num) = sum(thisTrial(indx_beg:indx_win_count))/win_size;
            indx_beg = indx_beg+1; indx_win_count = indx_win_count+1;
        end
    end
    
    % sacc anti
    for trialNum = 1:length(correctAntiTrials)
        thisTrial = units(cellNum).anti.neural.sacc.nspkCount(trialNum,:);
        indx_beg = 1;indx_win_count = win_size;
        for win_num = 1:length(thisTrial(1,:))-indx_win
            units(cellNum).anti.neural.sacc.spkCount_win(trialNum,win_num) = sum(thisTrial(indx_beg:indx_win_count))/win_size;
            indx_beg = indx_beg+1; indx_win_count = indx_win_count+1;
        end
    end
    
end

%% STATS %%
%%%%%%%%%%%
%%%%%%%%%%%
%% Eye kinematics
eyeKin = eyeKinematics_ProAnti(units); % extract and plot

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
    instr_win = units(cellNum).pro.neural.instr.ts_pst > 0 & units(cellNum).pro.neural.instr.ts_pst < 0.301;
    sacc_win = units(cellNum).pro.neural.sacc.ts_pst > -0.1 & units(cellNum).pro.neural.sacc.ts_pst < 0.201;
    base_win = units(cellNum).pro.neural.sacc.ts_pst > -0.301 & units(cellNum).pro.neural.sacc.ts_pst <-0.101;
    %get spks pro
    instr_spks_pro = units(cellNum).pro.neural.instr.rate_pst(instr_win);
    sacc_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(sacc_win);
    base_spks_pro = units(cellNum).pro.neural.sacc.rate_pst(base_win);
    
    %pst, mean and sem
    units(cellNum).pro.neural.instr.ts_pst_win = units(cellNum).pro.neural.instr.ts_pst(instr_win);
    units(cellNum).pro.neural.instr.rate_pst_win = units(cellNum).pro.neural.instr.rate_pst(instr_win);
    units(cellNum).pro.neural.instr.rate_mu = mean(instr_spks_pro);
    units(cellNum).pro.neural.instr.rate_sig = std(instr_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.sacc.ts_pst_win = units(cellNum).pro.neural.sacc.ts_pst(sacc_win);
    units(cellNum).pro.neural.sacc.rate_pst_win = units(cellNum).pro.neural.sacc.rate_pst(sacc_win);
    units(cellNum).pro.neural.sacc.rate_mu = mean(sacc_spks_pro);
    units(cellNum).pro.neural.sacc.rate_sig = std(sacc_spks_pro)/sqrt(ntrls_pro);
    units(cellNum).pro.neural.base.ts_pst_win = units(cellNum).pro.neural.sacc.ts_pst(base_win);
    units(cellNum).pro.neural.base.rate_pst_win = units(cellNum).pro.neural.sacc.rate_pst(base_win);
    units(cellNum).pro.neural.base.rate_mu = mean(base_spks_pro);
    units(cellNum).pro.neural.base.rate_sig = std(base_spks_pro)/sqrt(ntrls_pro);
    
    
    %get spks anti
    instr_spks_anti = units(cellNum).anti.neural.instr.rate_pst(instr_win);
    sacc_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(sacc_win);
    base_spks_anti = units(cellNum).anti.neural.sacc.rate_pst(base_win);
    %pst, mean and sem
    units(cellNum).anti.neural.instr.ts_pst_win = units(cellNum).anti.neural.instr.ts_pst(instr_win);
    units(cellNum).anti.neural.instr.rate_pst_win = units(cellNum).anti.neural.instr.rate_pst(instr_win);
    units(cellNum).anti.neural.instr.rate_mu = mean(instr_spks_anti);
    units(cellNum).anti.neural.instr.rate_sig = std(instr_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.sacc.ts_pst_win = units(cellNum).anti.neural.sacc.ts_pst(sacc_win);
    units(cellNum).anti.neural.sacc.rate_pst_win = units(cellNum).anti.neural.sacc.rate_pst(sacc_win);
    units(cellNum).anti.neural.sacc.rate_mu = mean(sacc_spks_anti);
    units(cellNum).anti.neural.sacc.rate_sig = std(sacc_spks_anti)/sqrt(ntrls_anti);
    units(cellNum).anti.neural.base.ts_pst_win = units(cellNum).anti.neural.sacc.ts_pst(base_win);
    units(cellNum).anti.neural.base.rate_pst_win = units(cellNum).anti.neural.sacc.rate_pst(base_win);
    units(cellNum).anti.neural.base.rate_mu = mean(base_spks_anti);
    units(cellNum).anti.neural.base.rate_sig = std(base_spks_anti)/sqrt(ntrls_anti);
    
    %% compare windows against baseline activity for pro and anti
    % pro
    [units(cellNum).stats.pro.pval.instrVSbase, units(cellNum).stats.pro.flags.instrVSbase] = signrank(units(cellNum).pro.neural.instr.nspk,units(cellNum).pro.neural.base.nspk);
    [units(cellNum).stats.pro.pval.saccVSbase, units(cellNum).stats.pro.flags.saccVSbase] = signrank(units(cellNum).pro.neural.sacc.nspk,units(cellNum).pro.neural.base.nspk);
    
    % anti
    [units(cellNum).stats.anti.pval.instrVSbase, units(cellNum).stats.anti.flags.instrVSbase] = signrank(units(cellNum).anti.neural.instr.nspk,units(cellNum).anti.neural.base.nspk);
    [units(cellNum).stats.anti.pval.saccVSbase, units(cellNum).stats.anti.flags.saccVSbase] = signrank(units(cellNum).anti.neural.sacc.nspk,units(cellNum).anti.neural.base.nspk);
    
    %% compare pro vs anti

    [units(cellNum).stats.instr.flags.proVsAnti_instr,units(cellNum).stats.instr.pval.proVsAnti_instr_win] = ttest2(units(cellNum).pro.neural.instr.nspk,units(cellNum).anti.neural.instr.nspk);
    [units(cellNum).stats.sacc.flags.proVsAnti_sacc,units(cellNum).stats.sacc.pval.proVsAnti_sacc_win] = ttest2(units(cellNum).pro.neural.sacc.nspk,units(cellNum).anti.neural.sacc.nspk);
    
    %% Compute change in FR from baseline
    % pro
    units(cellNum).pro.neural.instr.delta_rate = units(cellNum).pro.neural.instr.rate_pst_win - units(cellNum).pro.neural.instr.rate_mu;
    units(cellNum).pro.neural.sacc.delta_rate = units(cellNum).pro.neural.sacc.rate_pst_win - units(cellNum).pro.neural.sacc.rate_mu;
    
    % anti
    units(cellNum).anti.neural.instr.delta_rate = units(cellNum).anti.neural.instr.rate_pst_win - units(cellNum).anti.neural.instr.rate_mu;
    units(cellNum).anti.neural.sacc.delta_rate = units(cellNum).anti.neural.sacc.rate_pst_win - units(cellNum).anti.neural.sacc.rate_mu;
    
    %% noramlized FR - z-scored
    % pro
    % aligned to trial onset
    units(cellNum).pro.neural.sacc.norm_rate_pst = (units(cellNum).pro.neural.instr.rate_pst - mean(units(cellNum).pro.neural.instr.rate_pst))/...
        std(units(cellNum).pro.neural.instr.rate_pst);
    
    % aligned to saccade
    units(cellNum).pro.neural.sacc.norm_rate_pst = (units(cellNum).pro.neural.sacc.rate_pst - mean(units(cellNum).pro.neural.sacc.rate_pst))/...
        std(units(cellNum).pro.neural.sacc.rate_pst);
    
    %anti
    % aligned to trial onset
    units(cellNum).anti.neural.sacc.norm_rate_pst = (units(cellNum).anti.neural.instr.rate_pst - mean(units(cellNum).anti.neural.instr.rate_pst))/...
        std(units(cellNum).anti.neural.instr.rate_pst);
    
    % aligned to saccade
    units(cellNum).anti.neural.sacc.norm_rate_pst = (units(cellNum).anti.neural.sacc.rate_pst - mean(units(cellNum).anti.neural.sacc.rate_pst))/...
        std(units(cellNum).anti.neural.sacc.rate_pst);
    
    
    %% statistical test to compare if pro == anti - aligned to instr using spk count
    % H0:pro=anti  versus HA:pro?anti
    correctProTrials = units(cellNum).pro.indx_correctProTrials; % index pro trials
    correctAntiTrials = units(cellNum).anti.indx_correctAntiTrials;
    ntrls_pro = length(correctProTrials);
    ntrls_anti = length(correctAntiTrials);
    signif_criteria = 1.96; %two tailed test
    
    % instr
    for i=1:length(units(cellNum).pro.neural.instr.spkCount(1,:))
        p=(ntrls_pro*units(cellNum).pro.neural.instr.pbDist(i)+ntrls_anti*units(cellNum).anti.neural.instr.pbDist(i))/(ntrls_pro+ntrls_anti);
        units(cellNum).stats.instr.pval.pbDist_testStat(i) = (units(cellNum).pro.neural.instr.pbDist(1,i)-units(cellNum).anti.neural.instr.pbDist(1,i))/(sqrt(p*(1-p)*(1/ntrls_pro + 1/ntrls_anti)));
        if units(cellNum).stats.instr.pval.pbDist_testStat(i) > signif_criteria | units(cellNum).stats.instr.pval.pbDist_testStat(i) < -signif_criteria
            units(cellNum).stats.instr.flags.pbDist(i) = 1;
        else
            units(cellNum).stats.instr.flags.pbDist(i) = 0;
        end
    end
    
    % sacc
    for i=1:length(units(cellNum).anti.neural.sacc.spkCount(1,:))
        p=(ntrls_pro*units(cellNum).pro.neural.sacc.pbDist(i)+ntrls_anti*units(cellNum).anti.neural.sacc.pbDist(i))/(ntrls_pro+ntrls_anti);
        units(cellNum).stats.sacc.pval.pbDist_testStat(i) = (units(cellNum).pro.neural.sacc.pbDist(1,i)-units(cellNum).anti.neural.sacc.pbDist(1,i))/(sqrt(p*(1-p)*(1/ntrls_pro + 1/ntrls_anti)));
        if units(cellNum).stats.sacc.pval.pbDist_testStat(i) > signif_criteria | units(cellNum).stats.sacc.pval.pbDist_testStat(i) < -signif_criteria
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
    
    
    %% Sliding window to test diff pro vs anti
    % detect saccade related using sliding window (5 bins) 200 ms before
    % sacc onset to 200 ms after sacc onset and store time at > or < than 2*std
    
    win_size = 5; % num of bins to take in sliding window. 1 bin = 10 ms
    indx_win = win_size; win_size_prev = 10;
    t = units(cellNum).pro.neural.sacc.ts_pst;
    indx_beg = find(t>-0.2,1); indx_end = find(t>0.2,1);
    
    % pro
    rate_pst = units(cellNum).pro.neural.sacc.rate_pst;
    exc_pro=0; sup_pro=0;exc_anti=0; sup_anti=0;
    
    for indx = indx_beg:indx_end
        r_mu = mean(rate_pst(round(indx-win_size/2):round(indx+win_size/2)));
        r_thresh1 = mean(rate_pst(indx-win_size_prev:indx))+...
            2*std(rate_pst(indx-win_size_prev:indx));
        r_thresh2 = mean(rate_pst(indx-win_size_prev:indx))-...
            2*std(rate_pst(indx-win_size_prev:indx));
        if (r_thresh1>1 && r_mu>r_thresh1), exc_pro=exc_pro+1;else exc_pro=0;end
        if (r_thresh2>1 && r_mu<r_thresh2), sup_pro=sup_pro+1;else sup_pro=0;end
        if exc_pro==3, events.rise.t_on = t(indx); events.rise.type_on = 'exc';
            break;
        elseif sup_pro==3, units(cellNum).pro.events.rise.t_on = t(indx); units(cellNum).pro.events.rise.type_on = 'sup';
            break;
        end
    end
    
    % anti
    rate_pst = units(cellNum).anti.neural.sacc.rate_pst;
    exc_anti=0; sup_anti=0;
    
    for indx = indx_beg:indx_end
        r_mu = mean(rate_pst(round(indx-win_size/2):round(indx+win_size/2)));
        r_thresh1 = mean(rate_pst(indx-win_size_prev:indx))+...
            2*std(rate_pst(indx-win_size_prev:indx));
        r_thresh2 = mean(rate_pst(indx-win_size_prev:indx))-...
            2*std(rate_pst(indx-win_size_prev:indx));
        if (r_thresh1>1 && r_mu>r_thresh1), exc_anti=exc_anti+1;else exc_anti=0;end
        if (r_thresh2>1 && r_mu<r_thresh2), sup_anti=sup_anti+1;else sup_anti=0;end
        if exc_anti==3, events.rise.t_on = t(indx); events.rise.type_on = 'exc';
            break;
        elseif sup_anti==3, units(cellNum).pro.events.rise.t_on = t(indx); units(cellNum).pro.events.rise.type_on = 'sup';
            break;
        end
    end
    
end

%% coefficient of variation and CV2
% CV = mean_isi/std_isi
% CV2 =  2x|ISIn+1-ISIn|/(ISIn+1 + ISIn) Holt et al. 199


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
    units(cellNum).pro.meanCV =  mean([units(cellNum).pro.neural.trial.cv_isi]);
    
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
    
    units(cellNum).anti.meanCV =  mean([units(cellNum).anti.neural.trial.cv_isi]);
    
end

end
