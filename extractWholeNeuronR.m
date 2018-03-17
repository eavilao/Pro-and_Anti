function units = extractWholeNeuronR(wholeNeuronResults)
% extractWholeNeuronResults
% This function extracts data from wholeNeuronResults from Pete's code and
% performs stats on simple spikes.

% In: wholeNeuronResults.m
% Out:  units containing spike times (Simple and Complex spikes). Every
% struct inside units is a neuron.
% aligned to trial onset and eye stats.

% Associated functions: SpikeTimes2Rate.m, Spiketimes2RateTrial.m,
%                       plot_ProAnti.m, eyeKinematics_ProAnti.m

% TODO
% - eye kin plots to separate code

%% Create structure with relevant data (units)
%uiopen;
tsmooth = 0.050;

for cellNum = 1:length(wholeNeuronResults);
    for trialNum = 1:length(wholeNeuronResults(cellNum).allStableTrials);
        if ~isempty(wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes)
            
            %behavioral
            units(cellNum).trial.behav(trialNum).saccadeOnset =  wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeTime; % saccade onset from trial start (?)
            units(cellNum).trial.behav(trialNum).goCueTime = wholeNeuronResults(cellNum).allStableTrials(trialNum).goCueTime; %go cue from trial start
            units(cellNum).trial.behav(trialNum).saccAmplitude = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeAmplitude;
            units(cellNum).trial.behav(trialNum).saccDuration = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeDuration;
            units(cellNum).trial.behav(trialNum).saccPeakVel = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadePeakVel;
            units(cellNum).trial.behav(trialNum).reactionTime = wholeNeuronResults(cellNum).allStableTrials(trialNum).reactionTime;
            units(cellNum).trial.behav(trialNum).reward = wholeNeuronResults(cellNum).allStableTrials(trialNum).relativeRewardTime;
            
            %neural data      
            spks =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{1}; % contains spike times for SS aligned to trial onset
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
            if ismember(wholeNeuronResults(cellNum).allStableTrials(trialNum).conditionCode, [1 4 5 8 10 11 14 15]);
                units(cellNum).trial.behav(trialNum).condition = 'Prosaccade';
            else
                units(cellNum).trial.behav(trialNum).condition = 'Antisaccade';
            end
        end
    end
    % Indexes to select Pro and Anti trials
    units(cellNum).pro.indx_correctProTrials = wholeNeuronResults(cellNum).selectedTrials.corProTrials;
    units(cellNum).anti.indx_correctAntiTrials = wholeNeuronResults(cellNum).selectedTrials.corAntiTrials;
end

clear cellNum trialNum %wholeNeuronResults

%% Spike times to rate for pro and antisaccades and behav
binwidth = 0.01;
timepoints_instr = -0.1:binwidth:1;
timepoints_sacc = -0.8:binwidth:0.3;
analyse_sacc_win = 0;

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
    [units(cellNum).pro.neural.rate_pst,units(cellNum).pro.neural.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,timepoints_instr,binwidth,analyse_sacc_win,id); % aligned to trial onset
    units(cellNum).pro.neural.rate_pst = smooth_pst(units(cellNum).pro.neural.rate_pst,binwidth,tsmooth);
    % sacc aligned
    [units(cellNum).pro.neural.sacc.rate_pst,units(cellNum).pro.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).pro.neural.trial,timepoints_sacc,binwidth,analyse_sacc_win,id); % aligned to saccade onset
    units(cellNum).pro.neural.sacc.rate_pst = smooth_pst(units(cellNum).pro.neural.sacc.rate_pst,binwidth,tsmooth);
    
    
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
    [units(cellNum).anti.neural.rate_pst,units(cellNum).anti.neural.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,timepoints_instr,binwidth,analyse_sacc_win,id); % aligned to trial onset
    units(cellNum).anti.neural.rate_pst = smooth_pst(units(cellNum).anti.neural.rate_pst,binwidth,tsmooth);
    % sacc aligned
    [units(cellNum).anti.neural.sacc.rate_pst,units(cellNum).anti.neural.sacc.ts_pst] = Spiketimes2Rate(units(cellNum).anti.neural.trial,timepoints_sacc,binwidth,analyse_sacc_win,id); % aligned to saccade onset
    units(cellNum).anti.neural.sacc.rate_pst = smooth_pst(units(cellNum).anti.neural.sacc.rate_pst,binwidth,tsmooth);
end

%% For every trial - Spike times to rate for pro and antisaccades and behav

for cellNum = 1:length(units)
    id=units(cellNum).id;
    %pro trials
    correctProTrials = units(cellNum).pro.indx_correctProTrials; % index pro trials
    for trialNum = 1:length(correctProTrials)
        [units(cellNum).pro.neural.spkCount(trialNum,:),units(cellNum).pro.neural.ts_spkCount(trialNum,:)] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum),timepoints_instr,binwidth,analyse_sacc_win,id); % spk counts
        analyse_sacc_win=1; 
        [units(cellNum).pro.neural.sacc.spkCount(trialNum,:),units(cellNum).pro.neural.sacc.ts_spkCount(trialNum,:)] = Spiketimes2RateTrial(units(cellNum).pro.neural.trial(trialNum),timepoints_sacc,binwidth,analyse_sacc_win,id);
        analyse_sacc_win=0;
    end
    
    for i=1:length(units(cellNum).pro.neural.spkCount(1,:))
        units(cellNum).pro.neural.pbDist(1,i)= sum(units(cellNum).pro.neural.spkCount(:,i))/length(correctProTrials); % compute probability of spk
        units(cellNum).pro.neural.pbDist_sem = sqrt(units(cellNum).pro.neural.pbDist(1,i)*(1-units(cellNum).pro.neural.pbDist(1,i))/length(correctProTrials));
    end
    
    % anti trials
    correctAntiTrials = units(cellNum).anti.indx_correctAntiTrials;
    for trialNum = 1:length(correctAntiTrials)
        [units(cellNum).anti.neural.spkCount(trialNum,:),units(cellNum).anti.neural.trial(trialNum).ts] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum),timepoints_instr,binwidth,analyse_sacc_win,id); % spk counts
        analyse_sacc_win=1;
        [units(cellNum).anti.neural.sacc.spkCount(trialNum,:),units(cellNum).anti.neural.sacc.ts_spkCount(trialNum,:)] = Spiketimes2RateTrial(units(cellNum).anti.neural.trial(trialNum),timepoints_sacc,binwidth,analyse_sacc_win,id);
        analyse_sacc_win=0;
        
    end
    for i=1:length(units(cellNum).anti.neural.spkCount(1,:))
        units(cellNum).anti.neural.pbDist(1,i)= sum(units(cellNum).anti.neural.spkCount(:,i))/length(correctAntiTrials);% compute probability of spk
        units(cellNum).anti.neural.pbDist_sem = sqrt(units(cellNum).pro.neural.pbDist(1,i)*(1-units(cellNum).pro.neural.pbDist(1,i))/length(correctProTrials));
    end
    
    %% statistical test to compare if pro == anti - aligned to instr
    % H0:pro=anti  versus HA:pro?anti
    ntrls_pro = length(correctProTrials);
    ntrls_anti = length(correctAntiTrials);
    signif_criteria = 1.96; %two tailed test
    
    % instr
    for i=1:length(units(cellNum).anti.neural.spkCount(1,:))
        p=(ntrls_pro*units(cellNum).pro.neural.pbDist(i)+ntrls_anti*units(cellNum).anti.neural.pbDist(i))/(ntrls_pro+ntrls_anti);
        units(cellNum).stats.instr.pval.pbDist_testStat(i) = (units(cellNum).pro.neural.pbDist(1,i)-units(cellNum).anti.neural.pbDist(1,i))/(sqrt(p*(1-p)*(1/ntrls_pro + 1/ntrls_anti)));
        if units(cellNum).stats.instr.pval.pbDist_testStat(i) > signif_criteria | units(cellNum).stats.pval.pbDist_testStat(i) < -signif_criteria
            units(cellNum).stats.instr.flags.pbDist(i) = 1;
        else
            units(cellNum).stats.instr.flags.pbDist(i) = 0;
        end
    end
    
     % sacc
    for i=1:length(units(cellNum).anti.neural.spkCount(1,:))
        
        
    end 
    
end




