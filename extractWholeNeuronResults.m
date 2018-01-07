% extractWholeNeuronResults
% This function extracts data from wholeNeuronResults from Pete's code and
% performs some stats.

% In: wholeNeuronResults.m
% Out:  trialData containing spike times (Simple and Complex spikes). Every
% struct inside trialData is a neuron. 
% aligned to trial onset and eye stats.

% Associated functions: SpikeTimes2Rate.m, Spiketimes2RateTrial.m,
%                       plot_ProAnti.m

%% Create structure with relevant data (trialData)
%uiopen;
for cellNum = 1:length(wholeNeuronResults);
    for trialNum = 1:length(wholeNeuronResults(cellNum).allStableTrials);
        if ~isempty(wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes)
            
            %behavioral
            trialData(cellNum).trial.behav(trialNum).saccadeOnset =  wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeTime; % saccade onset from trial start (?)
            trialData(cellNum).trial.behav(trialNum).goCueTime = wholeNeuronResults(cellNum).allStableTrials(trialNum).goCueTime; %go cue from trial start
            trialData(cellNum).trial.behav(trialNum).saccAmplitude = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeAmplitude;
            trialData(cellNum).trial.behav(trialNum).saccDuration = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeDuration;
            trialData(cellNum).trial.behav(trialNum).saccPeakVel = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadePeakVel;
            trialData(cellNum).trial.behav(trialNum).reactionTime = wholeNeuronResults(cellNum).allStableTrials(trialNum).reactionTime;
            
            %neural data
            trialData(cellNum).trial.neural(trialNum).tspk_SS =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{1}; % contains spike times for SS aligned to trial onset (?)
            trialData(cellNum).trial.neural(trialNum).tspk_CS =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{2}; % contains spike times for CS aligned to trial onset (?)
            trialData(cellNum).trial.neural(trialNum).tspk_SS_align_sacc =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{1}-trialData(cellNum).trial.behav(trialNum).saccadeOnset; % contains spike times for SS aligned to sacc onset
            trialData(cellNum).trial.neural(trialNum).tspk_CS_align_sacc =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{2}-trialData(cellNum).trial.behav(trialNum).saccadeOnset; % contains spike times for CS aligned to sacc onset
            
            
            % select condition type Pro or Antisaccade
            if ismember(wholeNeuronResults(cellNum).allStableTrials(trialNum).conditionCode, [1 4 5 8 10 11 14 15]);
                trialData(cellNum).trial.behav(trialNum).condition = 'Prosaccade';
            else
                trialData(cellNum).trial.behav(trialNum).condition = 'Antisaccade';
            end
        end
    end
    % Indexes to select Pro and Anti trials
    trialData(cellNum).indx_correctProTrials = wholeNeuronResults(cellNum).selectedTrials.corProTrials;
    trialData(cellNum).indx_correctAntiTrials = wholeNeuronResults(cellNum).selectedTrials.corAntiTrials;
end

clear cellNum trialNum %wholeNeuronResults

%% Spike times to rate for pro and antisaccades and behav
binwidth = 0.01;
timepoints = -4:binwidth:4;
analyse_sacc_win = 0;
for cellNum = 1:length(trialData)
    % Rate for all trials - not distinguishing pro and anti
    [trialData(cellNum).all.neural.nspk,trialData(cellNum).all.neural.ts] = Spiketimes2Rate(trialData(cellNum).trial.neural,timepoints,binwidth,analyse_sacc_win);
    
    %% Pro trials
    correctProTrials = trialData(cellNum).indx_correctProTrials; % index pro trials
    
    % behav
    for trialNum = 1:length(trialData(cellNum).indx_correctProTrials)
        trialData(cellNum).pro.behav.trial(trialNum).saccadeOnset = trialData(cellNum).trial.behav(correctProTrials(trialNum)).saccadeOnset;
        trialData(cellNum).pro.behav.trial(trialNum).goCueTime = trialData(cellNum).trial.behav(correctProTrials(trialNum)).goCueTime;
        trialData(cellNum).pro.behav.trial(trialNum).saccAmplitude = trialData(cellNum).trial.behav(correctProTrials(trialNum)).saccAmplitude;
        trialData(cellNum).pro.behav.trial(trialNum).saccDuration = trialData(cellNum).trial.behav(correctProTrials(trialNum)).saccDuration;
        trialData(cellNum).pro.behav.trial(trialNum).saccPeakVel = trialData(cellNum).trial.behav(correctProTrials(trialNum)).saccPeakVel;
        trialData(cellNum).pro.behav.trial(trialNum).reactionTime = trialData(cellNum).trial.behav(correctProTrials(trialNum)).reactionTime;
    end
    
    % neural
    trialData(cellNum).pro.neural.trial = trialData(cellNum).trial.neural(trialData(cellNum).indx_correctProTrials);
    [trialData(cellNum).pro.neural.nspk,trialData(cellNum).pro.neural.ts] = Spiketimes2Rate(trialData(cellNum).pro.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to trial onset
    % Pro saccade window
    analyse_sacc_win = 1;
    [trialData(cellNum).pro.neural.nspk_sacc,trialData(cellNum).pro.neural.ts_sacc] = Spiketimes2Rate(trialData(cellNum).pro.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to saccade onset
    analyse_sacc_win = 0;
    
    %% Anti trials
    correctAntiTrials = trialData(cellNum).indx_correctAntiTrials; % index pro trials
    % behav
    for trialNum = 1:length(trialData(cellNum).indx_correctAntiTrials)
        trialData(cellNum).anti.behav.trial(trialNum).saccadeOnset = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).saccadeOnset;
        trialData(cellNum).anti.behav.trial(trialNum).goCueTime = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).goCueTime;
        trialData(cellNum).anti.behav.trial(trialNum).saccAmplitude = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).saccAmplitude;
        trialData(cellNum).anti.behav.trial(trialNum).saccDuration = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).saccDuration;
        trialData(cellNum).anti.behav.trial(trialNum).saccPeakVel = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).saccPeakVel;
        trialData(cellNum).anti.behav.trial(trialNum).reactionTime = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).reactionTime;
    end
    %neural
    trialData(cellNum).anti.neural.trial = trialData(cellNum).trial.neural(trialData(cellNum).indx_correctAntiTrials);
    [trialData(cellNum).anti.neural.nspk,trialData(cellNum).anti.neural.ts] = Spiketimes2Rate(trialData(cellNum).anti.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to trial onset
    % Anti saccade window
    analyse_sacc_win = 1;
    [trialData(cellNum).anti.neural.nspk_sacc,trialData(cellNum).anti.neural.ts_sacc] = Spiketimes2Rate(trialData(cellNum).anti.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to saccade onset
    analyse_sacc_win = 0;
    
end
%% For every trial
analyse_sacc_win = 0;%
for cellNum = 1:length(trialData)
    %pro trials
    correctProTrials = trialData(cellNum).indx_correctProTrials; % index pro trials
    for trialNum = 1:length(correctProTrials)
        [trialData(cellNum).pro.neural.trial(trialNum).nspk,trialData(cellNum).pro.neural.trial(trialNum).ts] = Spiketimes2RateTrial(trialData(cellNum).pro.neural.trial(trialNum),timepoints,binwidth,analyse_sacc_win); % all trial
        
        % get spks in sacc window
        sacc_window_pro = trialData(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc;
        win_indx_pro = find(sacc_window_pro>=-0.1 & sacc_window_pro<=0.1); % get spks that happened in sacc window
        trialData(cellNum).pro.neural.trial(trialNum).tspk_saccWin = sacc_window_pro(win_indx_pro);
        analyse_sacc_win = 1;
        [trialData(cellNum).pro.neural.trial(trialNum).nspk_sacc,trialData(cellNum).pro.neural.trial(trialNum).ts_sacc] = Spiketimes2RateTrial(trialData(cellNum).pro.neural.trial(trialNum).tspk_saccWin,timepoints,binwidth,analyse_sacc_win);
        analyse_sacc_win = 0;%
    end
    % anti trials
    correctAntiTrials = trialData(cellNum).indx_correctAntiTrials;
    for trialNum = 1:length(correctAntiTrials)
        [trialData(cellNum).anti.neural.trial(trialNum).nspk,trialData(cellNum).anti.neural.trial(trialNum).ts] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum),timepoints,binwidth,analyse_sacc_win); % all trial
        
        % get spks in sacc window
        sacc_window_anti = trialData(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc;
        win_indx_anti = find(sacc_window_anti>=-0.1 & sacc_window_anti<=0.1); % get spks that happened in sacc window
        trialData(cellNum).anti.neural.trial(trialNum).tspk_saccWin = sacc_window_anti(win_indx_anti);
        analyse_sacc_win = 1;
        [trialData(cellNum).anti.neural.trial(trialNum).nspk_sacc,trialData(cellNum).anti.neural.trial(trialNum).ts_sacc] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum).tspk_saccWin,timepoints,binwidth,analyse_sacc_win);
        analyse_sacc_win = 0;%
    end
end
clear cellNum trialNum correctProTrials correctAntiTrials analyse_sacc_win
% save trialData.m

%%






