% extractWholeNeuronResults
% This function extracts data from wholeNeuronResults from Pete's code and
% performs some stats. 

% In: wholeNeuronResults.m
% Out:  trialData containing spike times (Simple and Complex spikes)
% aligned to trial onset and eye stats.

% Associated functions: SpikeTimes2Rate.m, Spiketimes2RateTrial.m,
%                       plot_ProAnti.m


%%%%%%%%%%%%%%    TODO Avg FR window - saccade, mean and sem. / running avg and stats for every time point 
%% Create structure with relevant data (trialData)
uiopen;
for cellNum = 1:length(wholeNeuronResults);
    for trialNum = 1:length(wholeNeuronResults(cellNum).allStableTrials);        
        if ~isempty(wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes)
            
            %behavioral
            trialData(cellNum).trial(trialNum).saccadeOnset =  wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeTime; % saccade onset from trial start (?)
            trialData(cellNum).trial(trialNum).goCueTime = wholeNeuronResults(cellNum).allStableTrials(trialNum).goCueTime; %go cue from trial start 
            trialData(cellNum).trial(trialNum).behav.saccAmplitude = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeAmplitude; 
            trialData(cellNum).trial(trialNum).behav.saccDuration = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeDuration;
            trialData(cellNum).trial(trialNum).behav.saccPeakVel = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadePeakVel;
            trialData(cellNum).trial(trialNum).reactionTime = wholeNeuronResults(cellNum).allStableTrials(trialNum).reactionTime;
            
            %neural data
            trialData(cellNum).trial(trialNum).tspk_SS =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{1}; % contains spike times for SS aligned to trial onset (?)
            trialData(cellNum).trial(trialNum).tspk_CS =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{2}; % contains spike times for CS aligned to trial onset (?)
            trialData(cellNum).trial(trialNum).tspk_SS_align_sacc =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{1}-trialData(cellNum).trial(trialNum).saccadeOnset; % contains spike times for SS aligned to sacc onset
            trialData(cellNum).trial(trialNum).tspk_CS_align_sacc =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{2}-trialData(cellNum).trial(trialNum).saccadeOnset; % contains spike times for CS aligned to sacc onset
            
            
            % select condition type Pro or Antisaccade
            if ismember(wholeNeuronResults(cellNum).allStableTrials(trialNum).conditionCode, [1 4 5 8 10 11 14 15]);
                trialData(cellNum).trial(trialNum).condition = 'Prosaccade';
            else
                trialData(cellNum).trial(trialNum).condition = 'Antisaccade';
            end
        end
    end
    % Indexes to select Pro and Anti trials
    trialData(cellNum).indx_correctProTrials = wholeNeuronResults(cellNum).selectedTrials.corProTrials;
    trialData(cellNum).indx_correctAntiTrials = wholeNeuronResults(cellNum).selectedTrials.corAntiTrials;
end

clear cellNum trialNum %wholeNeuronResults

%% Spike times to rate for pro and antisaccades
binwidth = 0.01;
timepoints = -4:binwidth:4;
analyse_sacc_win = 0; 
for i = 1:length(trialData)
    % Rate for all trials - not distinguishing pro and anti
    [trialData(i).all.nspk,trialData(i).all.ts] = Spiketimes2Rate(trialData(i).trial,timepoints,binwidth,analyse_sacc_win);
    
    % Pro trials
    trialData(i).pro.trial = trialData(i).trial(trialData(i).indx_correctProTrials);
    [trialData(i).pro.nspk,trialData(i).pro.ts] = Spiketimes2Rate(trialData(i).pro.trial,timepoints,binwidth,analyse_sacc_win); % aligned to trial onset
    % Pro saccade window
    analyse_sacc_win = 1;
    [trialData(i).pro.nspk_sacc,trialData(i).pro.ts_sacc] = Spiketimes2Rate(trialData(i).pro.trial,timepoints,binwidth,analyse_sacc_win); % aligned to saccade onset
    analyse_sacc_win = 0;
    
    % Anti trials
    trialData(i).anti.trial = trialData(i).trial(trialData(i).indx_correctAntiTrials);
    [trialData(i).anti.nspk,trialData(i).anti.ts] = Spiketimes2Rate(trialData(i).anti.trial,timepoints,binwidth,analyse_sacc_win); % aligned to trial onset
    % Anti saccade window
    analyse_sacc_win = 1;
    [trialData(i).anti.nspk_sacc,trialData(i).anti.ts_sacc] = Spiketimes2Rate(trialData(i).anti.trial,timepoints,binwidth,analyse_sacc_win); % aligned to saccade onset
    analyse_sacc_win = 0;
    
end
%% For every trial
analyse_sacc_win = 0;%
for cellNum = 1:length(trialData)
    %pro trials
    for proTrial = 1:length(trialData(cellNum).pro.trial)
        [trialData(cellNum).pro.trial(proTrial).nspk,trialData(cellNum).pro.trial(proTrial).ts] = Spiketimes2RateTrial(trialData(cellNum).pro.trial(proTrial),timepoints,binwidth,analyse_sacc_win); % all trial
        
        % get spks in sacc window
        sacc_aligned = trialData(cellNum).pro.trial(proTrial).tspk_SS -trialData(cellNum).pro.trial(proTrial).saccadeOnset;
        win_indx_pro = find(sacc_window>=-0.1 & sacc_window<=0.1); % get spks that happened in sacc window
        trialData(cellNum).pro.trial(proTrial).tspk_saccWin = sacc_window(win_indx_pro);
        analyse_sacc_win = 1;
        [trialData(cellNum).pro.trial(proTrial).nspk_sacc,trialData(cellNum).pro.trial(proTrial).ts_sacc] = Spiketimes2RateTrial(trialData(cellNum).pro.trial(proTrial).tspk_saccWin,timepoints,binwidth,analyse_sacc_win);
        analyse_sacc_win = 0;%
    end
    % anti trials
    for antiTrial = 1:length(trialData(cellNum).anti.trial)
        [trialData(cellNum).anti.trial(antiTrial).nspk,trialData(cellNum).anti.trial(antiTrial).ts] = Spiketimes2RateTrial(trialData(cellNum).anti.trial(antiTrial),timepoints,binwidth,analyse_sacc_win);
        
        % get spks in sacc window
        sacc_aligned = trialData(cellNum).anti.trial(antiTrial).tspk_SS -trialData(cellNum).anti.trial(antiTrial).saccadeOnset;
        win_indx = find(sacc_window>=-0.1 & sacc_window<=0.1); % get spks that happened in sacc window
        trialData(cellNum).anti.trial(antiTrial).tspk_saccWin = sacc_window(win_indx);
        analyse_sacc_win = 1;
        [trialData(cellNum).anti.trial(antiTrial).nspk_sacc,trialData(cellNum).anti.trial(antiTrial).ts_sacc] = Spiketimes2RateTrial(trialData(cellNum).anti.trial(antiTrial).tspk_saccWin,timepoints,binwidth,analyse_sacc_win);
        analyse_sacc_win = 0;%
    end
    
    
end

clear i j proTrial antiTrial analyse_sacc_win
% save trialData.m

%% Stats 
  
  



    
    

