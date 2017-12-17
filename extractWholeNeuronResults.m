% extractWholeNeuronResults
% This function extracts data from wholeNeuronResults from Pete's code and
% performs some stats. 

% In: wholeNeuronResults.m
% Out:  trialData containing spike times (Simple and Complex spikes)
% aligned to trial onset and eye stats.

% Associated functions: SpikeTimes2Rate.m, Spiketimes2RateTrial.m,
%                       plot_ProAnti.m


%%%%%%%%%%%%%%    TODO Avg FR window pre, post, saccade, mean and sem. / running avg and stats for every time point 
%% Create structure with relevant data (trialData)
uiopen;
for cellNum = 1:length(wholeNeuronResults);
    for trialNum = 1:length(wholeNeuronResults(cellNum).allStableTrials);        
        if ~isempty(wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes)
            %neural data
            trialData(cellNum).trial(trialNum).tspk_SS =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{1}; % contains spike times for SS aligned to trial onset (?)
            trialData(cellNum).trial(trialNum).tspk_CS =  wholeNeuronResults(cellNum).allStableTrials(trialNum).alignedSpikes{2}; % contains spike times for CS aligned to trial onset (?)
            
            %behavioral
            trialData(cellNum).trial(trialNum).saccadeOnset =  wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeTime; % saccade onset from trial start (?)
            trialData(cellNum).trial(trialNum).goCueTime = wholeNeuronResults(cellNum).allStableTrials(trialNum).goCueTime; %go cue from trial start 
            trialData(cellNum).trial(trialNum).behav.saccAmplitude = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeAmplitude; 
            trialData(cellNum).trial(trialNum).behav.saccDuration = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadeDuration;
            trialData(cellNum).trial(trialNum).behav.saccPeakVel = wholeNeuronResults(cellNum).allStableTrials(trialNum).saccadePeakVel;
            trialData(cellNum).trial(trialNum).reactionTime = wholeNeuronResults(cellNum).allStableTrials(trialNum).reactionTime;
            
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

for i = 1:length(trialData)
    % Rate for all trials - not distinguishing pro and anti
    [trialData(i).all.nspk,trialData(i).all.ts] = Spiketimes2Rate(trialData(i).trial,timepoints,binwidth);
    
    % Pro trials
    trialData(i).pro.trial = trialData(i).trial(trialData(i).indx_correctProTrials);
    [trialData(i).pro.nspk,trialData(i).pro.ts] = Spiketimes2Rate(trialData(i).pro.trial,timepoints,binwidth);
    
    % Anti trials
    trialData(i).anti.trial = trialData(i).trial(trialData(i).indx_correctAntiTrials);
    [trialData(i).anti.nspk,trialData(i).anti.ts] = Spiketimes2Rate(trialData(i).anti.trial,timepoints,binwidth);
end
%% For every trial
for i = 1:length(trialData)
    %pro trials
    for proTrial = 1:length(trialData(i).pro.trial)
        [trialData(i).pro.trial(proTrial).nspk,trialData(i).pro.trial(proTrial).ts] = Spiketimes2RateTrial(trialData(i).pro.trial(proTrial),timepoints,binwidth);
    end 
    % anti trials
    for antiTrial = 1:length(trialData(i).anti.trial)
        [trialData(i).anti.trial(antiTrial).nspk,trialData(i).anti.trial(antiTrial).ts] = Spiketimes2RateTrial(trialData(i).anti.trial(antiTrial),timepoints,binwidth);
    end  
end

clear i j proTrial antiTrial
% save trialData.m

%% Stats 

% Pro vs anti in saccade window (window size = +/- 100 ms around saccade onset)

% Pro vs anti 100 ms after targ onset, which is 100 ms before go cue time

% Running average from target onset to end of saccade, every 10 ms.








    
    

