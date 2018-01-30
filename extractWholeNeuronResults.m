function trialData = extractWholeNeuronResults(wholeNeuronResults)

% extractWholeNeuronResults
% This function extracts data from wholeNeuronResults from Pete's code and
% performs stats.

% In: wholeNeuronResults.m
% Out:  trialData containing spike times (Simple and Complex spikes). Every
% struct inside trialData is a neuron. 
% aligned to trial onset and eye stats.

% Associated functions: SpikeTimes2Rate.m, Spiketimes2RateTrial.m,
%                       plot_ProAnti.m, eyeKinematics_ProAnti.m

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
    trialData(cellNum).pro.indx_correctProTrials = wholeNeuronResults(cellNum).selectedTrials.corProTrials;
    trialData(cellNum).anti.indx_correctAntiTrials = wholeNeuronResults(cellNum).selectedTrials.corAntiTrials;
end

clear cellNum trialNum %wholeNeuronResults

%% Spike times to rate for pro and antisaccades and behav
binwidth = 0.01;
timepoints = -4:binwidth:4;
analyse_sacc_win = 0;
analyse_instr_win=0;

for cellNum = 1:length(trialData)
    % Rate for all trials - not distinguishing pro and anti 
    [trialData(cellNum).all.neural.nspk,trialData(cellNum).all.neural.ts] = Spiketimes2Rate(trialData(cellNum).trial.neural,timepoints,binwidth,analyse_sacc_win);
    
    % Pro trials
    correctProTrials = trialData(cellNum).pro.indx_correctProTrials; % index pro trials
    
    % behav
    for trialNum = 1:length(trialData(cellNum).pro.indx_correctProTrials)
        trialData(cellNum).pro.behav.trial(trialNum).saccadeOnset = trialData(cellNum).trial.behav(correctProTrials(trialNum)).saccadeOnset;
        trialData(cellNum).pro.behav.trial(trialNum).goCueTime = trialData(cellNum).trial.behav(correctProTrials(trialNum)).goCueTime;
        trialData(cellNum).pro.behav.trial(trialNum).saccAmplitude = trialData(cellNum).trial.behav(correctProTrials(trialNum)).saccAmplitude;
        trialData(cellNum).pro.behav.trial(trialNum).saccDuration = trialData(cellNum).trial.behav(correctProTrials(trialNum)).saccDuration;
        trialData(cellNum).pro.behav.trial(trialNum).saccPeakVel = trialData(cellNum).trial.behav(correctProTrials(trialNum)).saccPeakVel;
        trialData(cellNum).pro.behav.trial(trialNum).reactionTime = trialData(cellNum).trial.behav(correctProTrials(trialNum)).reactionTime;
    end
    
    % neural
    trialData(cellNum).pro.neural.trial = trialData(cellNum).trial.neural(trialData(cellNum).pro.indx_correctProTrials);
    [trialData(cellNum).pro.neural.nspk,trialData(cellNum).pro.neural.ts] = Spiketimes2Rate(trialData(cellNum).pro.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to trial onset
    % Pro saccade window
    analyse_sacc_win = 1;
    [trialData(cellNum).pro.neural.nspk_align_sacc,trialData(cellNum).pro.neural.ts_align_sacc] = Spiketimes2Rate(trialData(cellNum).pro.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to saccade onset
    analyse_sacc_win = 0;    
    
    % Anti trials
    correctAntiTrials = trialData(cellNum).anti.indx_correctAntiTrials; % index pro trials
    % behav
    for trialNum = 1:length(trialData(cellNum).anti.indx_correctAntiTrials)
        trialData(cellNum).anti.behav.trial(trialNum).saccadeOnset = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).saccadeOnset;
        trialData(cellNum).anti.behav.trial(trialNum).goCueTime = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).goCueTime;
        trialData(cellNum).anti.behav.trial(trialNum).saccAmplitude = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).saccAmplitude;
        trialData(cellNum).anti.behav.trial(trialNum).saccDuration = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).saccDuration;
        trialData(cellNum).anti.behav.trial(trialNum).saccPeakVel = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).saccPeakVel;
        trialData(cellNum).anti.behav.trial(trialNum).reactionTime = trialData(cellNum).trial.behav(correctAntiTrials(trialNum)).reactionTime;
    end
    %neural
    trialData(cellNum).anti.neural.trial = trialData(cellNum).trial.neural(trialData(cellNum).anti.indx_correctAntiTrials);
    [trialData(cellNum).anti.neural.nspk,trialData(cellNum).anti.neural.ts] = Spiketimes2Rate(trialData(cellNum).anti.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to trial onset
    % Anti saccade window
    analyse_sacc_win = 1;
    [trialData(cellNum).anti.neural.nspk_align_sacc,trialData(cellNum).anti.neural.ts_align_sacc] = Spiketimes2Rate(trialData(cellNum).anti.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to saccade onset
    analyse_sacc_win = 0;
    
end

%% For every trial
analyse_sacc_win = 0;
for cellNum = 1:length(trialData)
    %pro trials
    correctProTrials = trialData(cellNum).pro.indx_correctProTrials; % index pro trials
    for trialNum = 1:length(correctProTrials)
        [trialData(cellNum).pro.neural.trial(trialNum).nspk,trialData(cellNum).pro.neural.trial(trialNum).ts] = Spiketimes2RateTrial(trialData(cellNum).pro.neural.trial(trialNum),timepoints,binwidth,analyse_sacc_win); % all trial
        
        % get spks in sacc window
        sacc_window_pro = trialData(cellNum).pro.neural.trial(trialNum).tspk_SS_align_sacc;
        win_indx_pro = find(sacc_window_pro>=-0.1 & sacc_window_pro<=0.1); % get spks that happened in sacc window
        trialData(cellNum).pro.neural.trial(trialNum).tspk_saccWin = sacc_window_pro(win_indx_pro);
        analyse_sacc_win = 1;
        [trialData(cellNum).pro.neural.trial(trialNum).nspk_sacc,trialData(cellNum).pro.neural.trial(trialNum).ts_sacc] = Spiketimes2RateTrial(trialData(cellNum).pro.neural.trial(trialNum).tspk_saccWin,timepoints,binwidth,analyse_sacc_win); %redefine timepoints from -100 to 100
        analyse_sacc_win = 0;
        
        %get spikes instruction window
        instr_window_pro = trialData(cellNum).pro.neural.trial(trialNum).tspk_SS; 
        win_indx_instr = find(instr_window_pro>=-0.1 & instr_window_pro<=0.1); % get spks that happened in instr window
        trialData(cellNum).pro.neural.trial(trialNum).tspk_instrWin = instr_window_pro(win_indx_instr);
        analyse_sacc_win = 1;
        [trialData(cellNum).pro.neural.trial(trialNum).nspk_instr,trialData(cellNum).pro.neural.trial(trialNum).ts_instr] = Spiketimes2RateTrial(trialData(cellNum).pro.neural.trial(trialNum).tspk_instrWin,timepoints,binwidth,analyse_sacc_win); %redefine timepoints from -100 to 100
        analyse_sacc_win = 0;
        
        
    end
    % anti trials
    correctAntiTrials = trialData(cellNum).anti.indx_correctAntiTrials;
    for trialNum = 1:length(correctAntiTrials)
        [trialData(cellNum).anti.neural.trial(trialNum).nspk,trialData(cellNum).anti.neural.trial(trialNum).ts] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum),timepoints,binwidth,analyse_sacc_win); % all trial
        
        % get spks in sacc window
        sacc_window_anti = trialData(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc;
        win_indx_anti = find(sacc_window_anti>=-0.1 & sacc_window_anti<=0.1); % get spks that happened in sacc window
        trialData(cellNum).anti.neural.trial(trialNum).tspk_saccWin = sacc_window_anti(win_indx_anti);
        analyse_sacc_win = 1;
        [trialData(cellNum).anti.neural.trial(trialNum).nspk_sacc,trialData(cellNum).anti.neural.trial(trialNum).ts_sacc] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum).tspk_saccWin,timepoints,binwidth,analyse_sacc_win);
        analyse_sacc_win = 0;%
        
        % get spikes instruction window
        instr_window_anti = trialData(cellNum).anti.neural.trial(trialNum).tspk_SS; 
        win_indx_instr = find(instr_window_anti>=-0.1 & instr_window_anti<=instr_window_anti+0.1); % get spks that happened in instr window
        trialData(cellNum).anti.neural.trial(trialNum).tspk_instrWin = instr_window_anti(win_indx_instr);
        analyse_sacc_win = 1;
        [trialData(cellNum).anti.neural.trial(trialNum).nspk_instr,trialData(cellNum).anti.neural.trial(trialNum).ts_instr] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum).tspk_instrWin,timepoints,binwidth,analyse_sacc_win); %redefine timepoints from -100 to 100
        analyse_sacc_win = 0;
 
    end
end
clear cellNum trialNum correctProTrials correctAntiTrials analyse_sacc_win

%% Stats
%% Eye kinematics

eyeKin = eyeKinematics_ProAnti(trialData); % extract and plot

proAmp = vertcat(eyeKin(1,:).proAmp); 
antiAmp = vertcat(eyeKin(1,:).antiAmp); 
proDur = vertcat(eyeKin(1,:).proDur); 
antiDur = vertcat(eyeKin(1,:).antiDur); 
proPV = vertcat(eyeKin(1,:).proPV);
antiPV = vertcat(eyeKin(1,:).antiPV); 
proRT = vertcat(eyeKin(1,:).proRT)*1000; 
antiRT = vertcat(eyeKin(1,:).antiRT)*1000; 

% stats per cell
for cellNum = 1:length(trialData)
    %pro
    trialData(cellNum).pro.stats.eye.amp.mu = nanmean(eyeKin(cellNum).proAmp); trialData(cellNum).pro.stats.eye.amp.sig = nanstd(eyeKin(cellNum).proAmp); 
    trialData(cellNum).pro.stats.eye.dur.mu = nanmean(eyeKin(cellNum).proDur); trialData(cellNum).pro.stats.eye.dur.sig = nanstd(eyeKin(cellNum).proDur); 
    trialData(cellNum).pro.stats.eye.pv.mu = nanmean(eyeKin(cellNum).proPV); trialData(cellNum).pro.stats.eye.pv.sig = nanstd(eyeKin(cellNum).proPV); 
    trialData(cellNum).pro.stats.eye.rt.mu = nanmean(eyeKin(cellNum).proRT); trialData(cellNum).pro.stats.eye.rt.sig = nanstd(eyeKin(cellNum).proRT); 
    %anti
    trialData(cellNum).anti.stats.eye.amp.mu = nanmean(eyeKin(cellNum).antiAmp); trialData(cellNum).anti.stats.eye.amp.sig = nanstd(eyeKin(cellNum).antiAmp); 
    trialData(cellNum).anti.stats.eye.dur.mu = nanmean(eyeKin(cellNum).antiDur); trialData(cellNum).anti.stats.eye.dur.sig = nanstd(eyeKin(cellNum).antiDur); 
    trialData(cellNum).anti.stats.eye.pv.mu = nanmean(eyeKin(cellNum).antiPV); trialData(cellNum).anti.stats.eye.pv.sig = nanstd(eyeKin(cellNum).antiPV); 
    trialData(cellNum).anti.stats.eye.rt.mu = nanmean(eyeKin(cellNum).antiRT); trialData(cellNum).anti.stats.eye.rt.sig = nanstd(eyeKin(cellNum).antiRT); 
end


%% Neural
%% Instruction and saccade window
for cellNum = 1:length(trialData)
    ntrls_pro = length(trialData(cellNum).pro.neural.trial); 
    ntrls_anti = length(trialData(cellNum).anti.neural.trial); 
    % instruction window pro
    instr_win = trialData(cellNum).pro.neural.ts>=0 & trialData(cellNum).pro.neural.ts<=0.2; % 200 ms after trial onset
    trialData(cellNum).pro.neural.ts_rate_instr = trialData(cellNum).pro.neural.ts(instr_win);
    trialData(cellNum).pro.neural.rate_instr = trialData(cellNum).pro.neural.nspk(instr_win);
    trialData(cellNum).pro.neural.rate_mu_instr = mean(trialData(cellNum).pro.neural.rate_instr);
    trialData(cellNum).pro.neural.rate_sig_instr = std(trialData(cellNum).pro.neural.rate_instr)/sqrt(ntrls_pro); 
    % instruction window anti
    trialData(cellNum).anti.neural.ts_rate_instr = trialData(cellNum).anti.neural.ts(instr_win); 
    trialData(cellNum).anti.neural.rate_instr = trialData(cellNum).anti.neural.nspk(instr_win); 
    trialData(cellNum).anti.neural.rate_mu_instr = mean(trialData(cellNum).anti.neural.rate_instr); 
    trialData(cellNum).anti.neural.rate_sig_instr = std(trialData(cellNum).anti.neural.rate_instr)/sqrt(ntrls_anti); 
    
    % saccade window pro
    sacc_win = trialData(cellNum).pro.neural.ts_align_sacc>= -0.101 & trialData(cellNum).pro.neural.ts_align_sacc<=0.101; % 200 ms around sacc
    trialData(cellNum).pro.neural.ts_rate_sacc = trialData(cellNum).pro.neural.ts_align_sacc(sacc_win);
    trialData(cellNum).pro.neural.rate_sacc = trialData(cellNum).pro.neural.nspk(sacc_win); 
    trialData(cellNum).pro.neural.rate_mu_sacc = mean(trialData(cellNum).pro.neural.rate_sacc);  
    trialData(cellNum).pro.neural.rate_sig_sacc = std(trialData(cellNum).pro.neural.rate_sacc)/sqrt(ntrls_pro);
    % saccade window anti
    trialData(cellNum).anti.neural.ts_rate_sacc = trialData(cellNum).pro.neural.ts_align_sacc(sacc_win);
    trialData(cellNum).anti.neural.rate_sacc = trialData(cellNum).anti.neural.nspk(sacc_win);
    trialData(cellNum).anti.neural.rate_mu_sacc = mean(trialData(cellNum).anti.neural.rate_sacc);
    trialData(cellNum).anti.neural.rate_sig_sacc = std(trialData(cellNum).anti.neural.rate_sacc)/sqrt(ntrls_anti); 
    
    % baseline comparison pro (intertrial interval) --  check timing
    ntrls_all = length(trialData(cellNum).trial.neural);
    baseline_win = trialData(cellNum).all.neural.ts>=-0.3 & trialData(cellNum).all.neural.ts<=-0.1; % 200 ms before trial onset
    trialData(cellNum).all.neural.ts_rate_base = trialData(cellNum).all.neural.ts(baseline_win);
    trialData(cellNum).all.neural.rate_base = trialData(cellNum).all.neural.nspk(baseline_win);
    trialData(cellNum).all.neural.rate_mu_base = mean(trialData(cellNum).all.neural.rate_base);
    trialData(cellNum).all.neural.rate_sig_base = std(trialData(cellNum).all.neural.rate_base)/sqrt(ntrls_all); 
    
    % stats per neuron
    % compare windows against baseline activity for pro and anti
     % pro
    [trialData(cellNum).stats.flags.pro_instr_base, trialData(cellNum).stats.pval.pro_instr_base] = ttest(trialData(cellNum).all.neural.rate_base,trialData(cellNum).pro.neural.rate_instr); 
    [trialData(cellNum).stats.flags.pro_sacc_base, trialData(cellNum).stats.pval.pro_sacc_base] = ttest(trialData(cellNum).all.neural.rate_base,trialData(cellNum).pro.neural.rate_instr); 
    
    % anti
    [trialData(cellNum).stats.flags.anti_instr_base, trialData(cellNum).stats.pval.anti_instr_base] = ttest(trialData(cellNum).all.neural.rate_base,trialData(cellNum).anti.neural.rate_instr); 
    [trialData(cellNum).stats.flags.anti_sacc_base, trialData(cellNum).stats.pval.anti_sacc_base] = ttest(trialData(cellNum).all.neural.rate_base,trialData(cellNum).anti.neural.rate_instr); 
    
    % compare pro vs anti IN PROGRESS TODO compute unequal sample test on
    % trial spike times
    % instruction
    [trialData(cellNum).stats.pval.proVsAnti_instr,trialData(cellNum).stats.flags.proVsAnti_instr] = ranksum(trialData(cellNum).pro.neural.rate_instr, trialData(cellNum).anti.neural.rate_instr);
    % saccade
    [trialData(cellNum).stats.pval.proVsAnti_sacc,trialData(cellNum).stats.flags.proVsAnti_sacc] = ranksum(trialData(cellNum).pro.neural.rate_sacc, trialData(cellNum).anti.neural.rate_sacc);
    
end

end