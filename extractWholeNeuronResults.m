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

% TODO
% - normalize FR
% - discindx
% - pass eye kin plots to separate code

%% Create structure with relevant data (trialData)
%uiopen;
tsmooth = 0.025;

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
    [trialData(cellNum).all.neural.rate_pst,trialData(cellNum).all.neural.ts_pst] = Spiketimes2Rate(trialData(cellNum).trial.neural,timepoints,binwidth,analyse_sacc_win);
    trialData(cellNum).all.neural.rate_pst = smooth_pst(trialData(cellNum).all.neural.rate_pst,binwidth, tsmooth);
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
    [trialData(cellNum).pro.neural.rate_pst,trialData(cellNum).pro.neural.ts_pst] = Spiketimes2Rate(trialData(cellNum).pro.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to trial onset
    trialData(cellNum).pro.neural.rate_pst = smooth_pst(trialData(cellNum).pro.neural.rate_pst,binwidth,tsmooth);
    % Pro saccade window
    analyse_sacc_win = 1;
    [trialData(cellNum).pro.neural.sacc.align_rate_pst,trialData(cellNum).pro.neural.sacc.align_ts_pst] = Spiketimes2Rate(trialData(cellNum).pro.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to saccade onset
    trialData(cellNum).pro.neural.sacc.align_rate_pst = smooth_pst(trialData(cellNum).pro.neural.sacc.align_rate_pst,binwidth,tsmooth);
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
    [trialData(cellNum).anti.neural.rate_pst,trialData(cellNum).anti.neural.ts] = Spiketimes2Rate(trialData(cellNum).anti.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to trial onset
    trialData(cellNum).anti.neural.rate_pst = smooth_pst(trialData(cellNum).anti.neural.rate_pst,binwidth,tsmooth);
    % Anti saccade window
    analyse_sacc_win = 1;
    [trialData(cellNum).anti.neural.sacc.align_rate_pst,trialData(cellNum).anti.neural.sacc.align_ts_pst] = Spiketimes2Rate(trialData(cellNum).anti.neural.trial,timepoints,binwidth,analyse_sacc_win); % aligned to saccade onset
    trialData(cellNum).anti.neural.sacc.align_rate_pst = smooth_pst(trialData(cellNum).anti.neural.sacc.align_rate_pst,binwidth,tsmooth);
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
        win_indx_pro = find(sacc_window_pro>-0.1 & sacc_window_pro<0.1); % get spks that happened in sacc window
        trialData(cellNum).pro.neural.trial(trialNum).sacc.tspk = sacc_window_pro(win_indx_pro);
        analyse_sacc_win = 1;
        [trialData(cellNum).pro.neural.trial(trialNum).sacc.nspk,trialData(cellNum).pro.neural.trial(trialNum).sacc.ts] = Spiketimes2RateTrial(trialData(cellNum).pro.neural.trial(trialNum).sacc.tspk,timepoints,binwidth,analyse_sacc_win); %redefine timepoints from -100 to 100
        analyse_sacc_win = 0;
        
        %get spikes instruction window
        instr_window_pro = trialData(cellNum).pro.neural.trial(trialNum).tspk_SS; 
        win_indx_instr = find(instr_window_pro>0 & instr_window_pro<0.3); % get spks that happened in instr window
        trialData(cellNum).pro.neural.trial(trialNum).instr.tspk = instr_window_pro(win_indx_instr);
        analyse_sacc_win = 1;
        [trialData(cellNum).pro.neural.trial(trialNum).instr.nspk,trialData(cellNum).pro.neural.trial(trialNum).instr.ts] = Spiketimes2RateTrial(trialData(cellNum).pro.neural.trial(trialNum).instr.tspk,timepoints,binwidth,analyse_sacc_win); %redefine timepoints from -100 to 100
        analyse_sacc_win = 0;
        
        
        % get baseline spks anti
        baseline_win = trialData(cellNum).pro.neural.trial(trialNum).tspk_SS; % 200 ms before trial onset
        win_indx_base = find(baseline_win>-0.3 & baseline_win<-0.1);
        trialData(cellNum).pro.neural.trial(trialNum).base.tspk = baseline_win(win_indx_base);
        analyse_sacc_win = 1;
        [trialData(cellNum).pro.neural.trial(trialNum).base.nspk,trialData(cellNum).pro.neural.trial(trialNum).base.ts] = Spiketimes2RateTrial(trialData(cellNum).pro.neural.trial(trialNum).base.tspk,timepoints,binwidth,analyse_sacc_win); %redefine timepoints from -100 to 100
        analyse_sacc_win = 0;
        
        
    end
    % anti trials
    correctAntiTrials = trialData(cellNum).anti.indx_correctAntiTrials;
    for trialNum = 1:length(correctAntiTrials)
        [trialData(cellNum).anti.neural.trial(trialNum).nspk,trialData(cellNum).anti.neural.trial(trialNum).ts] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum),timepoints,binwidth,analyse_sacc_win); % all trial
        
        % get spks in sacc window
        sacc_window_anti = trialData(cellNum).anti.neural.trial(trialNum).tspk_SS_align_sacc;
        win_indx_anti = find(sacc_window_anti>-0.1 & sacc_window_anti<0.1); % get spks that happened in sacc window
        trialData(cellNum).anti.neural.trial(trialNum).sacc.tspk = sacc_window_anti(win_indx_anti);
        analyse_sacc_win = 1;
        [trialData(cellNum).anti.neural.trial(trialNum).sacc.nspk,trialData(cellNum).anti.neural.trial(trialNum).sacc.ts] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum).sacc.tspk,timepoints,binwidth,analyse_sacc_win);
        analyse_sacc_win = 0;%
        
        % get spikes instruction window
        instr_window_anti = trialData(cellNum).anti.neural.trial(trialNum).tspk_SS; 
        win_indx_instr = find(instr_window_anti>0 & instr_window_anti<0.3); % get spks that happened in instr window
        trialData(cellNum).anti.neural.trial(trialNum).instr.tspk = instr_window_anti(win_indx_instr);
        analyse_sacc_win = 1;
        [trialData(cellNum).anti.neural.trial(trialNum).instr.nspk,trialData(cellNum).anti.neural.trial(trialNum).instr.ts] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum).instr.tspk,timepoints,binwidth,analyse_sacc_win); %redefine timepoints from -100 to 100
        analyse_sacc_win = 0;
        
        % get baseline spks anti
        baseline_win = trialData(cellNum).anti.neural.trial(trialNum).tspk_SS; % 200 ms before trial onset
        win_indx_base = find(baseline_win>-0.3 & baseline_win<-0.1);
        trialData(cellNum).anti.neural.trial(trialNum).base.tspk = baseline_win(win_indx_base);
        analyse_sacc_win = 1;
        [trialData(cellNum).anti.neural.trial(trialNum).base.nspk,trialData(cellNum).anti.neural.trial(trialNum).base.ts] = Spiketimes2RateTrial(trialData(cellNum).anti.neural.trial(trialNum).base.tspk,timepoints,binwidth,analyse_sacc_win); %redefine timepoints from -100 to 100
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
instr_pro = []; instr_anti = [];
sacc_pro = []; sacc_anti = [];

for cellNum = 1:length(trialData)
    ntrls_pro = length(trialData(cellNum).pro.neural.trial); 
    ntrls_anti = length(trialData(cellNum).anti.neural.trial); 
    % instruction window pro
    instr_win = trialData(cellNum).pro.neural.ts_pst>0 & trialData(cellNum).pro.neural.ts_pst<0.3; % 200 ms after trial onset
    trialData(cellNum).pro.neural.instr.ts_pst = trialData(cellNum).pro.neural.ts_pst(instr_win);
    trialData(cellNum).pro.neural.instr.rate_pst = trialData(cellNum).pro.neural.rate_pst(instr_win);
    trialData(cellNum).pro.neural.instr.rate_mu = mean(trialData(cellNum).pro.neural.instr.rate_pst);
    trialData(cellNum).pro.neural.instr.rate_sig = std(trialData(cellNum).pro.neural.instr.rate_pst)/sqrt(ntrls_pro); 
    % instruction window anti
    trialData(cellNum).anti.neural.instr.ts_pst = trialData(cellNum).anti.neural.ts(instr_win); 
    trialData(cellNum).anti.neural.instr.rate_pst = trialData(cellNum).anti.neural.rate_pst(instr_win); 
    trialData(cellNum).anti.neural.instr.rate_mu = mean(trialData(cellNum).anti.neural.instr.rate_pst); 
    trialData(cellNum).anti.neural.instr.rate_sig = std(trialData(cellNum).anti.neural.instr.rate_pst)/sqrt(ntrls_anti); 
    
    % saccade window pro
    sacc_win = trialData(cellNum).pro.neural.sacc.align_ts_pst>= -0.101 & trialData(cellNum).pro.neural.sacc.align_ts_pst<=0.101; % 200 ms around sacc
    trialData(cellNum).pro.neural.sacc.ts_pst = trialData(cellNum).pro.neural.sacc.align_ts_pst(sacc_win);
    trialData(cellNum).pro.neural.sacc.rate_pst = trialData(cellNum).pro.neural.sacc.align_rate_pst(sacc_win); 
    trialData(cellNum).pro.neural.sacc.rate_mu = mean(trialData(cellNum).pro.neural.sacc.rate_pst);  
    trialData(cellNum).pro.neural.sacc.rate_sig = std(trialData(cellNum).pro.neural.sacc.rate_pst)/sqrt(ntrls_pro);
    % saccade window anti
    trialData(cellNum).anti.neural.sacc.ts_pst = trialData(cellNum).anti.neural.sacc.align_ts_pst(sacc_win);
    trialData(cellNum).anti.neural.sacc.rate_pst = trialData(cellNum).anti.neural.sacc.align_rate_pst(sacc_win);
    trialData(cellNum).anti.neural.sacc.rate_mu = mean(trialData(cellNum).anti.neural.sacc.rate_pst);
    trialData(cellNum).anti.neural.sacc.rate_sig = std(trialData(cellNum).anti.neural.sacc.rate_pst)/sqrt(ntrls_anti); 
    
    % baseline comparison pro (intertrial interval)
    ntrls_all = length(trialData(cellNum).trial.neural);
    baseline_win = trialData(cellNum).all.neural.ts_pst>=-0.3 & trialData(cellNum).all.neural.ts_pst<=-0.1; % 200 ms before trial onset
    trialData(cellNum).all.neural.base.ts_rate_pst = trialData(cellNum).all.neural.ts_pst(baseline_win);
    trialData(cellNum).all.neural.base.rate_pst = trialData(cellNum).all.neural.rate_pst(baseline_win);
    trialData(cellNum).all.neural.base.rate_mu = mean(trialData(cellNum).all.neural.base.rate_pst);
    trialData(cellNum).all.neural.base.rate_sig = std(trialData(cellNum).all.neural.base.rate_pst)/sqrt(ntrls_all); 
    
    %% Stats per neuron
    %% compare windows against baseline activity for pro and anti (concatenate all trials and make comparison)
    %pro
    for i = 1:length(trialData(cellNum).pro.neural.trial)    
    pro_base(i) = trialData(cellNum).pro.neural.trial(i).base.nspk;
    pro_instr(i)= trialData(cellNum).pro.neural.trial(i).instr.nspk;
    pro_sacc(i)= trialData(cellNum).pro.neural.trial(i).sacc.nspk;
    end
    %anti
    for i = 1:length(trialData(cellNum).anti.neural.trial)
    anti_base(i) = trialData(cellNum).anti.neural.trial(i).base.nspk;
    anti_instr(i)= trialData(cellNum).anti.neural.trial(i).instr.nspk;
    anti_sacc(i)= trialData(cellNum).anti.neural.trial(i).sacc.nspk;
    end
    
     % pro
    [trialData(cellNum).stats.flags.pro.instr_base, trialData(cellNum).stats.pval.pro.instr_base] = ttest((pro_base)',(pro_instr)'); 
    [trialData(cellNum).stats.flags.pro.sacc_base, trialData(cellNum).stats.pval.pro.sacc_base] = ttest((pro_base)',(pro_sacc)'); 
    
    % anti
    [trialData(cellNum).stats.flags.anti.instr_base, trialData(cellNum).stats.pval.anti.instr_base] = ttest((anti_base)',(anti_instr)'); 
    [trialData(cellNum).stats.flags.anti.sacc_base, trialData(cellNum).stats.pval.anti.sacc_base] = ttest((anti_base)',(anti_sacc)'); 
    
    % compare pro vs anti - compute unequal sample test on
    % matrix of trial rate between pro and anti Welchs' test for unequal
    % variances and unequal sample sizes. Per neuron
    
    % instr period pro vs anti
    [trialData(cellNum).stats.flags.proVsAnti_instr,trialData(cellNum).stats.pval.proVsAnti_instr] = ttest2((pro_instr)',(anti_instr)');
    
    % sacc pro vs anti
    [trialData(cellNum).stats.flags.proVsAnti_sacc,trialData(cellNum).stats.pval.proVsAnti_sacc] = ttest2((pro_sacc)',(anti_sacc)');
    
    %% noramlized FR - z-scored
    % aligned to trial onset
trialData(cellNum).pro.neural.sacc.norm_rate_pst = (trialData(cellNum).pro.neural.instr.rate_pst - mean(trialData(cellNum).pro.neural.instr.rate_pst))/...
        std(trialData(cellNum).pro.neural.instr.rate_pst);
    
    % aligned to saccade
    trialData(cellNum).pro.neural.sacc.norm_rate_pst = (trialData(cellNum).pro.neural.sacc.align_rate_pst - mean(trialData(cellNum).pro.neural.sacc.align_rate_pst))/...
        std(trialData(cellNum).pro.neural.sacc.align_rate_pst);
    
    
    % discindx
    
    
    
end




end