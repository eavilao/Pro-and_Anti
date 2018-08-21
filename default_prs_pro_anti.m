% default parameters used to run extractWholeNeuronR

prs.min_trial = 5; % min amount of trials to have to extract a cell

 
% Windows - Make them the same size to make fair comparisons. 
prs.baseline_win = [-0.3 -0.1];  % before saccade onset
prs.instruction_win = [0.1 0.301]; % [0 0.301];  % aligned to trial onset
% prs.saccade_win = [-0.101 0.201];   % aligned to saccade << default
prs.saccade_win = [0 0.201];   %[-0.050 0.150];

% extract spk timing
prs.tspk = [-0.1 0.2]; %100 ms before trial starts to reward +200 ms

% Prosaccade condition codes
prs.proConditions = [1 4 5 8 10 11 14 15];
prs.antiConditions = [ 2 16 7 6 12 3 13 9 ];

 
prs.binwidth = 0.01;
prs.timepoints_instr = -0.1:prs.binwidth:1;
prs.timepoints_sacc = -0.8:prs.binwidth:0.3;

% smoothing window for psth
prs.tsmooth = 0.050;

% probability of spk in pron and anti - bigger window (5 bins = 50 ms)
prs.win_size = 5; % num of bins to take in window. 1 bin = 10 ms

% statistical test to compare if pro == anti - aligned to instr using spk count
prs.signif_criteria = 1.96; %two tailed test

% Sliding window to test diff pro vs anti
% detect saccade related using sliding window (5 bins) 200 ms before
% sacc onset to 200 ms after sacc onset and store time at > or < than 2*std
 prs.slide_win_size = 5; % num of bins to take in sliding window. 1 bin = 10 ms
 prs.slide_win_size_prev = 10; %num of prev bins to take to compare to