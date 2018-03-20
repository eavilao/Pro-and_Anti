% % %recordingListEricPvAElectrophys
% % this is the final list and taes over from
% % recordingListEricSaccadeElectrophys2015 and is the same but now has
% % neurons without a CS pause and a simple spike modualtion removed by
% % commenting and reordered to match Final List PvA.rtf numbering
% % this is now a new version of recordingListEricPvAElectrophysLoc which had
% % the depths and locations of teh cells on. Now this list will be trimmed to
% % only include those neurons that have been resorted (all other commented out) by Eric or Nico
% % manually, this will be represented by a 1 in the 7th column of the cell
% % array .fileList
% % 
% % this is a new list just for Nicos 2016 recordings.
% % 
% % 
% % PH 17/10/16, for the august files and onwards the eye channels seem to
% % have moved from CAI17 and 18 to 1 and 2, adding another field to the list
% % below to indicate this and then need to update the markStablePeriods
% % function to take an optional input argument
% % VERMIS NEURONS
% % %%
neuronList(1).neuronName = '16_03_03_N1';    %
neuronList(1).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160303_0004.mat'          3              [1 2]     []          []      NaN    1
'F160303_0005.mat'          3              [1 2]     []          []      NaN    1
'F160303_0006.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(1).group = 'Pause_Burst'; %not real 
neuronList(1).area = 'vermis';
neuronList(1).depth = []; %not real
neuronList(1).gridLocation = [ ]; %not real
neuronList(1).eyeChannels = [17 18]; %if blank then they are on 17 and 18 (x and y), if on other cahnnels then should be [x y] here
%%%%

neuronList(2).neuronName = '16_03_01_N3';    %
neuronList(2).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160301_0006.mat'          3              [1 2]     []          []      NaN    1
'F160301_0007.mat'          3              [1 2]     []          []      NaN    1
'F160301_0008.mat'          3              [1 2]     []          []      NaN    1
%'F160301_0009.mat'          3              [1 2]     []          []
%NaN    1 %not sorted
};
neuronList(2).group = 'Pause_Burst'; %not real 
neuronList(2).area = 'vermis';
neuronList(2).depth = []; %not real
neuronList(2).gridLocation = [ ]; %not real
neuronList(2).eyeChannels = [17 18];
%%%

neuronList(3).neuronName = '16_02_25_N1';    %
neuronList(3).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160225_0002.mat'          3              [1 2]     []          []      NaN    1
'F160225_0003.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(3).group = 'Pause_Burst'; %not real 
neuronList(3).area = 'vermis';
neuronList(3).depth = []; %not real
neuronList(3).gridLocation = [ ]; %not real
neuronList(3).eyeChannels = [17 18];
% % % 
% % % %%%
% % % 
neuronList(4).neuronName = '16_02_25_N4';    %
neuronList(4).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160225_0007.mat'          3              [1 2]     []          []      NaN    1
'F160225_0008.mat'          3              [1 2]     []          []      NaN    1
'F160225_0009.mat'          3              [1 2]     []          []      NaN    1
'F160225_0010.mat'          3              [1 2]     []          []      NaN    1
% % % NaN    0 %find as note says they are good but I have no spike sorting
% % % files
% 'F160225_0011.mat'          3              [1 2]     []          []      NaN    0
% 'F160225_0012.mat'          3              [1 2]     []          []      NaN    0
};
neuronList(4).group = 'Pause_Burst'; %not real 
neuronList(4).area = 'vermis';
neuronList(4).monkey = 'moshe';
neuronList(4).depth = []; %not real
neuronList(4).gridLocation = [ ]; %not real
neuronList(4).eyeChannels = [17 18];
% % 
neuronList(5).neuronName = '16_02_24_N2';    %
neuronList(5).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160224_0003.mat'          3              [1 2]     []          []      NaN    1
'F160224_0004.mat'          3              [1 2]     []          []      NaN    0
%'F160224_0005.mat'          3              [1 2]     []          []      NaN    0
};
neuronList(5).group = 'Pause_Burst'; %not real 
neuronList(5).area = 'vermis';
neuronList(5).depth =[]; %not real
neuronList(5).gridLocation = [ ]; %not real
neuronList(5).eyeChannels = [17 18];
% % % 
neuronList(6).neuronName = '16_02_24_N5';    %
neuronList(6).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160224_0009.mat'          3              [1 2]     []          []      NaN    1
'F160224_0010.mat'          3              [1 2]     []          []      NaN    1
'F160224_0011.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(6).group = 'Pause_Burst'; %not real 
neuronList(6).area = 'vermis';
neuronList(6).depth = []; %not real
neuronList(6).gridLocation = [ ]; %not real
neuronList(6).eyeChannels = [17 18];
% % 
neuronList(7).neuronName = '16_02_18_N4';    % stimlist file not big
%enough, last file says 16:30 stim list is 15:55, removed last two files
neuronList(7).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160218_0006.mat'          3              [1 2]     []          []      NaN    1
'F160218_0007.mat'          3              [1 2]     []          []      NaN    1
'F160218_0008.mat'          3              [1 2]     []          []      NaN    1
'F160218_0009.mat'          3              [1 2]     []          []      NaN    1
'F160218_0010.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(7).group = 'Pause_Burst'; %not real 
neuronList(7).area = 'vermis';
neuronList(7).depth = []; %not real
neuronList(7).gridLocation = [ ]; %not real
neuronList(7).eyeChannels = [];
% % 
% % 
neuronList(8).neuronName = '16_02_11_N6';    %
neuronList(8).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160211_0009.mat'          3              [1 2]     []          []      NaN    1
'F160211_0010.mat'          3              [1 2]     []          []      NaN    1
'F160211_0011.mat'          3              [1 2]     []          []      NaN    1
'F160211_0012.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(8).group = 'Pause_Burst'; %not real 
neuronList(8).area = 'vermis';
neuronList(8).depth = []; %not real
neuronList(8).gridLocation = [ ]; %not real
neuronList(8).eyeChannels = [];
% % 
% % %september update
% % 
neuronList(9).neuronName = '16_08_17_N1';    %
neuronList(9).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160817_0002.mat'          3              [1 2]     []          []      NaN    1
'F160817_0003.mat'          3              [1 2]     []          []      NaN    1
'F160817_0004.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(9).group = 'Pause_Burst'; %not real 
neuronList(9).area = 'vermis';
neuronList(9).depth = []; %not real
neuronList(9).gridLocation = [ ]; %not real
neuronList(9).eyeChannels = [1 2];
%%%%
neuronList(10).neuronName = '16_08_17_N2';    %
neuronList(10).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160817_0006.mat'          3              [1 2]     []          []      NaN    1
'F160817_0007.mat'          3              [1 2]     []          []      NaN    1
  };
neuronList(10).group = 'Pause_Burst'; %not real 
neuronList(10).area = 'vermis';
neuronList(10).depth = []; %not real
neuronList(10).gridLocation = [ ]; %not real
neuronList(10).eyeChannels = [1 2];
%%%%
neuronList(11).neuronName = '16_08_17_N3';    %
neuronList(11).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160817_0014.mat'          3              [1 2]     []          []      NaN    1
};

neuronList(11).group = 'Pause_Burst'; %not real 
neuronList(11).area = 'vermis';
neuronList(11).depth =  []; %not real
neuronList(11).gridLocation = [ ]; %not real
neuronList(11).eyeChannels = [1 2];
%%%
neuronList(12).neuronName = '16_08_24_N1';    % error in trail extraction
%with too great a stim file num, only 59 trials in the stimList file
neuronList(12).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160824_0002.mat'          3              [1 2]     []          []      NaN    1
'F160824_0003.mat'          3              [1 2]     []          []      NaN    1
'F160824_0004.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(12).group = 'Pause_Burst'; %not real 
neuronList(12).area = 'vermis';
neuronList(12).depth = []; %not real
neuronList(12).gridLocation = [ ]; %not real
neuronList(12).eyeChannels = [1 2];
%
neuronList(13).neuronName = '16_08_23_N1';    % error in calibration as
%looking for a target code that doesn't exist (26)
neuronList(13).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160823_0002.mat'          3              [1 2]     []          []      NaN    1
'F160823_0003.mat'          3              [1 2]     []          []      NaN    1
'F160823_0004.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(13).group = 'Pause_Burst'; %not real 
neuronList(13).area = 'vermis';
neuronList(13).depth = []; %not real
neuronList(13).gridLocation = [ ]; %not real
neuronList(13).eyeChannels = [1 2];
% %%%%
neuronList(14).neuronName = '16_08_23_N2';    %error in calibration as
% looking for a target code that doesn't exist (26)
neuronList(14).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160823_0005.mat'          3              [1 2]     []          []      NaN    1
'F160823_0006.mat'          3              [1 2]     []          []      NaN    1
'F160823_0007.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(14).group = 'Pause_Burst'; %not real 
neuronList(14).area = 'vermis';
neuronList(14).depth = []; %not real
neuronList(14).gridLocation = [ ]; %not real
neuronList(14).eyeChannels = [1 2];
%%%
% error in calibration as
% looking for a target code that doesn't exist (17)
neuronList(15).neuronName = '16_08_23_N3';    % %I renamed this from N2.2 don't know why it was called that
neuronList(15).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160823_0009.mat'          3              [1 2]     []          []      NaN    1
'F160823_0010.mat'          3              [1 2]     []          []      NaN    1
'F160823_0011.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(15).group = 'Pause_Burst'; %not real 
neuronList(15).area = 'vermis';
neuronList(15).depth =  []; %not real
neuronList(15).gridLocation = [ ]; %not real
neuronList(15).eyeChannels = [1 2];
%%%
neuronList(16).neuronName = '16_08_30_N1';    % 
neuronList(16).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160830_0002.mat'          3              [1 2]     []          []      NaN    1
'F160830_0003.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(16).group = 'Pause_Burst'; %not real 
neuronList(16).area = 'vermis';
neuronList(16).depth = []; %not real
neuronList(16).gridLocation = [ ]; %not real
neuronList(16).eyeChannels = [1 2];
%%%%
neuronList(17).neuronName = '16_08_30_N2';    % 
neuronList(17).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160830_0004.mat'          3              [1 2]     []          []      NaN    1
'F160830_0005.mat'          3              [1 2]     []          []      NaN    1
'F160830_0006.mat'          3              [1 2]     []          []      NaN    1
'F160830_0007.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(17).group = 'Pause_Burst'; %not real 
neuronList(17).area = 'vermis';
neuronList(17).depth = []; %not real
neuronList(17).gridLocation = [ ]; %not real
neuronList(17).eyeChannels = [1 2];
%%%%
neuronList(18).neuronName = '16_08_25_N2';    % 
neuronList(18).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160825_0014.mat'          3              [1 2]     []          []      NaN    1
'F160825_0015.mat'          3              [1 2]     []          []      NaN    1
'F160825_0016.mat'          3              [1 2]     []          []      NaN    1
'F160825_0017.mat'          3              [1 2]     []          []      NaN    1
};
neuronList(18).group = 'Pause_Burst'; %not real 
neuronList(18).area = 'vermis';
neuronList(18).depth =[]; %not real
neuronList(18).gridLocation = [ ]; %not real
neuronList(18).eyeChannels = [1 2];
% % %%%%
% % % %%%
neuronList(10).neuronName = 'Mickey_13_06_28_N1';    %
neuronList(10).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'Mickey_13_06_28_0002.mat'          1              [1 2]     []          []      NaN    1
'Mickey_13_06_28_0003.mat'          1              [1 2]     []          []      NaN    1
'Mickey_13_06_28_0004.mat'          1              [1 2]     []          []      NaN    1
};
neuronList(10).group = 'Pause_Burst'; 
neuronList(10).area = 'vermis';
neuronList(10).depth = [];
neuronList(10).gridLocation = [ ];
%%
neuronList(11).neuronName = 'Mickey_14_01_31_N3'; 
neuronList(11).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'Mickey_14_01_31_0003.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0004.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0005.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0006.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0007.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0008.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0009.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0010.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0011.mat'          1              [1 2]     []          []      NaN    0
'Mickey_14_01_31_0012.mat'          1              [1 2]     []          []      NaN    0
  };
neuronList(11).group = 'Burst_Pause'; 
neuronList(11).area = 'vermis';
neuronList(11).depth = 27.9;
neuronList(11).gridLocation = [7 14];
% %%%
neuronList(12).neuronName = 'Mickey_14_02_10_N1';    % DonePC 
neuronList(12).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'Mickey_14_02_10_0002.mat'          1              [1 2]     []          []      NaN    1
'Mickey_14_02_10_0003.mat'          1              [1 2]     []          []      NaN    1
'Mickey_14_02_10_0004.mat'          1              [1 2]     []          []      NaN    1
 };
neuronList(12).group = 'Burst'; 
neuronList(12).area = 'vermis';
neuronList(12).depth = 26.9;
neuronList(12).gridLocation = [7 14];
% % 
neuronList(20).neuronName = 'Moshe_14_11_03_N2';   % Include 
neuronList(20).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag 
'Moshe_14_11_03_0004.mat'          3              [1 2]     []          []      NaN     1
'Moshe_14_11_03_0005.mat'          3              [1 2]     []          []      NaN     1
};
neuronList(20).group = ''; 
neuronList(20).area = 'vermis';
neuronList(20).depth = 13.3;
neuronList(20).gridLocation = [21 5];
% 
neuronList(24).neuronName = 'Moshe_14_11_05_N4';  %Include   
neuronList(24).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'Moshe_14_11_05_0008.mat'          3              [1 2]     []          []      NaN     1
'Moshe_14_11_05_0009.mat'          3              [1 2]     []          []      NaN     1
'Moshe_14_11_05_0010.mat'          3              [1 2]     []          []      NaN     1
};
neuronList(24).group = ''; 
neuronList(24).area = 'vermis';
neuronList(24).depth = 14.3;
neuronList(24).gridLocation = [22 6];
% 
neuronList(26).neuronName = 'Moshe_14_11_08_N2';  % Include  
neuronList(26).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
%'Moshe_14_11_08_0008.mat'          2              [1 2]     []          []      NaN     0
'Moshe_14_11_08_0009.mat'          2              [1 2]     []          []      NaN     1
'Moshe_14_11_08_0010.mat'          2              [1 2]     []          []      NaN     1
'Moshe_14_11_08_0011.mat'          2              [1 2]     []          []      NaN     1
'Moshe_14_11_08_0012.mat'          2              [1 2]     []          []      NaN     1
};
neuronList(26).group = ''; 
neuronList(26).area = 'vermis';
neuronList(26).depth = 14;
neuronList(26).gridLocation = [25 5];

% 
% % 
% %
% %   ___   ___  __ ______   _   _ ______ _    _ _____   ____  _   _  _____ 
% %  |__ \ / _ \/_ |____  | | \ | |  ____| |  | |  __ \ / __ \| \ | |/ ____|
% %     ) | | | || |   / /  |  \| | |__  | |  | | |__) | |  | |  \| | (___  
% %    / /| | | || |  / /   | . ` |  __| | |  | |  _  /| |  | | . ` |\___ \ 
% %   / /_| |_| || | / /    | |\  | |____| |__| | | \ \| |__| | |\  |____) |
% %  |____|\___/ |_|/_/     |_| \_|______|\____/|_|  \_\\____/|_| \_|_____/ 
% %                                                                         
% %  
%  
% %new in 2017

neuronList(29).neuronName = '17_01_03_N1';    % 
neuronList(29).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170103-0003.mat'          6              [1 2]     []          []      NaN    0
'F170103-0004.mat'          6              [1 2]     []          []      NaN    0
'F170103-0005.mat'          6              [1 2]     []          []      NaN    0
'F170103-0006.mat'          6              [1 2]     []          []      NaN    0
'F170103-0007.mat'          6              [1 2]     []          []      NaN    0
'F170103-0008.mat'          6              [1 2]     []          []      NaN    0
'F170103-0009.mat'          6              [1 2]     []          []      NaN    0
'F170103-0010.mat'          6              [1 2]     []          []      NaN    0
'F170103-0011.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(29).group = 'Pause_Burst'; %not real 
neuronList(29).area = 'vermis';
neuronList(29).depth = 32.0; % real
neuronList(29).gridLocation = [23,7]; % real
neuronList(29).eyeChannels = [1 2];
% 


neuronList(30).neuronName = '17_01_04_N1';    
neuronList(30).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170104-0004.mat'          6              [1 2]     []          []      NaN    0
'F170104-0005.mat'          6              [1 2]     []          []      NaN    0
'F170104-0006.mat'          6              [1 2]     []          []      NaN    0
'F170104-0007.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(30).group = 'Pause_Burst'; %not real 
neuronList(30).area = 'vermis';
neuronList(30).depth = 30.7; % real
neuronList(30).gridLocation = [23 7]; % real
neuronList(30).eyeChannels = [1 2];
% 

neuronList(31).neuronName = '17_01_04_N2';    
neuronList(31).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170104-0011.mat'          6              [1 2]     []          []      NaN    0
'F170104-0012.mat'          6              [1 2]     []          []      NaN    0
'F170104-0013.mat'          6              [1 2]     []          []      NaN    0
'F170104-0014.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(31).group = 'Pause_Burst'; %not real 
neuronList(31).area = 'vermis';
neuronList(31).depth = 32.3; % real
neuronList(31).gridLocation = [23 7]; % real
neuronList(31).eyeChannels = [1 2];
% % 
% % 
neuronList(32).neuronName = '17_01_05_N1';    
neuronList(32).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170105-0005.mat'          6              [1 2]     []          []      NaN    0
'F170105-0006.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(32).group = 'Pause_Burst'; %not real 
neuronList(32).area = 'vermis';
neuronList(32).depth = 33.05 ; % real
neuronList(32).gridLocation = [23 7]; % real
neuronList(32).eyeChannels = [1 2];
% 
% % 
% % 
neuronList(33).neuronName = '17_01_05_N2';    
neuronList(33).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170105-0008.mat'          6              [1 2]     []          []      NaN    0
'F170105-0009.mat'          6              [1 2]     []          []      NaN    0
'F170105-0010.mat'          6              [1 2]     []          []      NaN    0
'F170105-0011.mat'          6              [1 2]     []          []      NaN    0
'F170105-0012.mat'          6              [1 2]     []          []      NaN    0
'F170105-0013.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(33).group = 'Pause_Burst'; %not real 
neuronList(33).area = 'vermis';
neuronList(33).depth = 33.25 ; % real
neuronList(33).gridLocation = [23 7]; % real
neuronList(33).eyeChannels = [1 2];

% % 
neuronList(34).neuronName = '17_01_05_N3';    
neuronList(34).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170105-0014.mat'          6              [1 2]     []          []      NaN    0
'F170105-0015.mat'          6              [1 2]     []          []      NaN    0
'F170105-0016.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(34).group = 'Pause_Burst'; %not real 
neuronList(34).area = 'vermis';
neuronList(34).depth = 33.45 ; % real
neuronList(34).gridLocation = [23 7]; % real
neuronList(34).eyeChannels = [1 2];




neuronList(35).neuronName = '17_01_06_N1';    % 
neuronList(35).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170106-0002.mat'          6              [1 2]     []          []      NaN    0
'F170106-0003.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(35).group = 'Pause_Burst'; %not real 
neuronList(35).area = 'vermis';
neuronList(35).depth = 32.6; % real
neuronList(35).gridLocation = [23 7]; % real
neuronList(35).eyeChannels = [1 2];

%  % 
% % 
neuronList(36).neuronName = '17_01_06_N2';    % 
neuronList(36).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170106-0004.mat'          6              [1 2]     []          []      NaN    0
'F170106-0005.mat'          6              [1 2]     []          []      NaN    0
'F170106-0006.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(36).group = 'Pause_Burst'; %not real 
neuronList(36).area = 'vermis';
neuronList(36).depth = 36.18; % real
neuronList(36).gridLocation = [23 7]; % real
neuronList(36).eyeChannels = [1 2];
% %  
% % 
neuronList(37).neuronName = '17_01_09_N1';    % 
neuronList(37).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170109-0003.mat'          6              [1 2]     []          []      NaN    0
'F170109-0004.mat'          6              [1 2]     []          []      NaN    0
'F170109-0005.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(37).group = 'Pause_Burst'; %not real 
neuronList(37).area = 'vermis';
neuronList(37).depth = 29.9; % real
neuronList(37).gridLocation = [23 8]; %real
neuronList(37).eyeChannels = [1 2];
% % 
% 
neuronList(38).neuronName = '17_01_18_N1';    % not enough trials
neuronList(38).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170118-0008.mat'          6              [1 2]     []          []      NaN    0
'F170118-0009.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(38).group = 'Pause_Burst'; %not real 
neuronList(38).area = 'vermis';
neuronList(38).depth = 35.1; %  real
neuronList(38).gridLocation = [23 5]; %  real
neuronList(38).eyeChannels = [1 2];
% 
% % only one file, but with longer file length ~600 mb
neuronList(39).neuronName = '17_01_20_N2';    % 
neuronList(39).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170120-0005.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(39).group = 'Pause_Burst'; %not real 
neuronList(39).area = 'vermis';
neuronList(39).depth = 31.6; % real
neuronList(39).gridLocation = [23 5]; %real
neuronList(39).eyeChannels = [1 2];

% 
neuronList(40).neuronName = '17_01_27_N1';    % one direction with no anti
neuronList(40).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170127-0002.mat'          6              [1 2]     []          []      NaN    0
'F170127-0003.mat'          6              [1 2]     []          []      NaN    0
'F170127-0004.mat'          6              [1 2]     []          []      NaN    0
'F170127-0005.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(40).group = 'Pause_Burst'; %not real 
neuronList(40).area = 'vermis';
neuronList(40).depth = 30.65; %real
neuronList(40).gridLocation = [23 5]; %not real
neuronList(40).eyeChannels = [1 2];
% 
% 
neuronList(41).neuronName = '17_01_27_N2';    
neuronList(41).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170127-0020.mat'          6              [1 2]     []          []      NaN    0
'F170127-0021.mat'          6              [1 2]     []          []      NaN    0
'F170127-0022.mat'          6              [1 2]     []          []      NaN    0
'F170127-0023.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(41).group = 'Pause_Burst'; %not real 
neuronList(41).area = 'vermis';
neuronList(41).depth =  []; %not real
neuronList(41).gridLocation = [ ]; %not real
neuronList(41).eyeChannels = [1 2];
% 
% 
% % % %at least 1 dir with no cor anti
neuronList(42).neuronName = '17_01_27_N3';    %originally labelled N1.2 but changed to prevent errors 
neuronList(42).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170127-0006.mat'          6              [1 2]     []          []      NaN    0
'F170127-0007.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(42).group = 'Pause_Burst'; %not real 
neuronList(42).area = 'vermis';
neuronList(42).depth = []; %not real
neuronList(42).gridLocation = [ ]; %not real
neuronList(42).eyeChannels = [1 2];
% 
% 
neuronList(43).neuronName = '17_01_31_N1';    
neuronList(43).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170131-0007.mat'          6              [1 2]     []          []      NaN    0
'F170131-0009.mat'          6              [1 2]     []          []      NaN    0
'F170131-0010.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(43).group = 'Pause_Burst'; %not real 
neuronList(43).area = 'vermis';
neuronList(43).depth = 36.02; %real
neuronList(43).gridLocation = [23 5]; %real
neuronList(43).eyeChannels = [1 2];
% 
% 
neuronList(44).neuronName = '17_02_01_N1';    
neuronList(44).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170201-0006.mat'          6              [1 2]     []          []      NaN    0
'F170201-0007.mat'          6              [1 2]     []          []      NaN    0
'F170201-0008.mat'          6              [1 2]     []          []      NaN    0
'F170201-0009.mat'          6              [1 2]     []          []      NaN    0
'F170201-0010.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(44).group = 'Pause_Burst'; %not real 
neuronList(44).area = 'vermis';
neuronList(44).depth = 33.90; %real
neuronList(44).gridLocation = [23 5]; %real
neuronList(44).eyeChannels = [1 2];


neuronList(45).neuronName = '17_02_02_N1';    
neuronList(45).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170202-0006.mat'          5              [1 2]     []          []      NaN    0
'F170202-0007.mat'          5              [1 2]     []          []      NaN    0
'F170202-0008.mat'          5              [1 2]     []          []      NaN    0
'F170202-0009.mat'          5              [1 2]     []          []      NaN    0
'F170202-0010.mat'          5              [1 2]     []          []      NaN    0
'F170202-0011.mat'          5              [1 2]     []          []      NaN    0
'F170202-0012.mat'          5              [1 2]     []          []      NaN    0
'F170202-0013.mat'          5              [1 2]     []          []      NaN    0
'F170202-0014.mat'          5              [1 2]     []          []      NaN    0
};
neuronList(45).group = 'Pause_Burst'; %not real 
neuronList(45).area = 'vermis';
neuronList(45).depth = 37.69; %real
neuronList(45).gridLocation = [24 6]; %real
neuronList(45).eyeChannels = [1 2];



neuronList(46).neuronName = '17_02_03_N1'; % MERGED FROM N1 and N1.2   
neuronList(46).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170203-0003.mat'          6              [1 2]     []          []      NaN    0
'F170203-0004.mat'          6              [1 2]     []          []      NaN    0
'F170203-0005.mat'          6              [1 2]     []          []      NaN    0
'F170203-0006.mat'          6              [1 2]     []          []      NaN    0
'F170203-0007.mat'          6              [1 2]     []          []      NaN    0
'F170203-0008.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(46).group = 'Pause_Burst'; %not real 
neuronList(46).area = 'vermis';
neuronList(46).depth = 36.14; %real
neuronList(46).gridLocation = [22 6]; %real
neuronList(46).eyeChannels = [1 2];


neuronList(47).neuronName = '17_02_03_N2';    
neuronList(47).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170203-0009.mat'          6              [1 2]     []          []      NaN    0
'F170203-0010.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(47).group = 'Pause_Burst'; %not real 
neuronList(47).area = 'vermis';
neuronList(47).depth = 39.65; %real
neuronList(47).gridLocation = [22 6]; %real
neuronList(47).eyeChannels = [1 2];


% 
% %
% 
% %   ___   ___  __   __    _   _ ______ _    _ _____   ____  _   _  _____ 
% %  |__ \ / _ \/_ | / /   | \ | |  ____| |  | |  __ \ / __ \| \ | |/ ____|
% %     ) | | | || |/ /_   |  \| | |__  | |  | | |__) | |  | |  \| | (___  
% %    / /| | | || | '_ \  | . ` |  __| | |  | |  _  /| |  | | . ` |\___ \ 
% %   / /_| |_| || | (_) | | |\  | |____| |__| | | \ \| |__| | |\  |____) |
% %  |____|\___/ |_|\___/  |_| \_|______|\____/|_|  \_\\____/|_| \_|_____/ 
% %     
% 

neuronList(48).neuronName = '160203N1';    
neuronList(48).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160203-0008.mat'          1              [1 2]     []          []      NaN    0
'F160203-0009.mat'          1              [1 2]     []          []      NaN    0

};
neuronList(48).group = 'Pause_Burst'; %not real 
neuronList(48).area = 'vermis';
neuronList(48).depth = 45.6; %  from edge of chamber (real) 
neuronList(48).gridLocation = [24 8]; %  real
neuronList(48).eyeChannels = [1 2];


neuronList(49).neuronName = '160203N2';    
neuronList(49).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160203-0011.mat'          1              [1 2]     []          []      NaN    0
'F160203-0012.mat'          1              [1 2]     []          []      NaN    0
'F160203-0013.mat'          1              [1 2]     []          []      NaN    0
};
neuronList(49).group = 'Pause_Burst'; %not real 
neuronList(49).area = 'vermis';
neuronList(49).depth = 49.2; %  from edge of chamber (real) 
neuronList(49).gridLocation = [24 8]; %  real
neuronList(49).eyeChannels = [1 2];
% 
% 
neuronList(50).neuronName = '16_02_03_N3';    
neuronList(50).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160203-0014.mat'          1              [1 2]     []          []      NaN    0
'F160203-0015.mat'          1              [1 2]     []          []      NaN    0
'F160203-0016.mat'          1              [1 2]     []          []      NaN    0
'F160203-0017.mat'          1              [1 2]     []          []      NaN    0
};
neuronList(50).group = 'Pause_Burst'; %not real 
neuronList(50).area = 'vermis';
neuronList(50).depth = 49.4; %  from edge of chamber (real) 
neuronList(50).gridLocation = [24 8]; %  real
neuronList(50).eyeChannels = [1 2];
% 
% 
neuronList(51).neuronName = '160204N1';    
neuronList(51).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160204-0003.mat'          3              [1 2]     []          []      NaN    0
'F160204-0004.mat'          3              [1 2]     []          []      NaN    0
};
neuronList(51).group = 'Pause_Burst'; %not real 
neuronList(51).area = 'vermis';
neuronList(51).depth =  [] ; % not real - driver depth missing  
neuronList(51).gridLocation = [24 8]; %  real
neuronList(51).eyeChannels = [1 2];


neuronList(52).neuronName = '160204N2';    
neuronList(52).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160204-0001.mat'          3              [1 2]     []          []      NaN    0
'F160204-0002.mat'          3              [1 2]     []          []      NaN    0
};
neuronList(52).group = 'Pause_Burst'; %not real 
neuronList(52).area = 'vermis';
neuronList(52).depth = 52.14; %  from edge of chamber (real) 
neuronList(52).gridLocation = [24 8]; %  real
neuronList(52).eyeChannels = [1 2];


neuronList(53).neuronName = '160208N1';    
neuronList(53).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160208-0001.mat'          3              [1 2]     []          []      NaN    0
};
neuronList(53).group = 'Pause_Burst'; %not real 
neuronList(53).area = 'vermis';
neuronList(53).depth = 49.7; %  from edge of chamber (real) 
neuronList(53).gridLocation = [24 8]; %  real
neuronList(53).eyeChannels = [1 2];
% 
% 
neuronList(54).neuronName = '16_02_08_N2';    
neuronList(54).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160208-0003.mat'          3              [1 2]     []          []      NaN    0
'F160208-0004.mat'          3              [1 2]     []          []      NaN    0
'F160208-0005.mat'          3              [1 2]     []          []      NaN    0
'F160208-0006.mat'          3              [1 2]     []          []      NaN    0
'F160208-0007.mat'          3              [1 2]     []          []      NaN    0
};
neuronList(54).group = 'Pause_Burst'; %not real 
neuronList(54).area = 'vermis';
neuronList(54).depth = 54.4; %  from edge of chamber (real) 
neuronList(54).gridLocation = [24 8]; %  real
neuronList(54).eyeChannels = [1 2];
% 
% 
neuronList(55).neuronName = '160210N1';    
neuronList(55).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160210-0002.mat'          3              [1 2]     []          []      NaN    0
'F160210-0003.mat'          3              [1 2]     []          []      NaN    0
};
neuronList(55).group = 'Pause_Burst'; %not real 
neuronList(55).area = 'vermis';
neuronList(55).depth = 48.3; %  from edge of chamber (real) 
neuronList(55).gridLocation = [24 6]; %  real
neuronList(55).eyeChannels = [1 2];


neuronList(56).neuronName = '160210N2';    
neuronList(56).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160210-0004.mat'          3              [1 2]     []          []      NaN    0
 
};
neuronList(56).group = 'Pause_Burst'; %not real 
neuronList(56).area = 'vermis';
neuronList(56).depth =  [] ; % not real - driver depth missing  
neuronList(56).gridLocation = [24 6]; %  real
neuronList(56).eyeChannels = [1 2];


neuronList(57).neuronName = '160210N3';    
neuronList(57).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160210-0005.mat'          3              [1 2]     []          []      NaN    0
};
neuronList(57).group = 'Pause_Burst'; %not real 
neuronList(57).area = 'vermis';
neuronList(57).depth =  [] ; % not real - driver depth missing  
neuronList(57).gridLocation = [24 6]; %  real
neuronList(57).eyeChannels = [1 2];


neuronList(58).neuronName = '160211N2';    
neuronList(58).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160211-0002.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(58).group = 'Pause_Burst'; %not real 
neuronList(58).area = 'vermis';
neuronList(58).depth =  41.48 ; %   real -  
neuronList(58).gridLocation = [24 6]; %  real
neuronList(58).eyeChannels = [1 2];


neuronList(59).neuronName = '160211N3';    
neuronList(59).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160211-0003.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(59).group = 'Pause_Burst'; %not real 
neuronList(59).area = 'vermis';
neuronList(59).depth =  [] ; % not real - driver depth missing  
neuronList(59).gridLocation = [24 6]; %  real
neuronList(59).eyeChannels = [1 2];


neuronList(60).neuronName = '160211N4';    
neuronList(60).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160211-0005.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(60).group = 'Pause_Burst'; %not real 
neuronList(60).area = 'vermis';
neuronList(60).depth =  [] ; % not real - driver depth missing  
neuronList(60).gridLocation = [24 6]; %  real
neuronList(60).eyeChannels = [1 2];


neuronList(61).neuronName = '160211N5';    
neuronList(61).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160211-0006.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(61).group = 'Pause_Burst'; %not real 
neuronList(61).area = 'vermis';
neuronList(61).depth =  45.98 ; % real - 
neuronList(61).gridLocation = [24 6]; %  real
neuronList(61).eyeChannels = [1 2];
% 
% 
neuronList(62).neuronName = '16_02_11_N6';    
neuronList(62).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160211-0009.mat'          3                 [1 2]     []          []      NaN    0
'F160211-0010.mat'          3                 [1 2]     []          []      NaN    0
'F160211-0011.mat'          3                 [1 2]     []          []      NaN    0
'F160211-0012.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(62).group = 'Pause_Burst'; %not real 
neuronList(62).area = 'vermis';
neuronList(62).depth =  45.98 ; % real - 
neuronList(62).gridLocation = [24 6]; %  real
neuronList(62).eyeChannels = [1 2];
% 
% 
neuronList(63).neuronName = '160212N1';    
neuronList(63).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160212-0005.mat'          3                 [1 2]     []          []      NaN    0
'F160212-0006.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(63).group = 'Pause_Burst'; %not real 
neuronList(63).area = 'vermis';
neuronList(63).depth =  48.03 ; % real - 
neuronList(63).gridLocation = [24 6]; %  real
neuronList(63).eyeChannels = [1 2];



neuronList(64).neuronName = '160215N3';    
neuronList(64).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160215-0004.mat'          1                 [1 2]     []          []      NaN    0
};
neuronList(64).group = 'Pause_Burst'; %not real 
neuronList(64).area = 'vermis';
neuronList(64).depth =  43.28 ; % real - 
neuronList(64).gridLocation = [24 8]; %  real
neuronList(64).eyeChannels = [1 2];



neuronList(65).neuronName = '160215N4';    
neuronList(65).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F165215-0005.mat'          1                 [1 2]     []          []      NaN    0
};
neuronList(65).group = 'Pause_Burst'; %not real 
neuronList(65).area = 'vermis';
neuronList(65).depth =  51.28 ; % real - 
neuronList(65).gridLocation = [24 9]; %  real
neuronList(65).eyeChannels = [1 2];



neuronList(66).neuronName = '160218N1';    
neuronList(66).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160218-0002.mat'          3                 [1 2]     []          []      NaN    0
'F160218-0003.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(66).group = 'Pause_Burst'; %not real 
neuronList(66).area = 'vermis';
neuronList(66).depth =  []  ;% not real - 
neuronList(66).gridLocation = [24 7]; %  real
neuronList(66).eyeChannels = [1 2];
% 
% 
% 
neuronList(67).neuronName = '16_02_18_N4';    
neuronList(67).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160218-0006.mat'          3                 [1 2]     []          []      NaN    0
'F160218-0007.mat'          3                 [1 2]     []          []      NaN    0
'F160218-0008.mat'          3                 [1 2]     []          []      NaN    0
'F160218-0009.mat'          3                 [1 2]     []          []      NaN    0
'F160218-0010.mat'          3                 [1 2]     []          []      NaN    0
'F160218-0011.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(67).group = 'Pause_Burst'; %not real 
neuronList(67).area = 'vermis';
neuronList(67).depth = [] ; % not real - 
neuronList(67).gridLocation = [24 7]; %  real
neuronList(67).eyeChannels = [1 2];
% 
% 
% 
neuronList(68).neuronName = '160222N1';    
neuronList(68).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160222-0002.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(68).group = 'Pause_Burst'; %not real 
neuronList(68).area = 'vermis';
neuronList(68).depth =  44.73 ; % real
neuronList(68).gridLocation = [21 7]; %  real
neuronList(68).eyeChannels = [1 2];
% 
% 
neuronList(69).neuronName = '16_02_24_N2';    
neuronList(69).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160224-0003.mat'          3                 [1 2]     []          []      NaN    0
'F160224-0004.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(69).group = 'Pause_Burst'; %not real 
neuronList(69).area = 'vermis';
neuronList(69).depth =  38.08;  % real 
neuronList(69).gridLocation = [23 7]; %  real
neuronList(69).eyeChannels = [1 2];
% 
% 
% 
neuronList(70).neuronName = '160224N4';    
neuronList(70).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160224-0007.mat'          3                 [1 2]     []          []      NaN    0
'F160224-0008.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(70).group = 'Pause_Burst'; %not real 
neuronList(70).area = 'vermis';
neuronList(70).depth =  39.50 ; % real 
neuronList(70).gridLocation = [23 7]; %  real
neuronList(70).eyeChannels = [1 2];
% 
% 
% 
neuronList(71).neuronName = '16_02_24_N5';    
neuronList(71).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160224-0009.mat'          3                 [1 2]     []          []      NaN    0
'F160224-0010.mat'          3                 [1 2]     []          []      NaN    0
'F160224-0011.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(71).group = 'Pause_Burst'; %not real 
neuronList(71).area = 'vermis';
neuronList(71).depth =  [] ;  % not real 
neuronList(71).gridLocation = [23 7]; %  real
neuronList(71).eyeChannels = [1 2];
% 
% 
neuronList(72).neuronName = '16_02_25_N1';    
neuronList(72).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160225-0002.mat'          3                 [1 2]     []          []      NaN    0
'F160225-0003.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(72).group = 'Pause_Burst'; %not real 
neuronList(72).area = 'vermis';
neuronList(72).depth =  30.98 ;  % real 
neuronList(72).gridLocation = [23 8]; %  real
neuronList(72).eyeChannels = [1 2];
% 
% 
neuronList(73).neuronName = '16_02_25_N4';    
neuronList(73).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160225-0007.mat'          3                 [1 2]     []          []      NaN    0
'F160225-0008.mat'          3                 [1 2]     []          []      NaN    0
'F160225-0009.mat'          3                 [1 2]     []          []      NaN    0
'F160225-0010.mat'          3                 [1 2]     []          []      NaN    0
'F160225-0011.mat'          3                 [1 2]     []          []      NaN    0
'F160225-0012.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(73).group = 'Pause_Burst'; %not real 
neuronList(73).area = 'vermis';
neuronList(73).depth =  32.08   ;% real 
neuronList(73).gridLocation = [23 8]; %  real
neuronList(73).eyeChannels = [1 2];
% 
% 
neuronList(74).neuronName = '160229N1';    
neuronList(74).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160229-0002.mat'          3                 [1 2]     []          []      NaN    0
'F160229-0003.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(74).group = 'Pause_Burst'; %not real 
neuronList(74).area = 'vermis';
neuronList(74).depth =  36.68  ; % real 
neuronList(74).gridLocation = [23 8]; %  real
neuronList(74).eyeChannels = [1 2];



neuronList(75).neuronName = '160229N2';    
neuronList(75).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160229-0004.mat'          3                 [1 2]     []          []      NaN    0
'F160229-0005.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(75).group = 'Pause_Burst'; %not real 
neuronList(75).area = 'vermis';
neuronList(75).depth =  36.98 ;  % real 
neuronList(75).gridLocation = [23 8]; %  real
neuronList(75).eyeChannels = [1 2];


neuronList(76).neuronName = '160229N5';    
neuronList(76).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160229-0008.mat'          3                 [1 2]     []          []      NaN    0
'F160229-0009.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(76).group = 'Pause_Burst'; %not real 
neuronList(76).area = 'vermis';
neuronList(76).depth =  39.48 ;  % real 
neuronList(76).gridLocation = [23 8]; %  real
neuronList(76).eyeChannels = [1 2];
% 
% 
% 
neuronList(77).neuronName = '16_03_01_N3';    
neuronList(77).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160301-0006.mat'          3                 [1 2]     []          []      NaN    0
'F160301-0007.mat'          3                 [1 2]     []          []      NaN    0
'F160301-0008.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(77).group = 'Pause_Burst'; %not real 
neuronList(77).area = 'vermis';
neuronList(77).depth =  36.18  ; % real 
neuronList(77).gridLocation = [23 8]; %  real
neuronList(77).eyeChannels = [1 2];
% 
% 
neuronList(78).neuronName = '160301N4';    
neuronList(78).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160301-0011.mat'          3                 [1 2]     []          []      NaN    0
'F160301-0012.mat'          3                 [1 2]     []          []      NaN    0
'F160301-0013.mat'          3                 [1 2]     []          []      NaN    0
'F160301-0014.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(78).group = 'Pause_Burst'; %not real 
neuronList(78).area = 'vermis';
neuronList(78).depth =  45.58 ;  % real 
neuronList(78).gridLocation = [23 8]; %  real
neuronList(78).eyeChannels = [1 2];


neuronList(79).neuronName = '160302N2';    
neuronList(79).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160302-0004.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(79).group = 'Pause_Burst'; %not real 
neuronList(79).area = 'vermis';
neuronList(79).depth =  40.33;   % real 
neuronList(79).gridLocation = [23 8]; %  real
neuronList(79).eyeChannels = [1 2];


neuronList(80).neuronName = '160302N3';    
neuronList(80).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160302-0006.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(80).group = 'Pause_Burst'; %not real 
neuronList(80).area = 'vermis';
neuronList(80).depth =  46.03 ;  % real 
neuronList(80).gridLocation = [23 8]; %  real
neuronList(80).eyeChannels = [1 2];
% 
% 
neuronList(81).neuronName = '16_03_03_N1';    
neuronList(81).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160303-0004.mat'          3                 [1 2]     []          []      NaN    0
'F160303-0005.mat'          3                 [1 2]     []          []      NaN    0
'F160303-0006.mat'          3                 [1 2]     []          []      NaN    0
};
neuronList(81).group = 'Pause_Burst'; %not real 
neuronList(81).area = 'vermis';
neuronList(81).depth =   40.88    ;% real 
neuronList(81).gridLocation = [23 9]; %  real
neuronList(81).eyeChannels = [1 2];
% 
% 
neuronList(82).neuronName = '16_08_17_N1';    
neuronList(82).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160817-0002.mat'          4                 [1 2]     []          []      NaN    0
'F160817-0003.mat'          4                 [1 2]     []          []      NaN    0
'F160817-0004.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(82).group = 'Pause_Burst'; %not real 
neuronList(82).area = 'vermis';
neuronList(82).depth =  32.12  ; % real 
neuronList(82).gridLocation = [23 4]; %  real
neuronList(82).eyeChannels = [1 2];
% 
% 
neuronList(83).neuronName = '16_08_17_N2';    
neuronList(83).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160817-0006.mat'          4                 [1 2]     []          []      NaN    0
'F160817-0007.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(83).group = 'Pause_Burst'; %not real 
neuronList(83).area = 'vermis';
neuronList(83).depth =  34.0 ;  % real 
neuronList(83).gridLocation = [23 4]; %  real
neuronList(83).eyeChannels = [1 2];
% 
% 
% 
neuronList(84).neuronName = '160817N3';    
neuronList(84).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160817-0014.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(84).group = 'Pause_Burst'; %not real 
neuronList(84).area = 'vermis';
neuronList(84).depth =  39.70  ; % real 
neuronList(84).gridLocation = [23 4]; %  real
neuronList(84).eyeChannels = [1 2];
% 
% 
neuronList(85).neuronName = '16_08_23_N1';    
neuronList(85).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160823-0002.mat'          4                 [1 2]     []          []      NaN    0
'F160823-0003.mat'          4                 [1 2]     []          []      NaN    0
'F160823-0004.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(85).group = 'Pause_Burst'; %not real 
neuronList(85).area = 'vermis';
neuronList(85).depth =  []  ; % not real 
neuronList(85).gridLocation = [23 4]; %  real
neuronList(85).eyeChannels = [1 2];
% 
% 
neuronList(86).neuronName = '16_08_23_N2';    
neuronList(86).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160823-0005.mat'          4                 [1 2]     []          []      NaN    0
'F160823-0006.mat'          4                 [1 2]     []          []      NaN    0
'F160823-0007.mat'          4                 [1 2]     []          []      NaN    0
'F160823-0009.mat'          4                 [1 2]     []          []      NaN    0
'F160823-0010.mat'          4                 [1 2]     []          []      NaN    0
'F160823-0011.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(86).group = 'Pause_Burst'; %not real 
neuronList(86).area = 'vermis';
neuronList(86).depth =  []  ; % not real 
neuronList(86).gridLocation = [23 4]; %  real
neuronList(86).eyeChannels = [1 2];
% 
% 
neuronList(87).neuronName = '160824N1';    
neuronList(87).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160874-0002.mat'          4                 [1 2]     []          []      NaN    0
'F160874-0003.mat'          4                 [1 2]     []          []      NaN    0
'F160874-0004.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(87).group = 'Pause_Burst'; %not real 
neuronList(87).area = 'vermis';
neuronList(87).depth =  33.74  ; %  real 
neuronList(87).gridLocation = [23 5]; %  real
neuronList(87).eyeChannels = [1 2];
% 
% 
neuronList(88).neuronName = '16_08_25_N1';    
neuronList(88).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160825-0004.mat'          4                 [1 2]     []          []      NaN    0
'F160825-0005.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(88).group = 'Pause_Burst'; %not real 
neuronList(88).area = 'vermis';
neuronList(88).depth =  41.03  ; %  real 
neuronList(88).gridLocation = [23 5]; %  real
neuronList(88).eyeChannels = [1 2];
% 
% 
neuronList(89).neuronName = '16_08_25_N2';    
neuronList(89).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160825-0014.mat'          4                 [1 2]     []          []      NaN    0
'F160825-0015.mat'          4                 [1 2]     []          []      NaN    0
'F160825-0016.mat'          4                 [1 2]     []          []      NaN    0
'F160825-0017.mat'          4                 [1 2]     []          []      NaN    0
'F160825-0018.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(89).group = 'Pause_Burst'; %not real 
neuronList(89).area = 'vermis';
neuronList(89).depth =  43.03  ; %  real 
neuronList(89).gridLocation = [23 5]; %  real
neuronList(89).eyeChannels = [1 2];
% 
% 
% 
neuronList(90).neuronName = '16_08_30_N1';    
neuronList(90).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160830-0002.mat'          4                 [1 2]     []          []      NaN    0
'F160830-0003.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(90).group = 'Pause_Burst'; %not real 
neuronList(90).area = 'vermis';
neuronList(90).depth =  [];   % not real 
neuronList(90).gridLocation = [22 6]; %  real
neuronList(90).eyeChannels = [1 2];
% 
% 
% 
neuronList(91).neuronName = '16_08_30_N2';    
neuronList(91).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160830-0004.mat'          4                 [1 2]     []          []      NaN    0
'F160830-0005.mat'          4                 [1 2]     []          []      NaN    0
'F160830-0006.mat'          4                 [1 2]     []          []      NaN    0
'F160830-0007.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(91).group = 'Pause_Burst'; %not real 
neuronList(91).area = 'vermis';
neuronList(91).depth =  [] ; % not real 
neuronList(91).gridLocation = [22 6]; %  real
neuronList(91).eyeChannels = [1 2];
% 
% 
neuronList(92).neuronName = '16_08_30_N3';    
neuronList(92).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160830-0008.mat'          4                 [1 2]     []          []      NaN    0
'F160830-0009.mat'          4                 [1 2]     []          []      NaN    0
'F160830-0010.mat'          4                 [1 2]     []          []      NaN    0
'F160830-0011.mat'          4                 [1 2]     []          []      NaN    0
'F160830-0012.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(92).group = 'Pause_Burst'; %not real 
neuronList(92).area = 'vermis';
neuronList(92).depth =  []  ; % not real 
neuronList(92).gridLocation = [22 6]; %  real
neuronList(92).eyeChannels = [1 2];
% 
% 
neuronList(93).neuronName = '16_08_31_N1';    
neuronList(93).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160831-0003.mat'          4                 [1 2]     []          []      NaN    0
'F160831-0004.mat'          4                 [1 2]     []          []      NaN    0
'F160831-0005.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(93).group = 'Pause_Burst'; %not real 
neuronList(93).area = 'vermis';
neuronList(93).depth =  37.80 ;  % real 
neuronList(93).gridLocation = [22 6]; %  real
neuronList(93).eyeChannels = [1 2];
% 
% 
neuronList(94).neuronName = '16_09_01_N1';    
neuronList(94).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F160941-0003.mat'          4                 [1 2]     []          []      NaN    0
'F160941-0004.mat'          4                 [1 2]     []          []      NaN    0
'F160941-0005.mat'          4                 [1 2]     []          []      NaN    0
'F160941-0006.mat'          4                 [1 2]     []          []      NaN    0
};
neuronList(94).group = 'Pause_Burst'; %not real 
neuronList(94).area = 'vermis';
neuronList(94).depth =  41.60  ; % real 
neuronList(94).gridLocation = [21 7]; %  real
neuronList(94).eyeChannels = [1 2];
% 
% 
neuronList(95).neuronName = '16_11_16_N1';    
neuronList(95).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F161116-0004.mat'          6                 [1 2]     []          []      NaN    0
'F161116-0005.mat'          6                 [1 2]     []          []      NaN    0
};
neuronList(95).group = 'Pause_Burst'; %not real 
neuronList(95).area = 'vermis';
neuronList(95).depth =  29.75 ;  % real 
neuronList(95).gridLocation = [23 7]; %  real
neuronList(95).eyeChannels = [1 2];
% 
% 
neuronList(96).neuronName = '16_11_16_N2';    
neuronList(96).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F161116-0010.mat'          6                 [1 2]     []          []      NaN    0
'F161116-0011.mat'          6                 [1 2]     []          []      NaN    0
'F161116-0012.mat'          6                 [1 2]     []          []      NaN    0
};
neuronList(96).group = 'Pause_Burst'; %not real 
neuronList(96).area = 'vermis';
neuronList(96).depth =  30.80  ; % real 
neuronList(96).gridLocation = [23 7]; %  real
neuronList(96).eyeChannels = [1 2];
% 
% 
% 
neuronList(97).neuronName = '16_11_16_N3';    
neuronList(97).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F161116-0018.mat'          6                 [1 2]     []          []      NaN    0
'F161116-0019.mat'          6                 [1 2]     []          []      NaN    0
'F161116-0020.mat'          6                 [1 2]     []          []      NaN    0
};
neuronList(97).group = 'Pause_Burst'; %not real 
neuronList(97).area = 'vermis';
neuronList(97).depth =  34.90 ;  % real 
neuronList(97).gridLocation = [23 7]; %  real
neuronList(97).eyeChannels = [1 2];


neuronList(98).neuronName = '16_11_16_N4';    
neuronList(98).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F161116-0021.mat'          6                 [1 2]     []          []      NaN    0
'F161116-0022.mat'          6                 [1 2]     []          []      NaN    0
'F161116-0023.mat'          6                 [1 2]     []          []      NaN    0
};
neuronList(98).group = 'Pause_Burst'; %not real 
neuronList(98).area = 'vermis';
neuronList(98).depth =  35.60 ;  % real 
neuronList(98).gridLocation = [23 7]; %  real
neuronList(98).eyeChannels = [1 2];

%  
