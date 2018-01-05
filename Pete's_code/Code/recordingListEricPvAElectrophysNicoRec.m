% %recordingListEricPvAElectrophys
% this is the final list and taes over from
% recordingListEricSaccadeElectrophys2015 and is the same but now has
% neurons without a CS pause and a simple spike modualtion removed by
% commenting and reordered to match Final List PvA.rtf numbering
%this is now a new version of recordingListEricPvAElectrophysLoc which had
%the depths and locations of teh cells on. Now this list will be trimmed to
%only include those neurons that have been resorted (all other commented out) by Eric or Nico
%manually, this will be represented by a 1 in the 7th column of the cell
%array .fileList

%this is a new list just for Nicos 2016 recordings.


%PH 17/10/16, for the august files and onwards the eye channels seem to
%have moved from CAI_17 and 18 to 1 and 2, adding another field to the list
%below to indicate this and then need to update the markStablePeriods
%function to take an optional input argument
%% VERMIS NEURONS
%%%
% neuronList(1).neuronName = '16_03_03_N1';    %
% neuronList(1).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160303_0004.mat'          3              [1 2]     []          []      NaN    1
% 'F160303_0005.mat'          3              [1 2]     []          []      NaN    1
% 'F160303_0006.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(1).group = 'Pause_Burst'; %not real 
% neuronList(1).area = 'vermis';
% neuronList(1).depth = 18; %not real
% neuronList(1).gridLocation = [6 12]; %not real
% neuronList(1).eyeChannels = [17 18]; %if blank then they are on 17 and 18 (x and y), if on other cahnnels then should be [x y] here
% %%%%

% neuronList(2).neuronName = '16_03_01_N3';    %
% neuronList(2).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160301_0006.mat'          3              [1 2]     []          []      NaN    1
% 'F160301_0007.mat'          3              [1 2]     []          []      NaN    1
% 'F160301_0008.mat'          3              [1 2]     []          []      NaN    1
% %'F160301_0009.mat'          3              [1 2]     []          []
% %NaN    1 %not sorted
% };
% neuronList(2).group = 'Pause_Burst'; %not real 
% neuronList(2).area = 'vermis';
% neuronList(2).depth = 18; %not real
% neuronList(2).gridLocation = [6 12]; %not real
% neuronList(2).eyeChannels = [17 18];
% %%%
% 
% neuronList(3).neuronName = '16_02_25_N1';    %
% neuronList(3).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160225_0002.mat'          3              [1 2]     []          []      NaN    1
% 'F160225_0003.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(3).group = 'Pause_Burst'; %not real 
% neuronList(3).area = 'vermis';
% neuronList(3).depth = 18; %not real
% neuronList(3).gridLocation = [6 12]; %not real
% neuronList(3).eyeChannels = [17 18];
% % 
% % %%%
% % 
% neuronList(4).neuronName = '16_02_25_N4';    %
% neuronList(4).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160225_0007.mat'          3              [1 2]     []          []      NaN    1
% 'F160225_0008.mat'          3              [1 2]     []          []      NaN    1
% 'F160225_0009.mat'          3              [1 2]     []          []      NaN    1
% % 'F160225_0010.mat'          3              [1 2]     []          []
% % NaN    0 %find as note says they are good but I have no spike sorting
% % files
% % 'F160225_0011.mat'          3              [1 2]     []          []      NaN    0
% % 'F160225_0012.mat'          3              [1 2]     []          []      NaN    0
% };
% neuronList(4).group = 'Pause_Burst'; %not real 
% neuronList(4).area = 'vermis';
% neuronList(4).depth = 18; %not real
% neuronList(4).gridLocation = [6 12]; %not real
% neuronList(4).eyeChannels = [17 18];
% 
% neuronList(5).neuronName = '16_02_24_N2';    %
% neuronList(5).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160224_0003.mat'          3              [1 2]     []          []      NaN    1
% 'F160224_0004.mat'          3              [1 2]     []          []      NaN    0
% %'F160224_0005.mat'          3              [1 2]     []          []      NaN    0
% };
% neuronList(5).group = 'Pause_Burst'; %not real 
% neuronList(5).area = 'vermis';
% neuronList(5).depth = 18; %not real
% neuronList(5).gridLocation = [6 12]; %not real
% neuronList(5).eyeChannels = [17 18];
% % 
% neuronList(6).neuronName = '16_02_24_N5';    %
% neuronList(6).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160224_0009.mat'          3              [1 2]     []          []      NaN    1
% 'F160224_0010.mat'          3              [1 2]     []          []      NaN    1
% 'F160224_0011.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(6).group = 'Pause_Burst'; %not real 
% neuronList(6).area = 'vermis';
% neuronList(6).depth = 18; %not real
% neuronList(6).gridLocation = [6 12]; %not real
% neuronList(6).eyeChannels = [17 18];
% 
% neuronList(7).neuronName = '16_02_18_N4';    % stimlist file not big
% %enough, last file says 16:30 stim list is 15:55, removed last two files
% neuronList(7).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160218_0006.mat'          3              [1 2]     []          []      NaN    1
% 'F160218_0007.mat'          3              [1 2]     []          []      NaN    1
% 'F160218_0008.mat'          3              [1 2]     []          []      NaN    1
% % 'F160218_0009.mat'          3              [1 2]     []          []      NaN    1
% % 'F160218_0010.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(7).group = 'Pause_Burst'; %not real 
% neuronList(7).area = 'vermis';
% neuronList(7).depth = 18; %not real
% neuronList(7).gridLocation = [6 12]; %not real
% neuronList(7).eyeChannels = [];
% 
% 
% neuronList(8).neuronName = '16_02_11_N6';    %
% neuronList(8).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160211_0009.mat'          3              [1 2]     []          []      NaN    1
% 'F160211_0010.mat'          3              [1 2]     []          []      NaN    1
% 'F160211_0011.mat'          3              [1 2]     []          []      NaN    1
% 'F160211_0012.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(8).group = 'Pause_Burst'; %not real 
% neuronList(8).area = 'vermis';
% neuronList(8).depth = 18; %not real
% neuronList(8).gridLocation = [6 12]; %not real
% neuronList(8).eyeChannels = [];
% 
% %september update
% 
% neuronList(9).neuronName = '16_08_17_N1';    %
% neuronList(9).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160817_0002.mat'          3              [1 2]     []          []      NaN    1
% 'F160817_0003.mat'          3              [1 2]     []          []      NaN    1
% 'F160817_0004.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(9).group = 'Pause_Burst'; %not real 
% neuronList(9).area = 'vermis';
% neuronList(9).depth = 18; %not real
% neuronList(9).gridLocation = [6 12]; %not real
% neuronList(9).eyeChannels = [1 2];
% %%%%
% neuronList(10).neuronName = '16_08_17_N2';    %
% neuronList(10).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160817_0006.mat'          3              [1 2]     []          []      NaN    1
% 'F160817_0007.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(10).group = 'Pause_Burst'; %not real 
% neuronList(10).area = 'vermis';
% neuronList(10).depth = 18; %not real
% neuronList(10).gridLocation = [6 12]; %not real
% neuronList(10).eyeChannels = [1 2];
% %%%%
% neuronList(11).neuronName = '16_08_17_N3';    %
% neuronList(11).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160817_0014.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(11).group = 'Pause_Burst'; %not real 
% neuronList(11).area = 'vermis';
% neuronList(11).depth = 18; %not real
% neuronList(11).gridLocation = [6 12]; %not real
% neuronList(11).eyeChannels = [1 2];
% %%%%
% % neuronList(12).neuronName = '16_08_24_N1';    % error in trail extraction
% % %with too great a stim file num, only 59 trials in the stimList file
% % neuronList(12).fileList = {
% % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % 'F160824_0002.mat'          3              [1 2]     []          []      NaN    1
% % 'F160824_0003.mat'          3              [1 2]     []          []      NaN    1
% % %'F160824_0004.mat'          3              [1 2]     []          []      NaN    1
% % };
% % neuronList(12).group = 'Pause_Burst'; %not real 
% % neuronList(12).area = 'vermis';
% % neuronList(12).depth = 18; %not real
% % neuronList(12).gridLocation = [6 12]; %not real
% % neuronList(12).eyeChannels = [1 2];
% %%%%
% % neuronList(13).neuronName = '16_08_23_N1';    % error in calibration as
% % %looking for a target code that doesn't exist (26)
% % neuronList(13).fileList = {
% % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % 'F160823_0002.mat'          3              [1 2]     []          []      NaN    1
% % 'F160823_0003.mat'          3              [1 2]     []          []      NaN    1
% % 'F160823_0004.mat'          3              [1 2]     []          []      NaN    1
% % };
% % neuronList(13).group = 'Pause_Burst'; %not real 
% % neuronList(13).area = 'vermis';
% % neuronList(13).depth = 18; %not real
% % neuronList(13).gridLocation = [6 12]; %not real
% % neuronList(13).eyeChannels = [1 2];
% %%%%
% % neuronList(14).neuronName = '16_08_23_N2';    %error in calibration as
% % % looking for a target code that doesn't exist (26)
% % neuronList(14).fileList = {
% % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % 'F160823_0005.mat'          3              [1 2]     []          []      NaN    1
% % 'F160823_0006.mat'          3              [1 2]     []          []      NaN    1
% % 'F160823_0007.mat'          3              [1 2]     []          []      NaN    1
% % };
% % neuronList(14).group = 'Pause_Burst'; %not real 
% % neuronList(14).area = 'vermis';
% % neuronList(14).depth = 18; %not real
% % neuronList(14).gridLocation = [6 12]; %not real
% % neuronList(14).eyeChannels = [1 2];
% %%%%
% %error in calibration as
% % % looking for a target code that doesn't exist (17)
% % neuronList(15).neuronName = '16_08_23_N3';    % %I renamed this from N2.2 don't know why it was called that
% % neuronList(15).fileList = {
% % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % 'F160823_0009.mat'          3              [1 2]     []          []      NaN    1
% % 'F160823_0010.mat'          3              [1 2]     []          []      NaN    1
% % 'F160823_0011.mat'          3              [1 2]     []          []      NaN    1
% % };
% % neuronList(15).group = 'Pause_Burst'; %not real 
% % neuronList(15).area = 'vermis';
% % neuronList(15).depth = 18; %not real
% % neuronList(15).gridLocation = [6 12]; %not real
% % neuronList(15).eyeChannels = [1 2];
% %%%%
% neuronList(16).neuronName = '16_08_30_N1';    % 
% neuronList(16).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160830_0002.mat'          3              [1 2]     []          []      NaN    1
% 'F160830_0003.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(16).group = 'Pause_Burst'; %not real 
% neuronList(16).area = 'vermis';
% neuronList(16).depth = 18; %not real
% neuronList(16).gridLocation = [6 12]; %not real
% neuronList(16).eyeChannels = [1 2];
% %%%%
% neuronList(17).neuronName = '16_08_30_N2';    % 
% neuronList(17).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160830_0004.mat'          3              [1 2]     []          []      NaN    1
% 'F160830_0005.mat'          3              [1 2]     []          []      NaN    1
% 'F160830_0006.mat'          3              [1 2]     []          []      NaN    1
% 'F160830_0007.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(17).group = 'Pause_Burst'; %not real 
% neuronList(17).area = 'vermis';
% neuronList(17).depth = 18; %not real
% neuronList(17).gridLocation = [6 12]; %not real
% neuronList(17).eyeChannels = [1 2];
% %%%%
% neuronList(18).neuronName = '16_08_25_N2';    % 
% neuronList(18).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F160825_0014.mat'          3              [1 2]     []          []      NaN    1
% 'F160825_0015.mat'          3              [1 2]     []          []      NaN    1
% 'F160825_0016.mat'          3              [1 2]     []          []      NaN    1
% 'F160825_0017.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(18).group = 'Pause_Burst'; %not real 
% neuronList(18).area = 'vermis';
% neuronList(18).depth = 18; %not real
% neuronList(18).gridLocation = [6 12]; %not real
% neuronList(18).eyeChannels = [1 2];
% %%%%
% % %%%
% % % neuronList(10).neuronName = 'Mickey_13_06_28_N1';    %
% % % neuronList(10).fileList = {
% % % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % % 'Mickey_13_06_28_0002.mat'          1              [1 2]     []          []      NaN    1
% % % 'Mickey_13_06_28_0003.mat'          1              [1 2]     []          []      NaN    1
% % % 'Mickey_13_06_28_0004.mat'          1              [1 2]     []          []      NaN    1
% % % };
% % % neuronList(10).group = 'Pause_Burst'; 
% % % neuronList(10).area = 'vermis';
% % % neuronList(10).depth = 18;
% % % neuronList(10).gridLocation = [6 12];
% % %%%
% % % neuronList(11).neuronName = 'Mickey_14_01_31_N3'; 
% % % neuronList(11).fileList = {
% % % % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % % % 'Mickey_14_01_31_0003.mat'          1              [1 2]     []          []      NaN    0
% % % % 'Mickey_14_01_31_0004.mat'          1              [1 2]     []          []      NaN    0
% % % % 'Mickey_14_01_31_0005.mat'          1              [1 2]     []          []      NaN    0
% % % % 'Mickey_14_01_31_0006.mat'          1              [1 2]     []          []      NaN    0
% % % % 'Mickey_14_01_31_0007.mat'          1              [1 2]     []          []      NaN    0
% % % 'Mickey_14_01_31_0008.mat'          1              [1 2]     []          []      NaN    0
% % % 'Mickey_14_01_31_0009.mat'          1              [1 2]     []          []      NaN    0
% % % 'Mickey_14_01_31_0010.mat'          1              [1 2]     []          []      NaN    0
% % % 'Mickey_14_01_31_0011.mat'          1              [1 2]     []          []      NaN    0
% % % 'Mickey_14_01_31_0012.mat'          1              [1 2]     []          []      NaN    0
% % %   };
% % % neuronList(11).group = 'Burst_Pause'; 
% % % neuronList(11).area = 'vermis';
% % % neuronList(11).depth = 27.9;
% % % neuronList(11).gridLocation = [7 14];
% % % % %%%
% % % neuronList(12).neuronName = 'Mickey_14_02_10_N1';    % DonePC 
% % % neuronList(12).fileList = {
% % % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % % 'Mickey_14_02_10_0002.mat'          1              [1 2]     []          []      NaN    1
% % % 'Mickey_14_02_10_0003.mat'          1              [1 2]     []          []      NaN    1
% % % 'Mickey_14_02_10_0004.mat'          1              [1 2]     []          []      NaN    1
% % %  };
% % % neuronList(12).group = 'Burst'; 
% % % neuronList(12).area = 'vermis';
% % % neuronList(12).depth = 26.9;
% % % neuronList(12).gridLocation = [7 14];
% % 
% % % neuronList(20).neuronName = 'Moshe_14_11_03_N2';   % Include 
% % % neuronList(20).fileList = {
% % % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag 
% % % 'Moshe_14_11_03_0004.mat'          3              [1 2]     []          []      NaN     1
% % % 'Moshe_14_11_03_0005.mat'          3              [1 2]     []          []      NaN     1
% % % };
% % % neuronList(20).group = ''; 
% % % neuronList(20).area = 'vermis';
% % % neuronList(20).depth = 13.3;
% % % neuronList(20).gridLocation = [21 5];
% % % % 
% % % neuronList(24).neuronName = 'Moshe_14_11_05_N4';  %Include   
% % % neuronList(24).fileList = {
% % % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % % 'Moshe_14_11_05_0008.mat'          3              [1 2]     []          []      NaN     1
% % % 'Moshe_14_11_05_0009.mat'          3              [1 2]     []          []      NaN     1
% % % 'Moshe_14_11_05_0010.mat'          3              [1 2]     []          []      NaN     1
% % % };
% % % neuronList(24).group = ''; 
% % % neuronList(24).area = 'vermis';
% % % neuronList(24).depth = 14.3;
% % % neuronList(24).gridLocation = [22 6];
% % % % 
% % % neuronList(26).neuronName = 'Moshe_14_11_08_N2';  % Include  
% % % neuronList(26).fileList = {
% % % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % % %'Moshe_14_11_08_0008.mat'          2              [1 2]     []          []      NaN     0
% % % 'Moshe_14_11_08_0009.mat'          2              [1 2]     []          []      NaN     1
% % % 'Moshe_14_11_08_0010.mat'          2              [1 2]     []          []      NaN     1
% % % 'Moshe_14_11_08_0011.mat'          2              [1 2]     []          []      NaN     1
% % % 'Moshe_14_11_08_0012.mat'          2              [1 2]     []          []      NaN     1
% % % };
% % % neuronList(26).group = ''; 
% % % neuronList(26).area = 'vermis';
% % % neuronList(26).depth = 14;
% % % neuronList(26).gridLocation = [25 5];
% %%
% neuronList(27).neuronName = '16_11_16_N1';    % 
% neuronList(27).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F161116_0004.mat'          6              [1 2]     []          []      NaN    0
% 'F161116_0005.mat'          6              [1 2]     []          []      NaN    0
% % 'F160830_0006.mat'          3              [1 2]     []          []      NaN    1
% % 'F160830_0007.mat'          3              [1 2]     []          []      NaN    1
% };
% neuronList(27).group = 'Pause_Burst'; %not real 
% neuronList(27).area = 'vermis';
% neuronList(27).depth = 18; %not real
% neuronList(27).gridLocation = [6 12]; %not real
% neuronList(27).eyeChannels = [1 2];
% 
% 
% neuronList(28).neuronName = '16_11_22_N1';    % 
% neuronList(28).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% %'F161116_0004.mat'          6              [1 2]     []          []      NaN    0
% 'F161122_0006.mat'          6              [1 2]     []          []      NaN    0
% % 'F160830_0006.mat'          3              [1 2]     []          []      NaN    1
%  'F161122_0007.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(28).group = 'Pause_Burst'; %not real 
% neuronList(28).area = 'vermis';
% neuronList(28).depth = 18; %not real
% neuronList(28).gridLocation = [6 12]; %not real
% neuronList(28).eyeChannels = [1 2];
% 
% %%
% %new in 2017
% neuronList(29).neuronName = '17_01_03_N1';    % 
% neuronList(29).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170103_0003.mat'          6              [1 2]     []          []      NaN    0
% 'F170103_0004.mat'          6              [1 2]     []          []      NaN    0
% 'F170103_0005.mat'          6              [1 2]     []          []      NaN    0
%  'F170103_0006.mat'          6              [1 2]     []          []      NaN    0
%  'F170103_0007.mat'          6              [1 2]     []          []      NaN    0
%  'F170103_0008.mat'          6              [1 2]     []          []      NaN    0
% 'F170103_0009.mat'          6              [1 2]     []          []      NaN    0
%  'F170103_0010.mat'          6              [1 2]     []          []      NaN    0
%  'F170103_0011.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(29).group = 'Pause_Burst'; %not real 
% neuronList(29).area = 'vermis';
% neuronList(29).depth = 18; %not real
% neuronList(29).gridLocation = [6 12]; %not real
% neuronList(29).eyeChannels = [1 2];
% 
% neuronList(30).neuronName = '17_01_06_N1';    % 
% neuronList(30).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170106_0002.mat'          6              [1 2]     []          []      NaN    0
% 'F170106_0003.mat'          6              [1 2]     []          []      NaN    0
% 
% };
% neuronList(30).group = 'Pause_Burst'; %not real 
% neuronList(30).area = 'vermis';
% neuronList(30).depth = 18; %not real
% neuronList(30).gridLocation = [6 12]; %not real
% neuronList(30).eyeChannels = [1 2];
% 
% neuronList(31).neuronName = '17_01_06_N2';    % 
% neuronList(31).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170106_0004.mat'          6              [1 2]     []          []      NaN    0
% 'F170106_0005.mat'          6              [1 2]     []          []      NaN    0
% 'F170106_0006.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(31).group = 'Pause_Burst'; %not real 
% neuronList(31).area = 'vermis';
% neuronList(31).depth = 18; %not real
% neuronList(31).gridLocation = [6 12]; %not real
% neuronList(31).eyeChannels = [1 2];
% 
% %no spi files
% % neuronList(32).neuronName = '17_01_09_N1';    % 
% % neuronList(32).fileList = {
% % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % 'F170109_0003.mat'          6              [1 2]     []          []      NaN    0
% % 'F170109_0004.mat'          6              [1 2]     []          []      NaN    0
% % };
% % neuronList(32).group = 'Pause_Burst'; %not real 
% % neuronList(32).area = 'vermis';
% % neuronList(32).depth = 18; %not real
% % neuronList(32).gridLocation = [6 12]; %not real
% % neuronList(32).eyeChannels = [1 2];
% 
% % neuronList(33).neuronName = '17_01_18_N1';    % not enough trials
% % neuronList(33).fileList = {
% % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % 'F170118_0008.mat'          6              [1 2]     []          []      NaN    0
% % 'F170118_0009.mat'          6              [1 2]     []          []      NaN    0
% % };
% % neuronList(33).group = 'Pause_Burst'; %not real 
% % neuronList(33).area = 'vermis';
% % neuronList(33).depth = 18; %not real
% % neuronList(33).gridLocation = [6 12]; %not real
% % neuronList(33).eyeChannels = [1 2];
% 
% 
% neuronList(34).neuronName = '17_01_27_N1';    % one direction with no anti
% neuronList(34).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170127_0002.mat'          6              [1 2]     []          []      NaN    0
% 'F170127_0003.mat'          6              [1 2]     []          []      NaN    0
% 'F170127_0004.mat'          6              [1 2]     []          []      NaN    0
% 'F170127_0005.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(34).group = 'Pause_Burst'; %not real 
% neuronList(34).area = 'vermis';
% neuronList(34).depth = 18; %not real
% neuronList(34).gridLocation = [6 12]; %not real
% neuronList(34).eyeChannels = [1 2];
% 
% %at least 1 dir with no cor anti
% neuronList(35).neuronName = '17_01_27_N3';    %originally labelled N1.2 but changed to prevent errors 
% neuronList(35).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170127_0006.mat'          6              [1 2]     []          []      NaN    0
% 'F170127_0007.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(35).group = 'Pause_Burst'; %not real 
% neuronList(35).area = 'vermis';
% neuronList(35).depth = 18; %not real
% neuronList(35).gridLocation = [6 12]; %not real
% neuronList(35).eyeChannels = [1 2];
% 
% neuronList(36).neuronName = '17_01_27_N2';    
% neuronList(36).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170127_0020.mat'          6              [1 2]     []          []      NaN    0
% 'F170127_0021.mat'          6              [1 2]     []          []      NaN    0
% 'F170127_0022.mat'          6              [1 2]     []          []      NaN    0
% 'F170127_0023.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(36).group = 'Pause_Burst'; %not real 
% neuronList(36).area = 'vermis';
% neuronList(36).depth = 18; %not real
% neuronList(36).gridLocation = [6 12]; %not real
% neuronList(36).eyeChannels = [1 2];
% 
% neuronList(37).neuronName = '17_01_05_N1';    
% neuronList(37).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170105_0005.mat'          6              [1 2]     []          []      NaN    0
% 'F170105_0006.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(37).group = 'Pause_Burst'; %not real 
% neuronList(37).area = 'vermis';
% neuronList(37).depth = 18; %not real
% neuronList(37).gridLocation = [6 12]; %not real
% neuronList(37).eyeChannels = [1 2];
% 
neuronList(38).neuronName = '17_01_05_N2';    
neuronList(38).fileList = {
% File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
'F170105_0008.mat'          6              [1 2]     []          []      NaN    0
'F170105_0009.mat'          6              [1 2]     []          []      NaN    0
'F170105_0010.mat'          6              [1 2]     []          []      NaN    0
'F170105_0011.mat'          6              [1 2]     []          []      NaN    0
'F170105_0012.mat'          6              [1 2]     []          []      NaN    0
'F170105_0013.mat'          6              [1 2]     []          []      NaN    0
};
neuronList(38).group = 'Pause_Burst'; %not real 
neuronList(38).area = 'vermis';
neuronList(38).depth = 18; %not real
neuronList(38).gridLocation = [6 12]; %not real
neuronList(38).eyeChannels = [1 2];

% neuronList(39).neuronName = '17_01_05_N3';    
% neuronList(39).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170105_0014.mat'          6              [1 2]     []          []      NaN    0
% 'F170105_0015.mat'          6              [1 2]     []          []      NaN    0
% 'F170105_0016.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(39).group = 'Pause_Burst'; %not real 
% neuronList(39).area = 'vermis';
% neuronList(39).depth = 18; %not real
% neuronList(39).gridLocation = [6 12]; %not real
% neuronList(39).eyeChannels = [1 2];
% 
% 
% neuronList(40).neuronName = '17_01_04_N1';    
% neuronList(40).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170104_0004.mat'          6              [1 2]     []          []      NaN    0
% 'F170104_0005.mat'          6              [1 2]     []          []      NaN    0
% 'F170104_0006.mat'          6              [1 2]     []          []      NaN    0
% 'F170104_0007.mat'          6              [1 2]     []          []      NaN    0
% 'F170104_0008.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(40).group = 'Pause_Burst'; %not real 
% neuronList(40).area = 'vermis';
% neuronList(40).depth = 18; %not real
% neuronList(40).gridLocation = [6 12]; %not real
% neuronList(40).eyeChannels = [1 2];
% 
% 
% % neuronList(41).neuronName = '17_01_04_N2';    
% % neuronList(41).fileList = {
% % % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% % 'F170104_0004.mat'          6              [1 2]     []          []      NaN    0
% % 'F170104_0005.mat'          6              [1 2]     []          []      NaN    0
% % 'F170104_0006.mat'          6              [1 2]     []          []      NaN    0
% % 'F170104_0007.mat'          6              [1 2]     []          []      NaN    0
% % 'F170104_0008.mat'          6              [1 2]     []          []      NaN    0
% % };
% % neuronList(41).group = 'Pause_Burst'; %not real 
% % neuronList(41).area = 'vermis';
% % neuronList(41).depth = 18; %not real
% % neuronList(41).gridLocation = [6 12]; %not real
% % neuronList(41).eyeChannels = [1 2];
% 
% 
% neuronList(42).neuronName = '17_02_02_N1';    
% neuronList(42).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170202_0006.mat'          5              [1 2]     []          []      NaN    0
% 'F170202_0007.mat'          5              [1 2]     []          []      NaN    0
% 'F170202_0008.mat'          5              [1 2]     []          []      NaN    0
% 'F170202_0009.mat'          5              [1 2]     []          []      NaN    0
% 'F170202_0010.mat'          5              [1 2]     []          []      NaN    0
% 'F170202_0011.mat'          5              [1 2]     []          []      NaN    0
% 'F170202_0012.mat'          5              [1 2]     []          []      NaN    0
% 'F170202_0013.mat'          5              [1 2]     []          []      NaN    0
% 'F170202_0014.mat'          5              [1 2]     []          []      NaN    0
% };
% neuronList(42).group = 'Pause_Burst'; %not real 
% neuronList(42).area = 'vermis';
% neuronList(42).depth = 18; %not real
% neuronList(42).gridLocation = [6 12]; %not real
% neuronList(42).eyeChannels = [1 2];