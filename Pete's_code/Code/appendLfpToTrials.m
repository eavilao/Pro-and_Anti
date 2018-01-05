function trialData = appendLfpToTrials(trialData,fileName,channelNum,windowSize)
%function for Eric's data that should load the lfp data ofa  file and take
%a snippit of given length aligned to trial onset then add it to the trial
%strcutre.



%newTrialdata = trialData;
numTrials = size(trialData,2);




%build variable strings
channelType = 'LFP';
dataString = ['C' channelType '_' sprintf('%03d',channelNum)];
dataFsString = [dataString '_KHz'];
dataStartTimeString = [dataString '_TimeBegin'];

%build variable list to load
varList = {dataString,dataFsString,dataStartTimeString,...
   };

%load data
load(fileName,varList{:})

%TODO replace with dynamic field names?
eval(['lfpData = double(' dataString ');']);
eval(['lfpFs = ' dataFsString '*1000;']);
eval(['lfpDataStartTime = ' dataStartTimeString ';']);

%get window size in samples
windowSizeSamples = [floor(windowSize(1)*lfpFs) ceil(windowSize(2)*lfpFs)];
lastDataSample = length(lfpData);


for trialNum = 1:numTrials
    
    
    disp(['processing trial number ' num2str(trialNum) ' of ' num2str(numTrials)])
    %get start and end time of this trial
    thisTrialStart = trialData(trialNum).trialStart;
    thisTrialEnd = trialData(trialNum).trialEnd;
    goCueTime = trialData(trialNum).bit3time;
    %if it is a 'valid' trial (ie both a start, end and go cue) then check if tehre
    %is a saccade(s) within it
    if ~isempty(thisTrialStart) && ~isempty(thisTrialEnd) && ~isempty(goCueTime)
        trialData(trialNum).lfpEpoch = nan(sum(windowSizeSamples)+1,1);
        %now convert trialStart to samples and get lfp sample
        trialStartSample = round(thisTrialStart*lfpFs);
        
        if trialStartSample-windowSizeSamples(1)>1 && trialStartSample+windowSizeSamples(2)<lastDataSample
            trialData(trialNum).lfpEpoch = lfpData(trialStartSample-windowSizeSamples(1):trialStartSample+windowSizeSamples(2))';
            
            
        end
    else
        trialData(trialNum).lfpEpoch = nan(sum(windowSizeSamples)+1,1);
    end
    
        %TODO should make spectrogram plots as welll.
end
end
