%load and edit Pete
function loadAndEditPete(fileName,spkChannelNum)
%file and channel to load and edit
%should detect if filename is actualyl cell and if so load all and concatenate 
%build variable strings
channelType = 'Spk';
dataString = ['C' channelType '_' sprintf('%03d',spkChannelNum)];
dataFsString = [dataString '_KHz'];
dataStartTimeString = [dataString '_TimeBegin'];

%build variable list to load
varList = {dataString,dataFsString,dataStartTimeString};


if iscell(fileName)
    numFiles = numel(fileName);
    stInit = cell(numFiles,1);
    %fltSignal is same as rawSignal in this case so only do after
    %fs is same for all so do after as are lpFiltFreq, hpFiltFreq and timesToExclude
    allFiles = struct('rawSignal',stInit,'spikeTimes',stInit,'startTime',stInit,'stopTime',stInit,'fileName',stInit);
    for fileNum = 1:numFiles
        
        
        load(fileName{fileNum},varList{:})
        eval(['spkData = double(' dataString ');']);
        eval(['spkFs = ' dataFsString '*1000;']);
        eval(['spkDataStartTime = ' dataStartTimeString ';']);
        load([fileName{fileNum}(1:end-4) '_0' num2str(spkChannelNum) '.spi'],'-mat');
        
        allFiles(fileNum).rawSignal = spkData;
        activeUnits = unique(S.NetId);
        activeUnits = activeUnits(activeUnits~=0);
        spikeTimes = cell(1,max(activeUnits));
        for unitNum = activeUnits'
            theseSpikes = S.Tm(S.NetId==unitNum);
            spikeTimes{unitNum} = theseSpikes;
        end
        clear spkData
        allFiles(fileNum).spikeTimes = spikeTimes;
       allFiles(fileNum).startTime = spkDataStartTime;
       allFiles(fileNum).stopTime = allFiles(fileNum).startTime+length(allFiles(fileNum).rawSignal)/spkFs;
       temp = strsplit(fileName{fileNum},filesep);
        allFiles(fileNum).fileName = temp{end};
    end
    
    
    %now we need to concatenate all
    %first find first file in time
    [startFileTime startFileNum] = min([allFiles.startTime]);
    
    %now find the total number of samples we need
    [lastFileTime lastFileNum] = max([allFiles.startTime]);
    
    %now need to add length of last file
    totalDataSamples = length(allFiles(lastFileNum).rawSignal)+ceil((lastFileTime-startFileTime)*spkFs);
    fullTimeVector = linspace(0,totalDataSamples/spkFs,totalDataSamples);
    fullDataVector = nan(size(fullTimeVector));
    numUnitsPerFile = cellfun(@(x) numel(x),{allFiles.spikeTimes},'uniformoutput',false);
    bigNumUnits = max([numUnitsPerFile{:}])
    allSpikeTimes = cell(1, bigNumUnits);
    [i fileOrder] = sort([allFiles.startTime]);
    for sortFileNum = fileOrder
        totalTimeFromStart = allFiles(sortFileNum).startTime-startFileTime;
        startSample = ceil(totalTimeFromStart*spkFs)+1;
        numSamples = ceil(length(allFiles(sortFileNum).rawSignal));
        fullDataVector(startSample:startSample+numSamples-1) = allFiles(sortFileNum).rawSignal;
       
        allSpikeTimes = vertcat(allSpikeTimes,cellfun(@(x) x+totalTimeFromStart,allFiles(sortFileNum).spikeTimes,'uniformoutput',false))
         
    end
    %tODO make tloopto include noie units as well
    spikeTimes{1} = vertcat(allSpikeTimes{:,1});
    spikeTimes{2}  =vertcat(allSpikeTimes{:,2});
    
    
     channelData.rawSignal =  fullDataVector;
    channelData.filtSignal =  fullDataVector;
    channelData.fullTimeVector = fullTimeVector;
    channelData.fs = spkFs;
    channelData.lpFiltFreq = NaN;
    channelData.hpFiltFreq = NaN;
    channelData.timesToExclude = [];
    %this is anew optional argument to allow the sorter to save the results
    %back to the original format
    %needs to be a cell (numFiles,2), first column is name (only name not
    %path), second is the start and stop time for this file.
    multFileInfo = cell(numFiles,2)
    for fileNum = 1:numFiles
    multFileInfo{fileNum,1} = allFiles(fileOrder(fileNum)).fileName;
     multFileInfo{fileNum,2}(1) = allFiles(fileOrder(fileNum)).startTime-startFileTime;
     multFileInfo{fileNum,2}(2) = allFiles(fileOrder(fileNum)).stopTime-startFileTime;
   
    end
    [hSpikeGui] = launchSpikeGui(channelData,spikeTimes,'multipleFileInfo',multFileInfo,'colourScheme',2)
else
    %
    
    load(fileName,varList{:})
    eval(['spkData = double(' dataString ');']);
    eval(['spkFs = ' dataFsString '*1000;']);
    eval(['spkDataStartTime = ' dataStartTimeString ';']);
    load([fileName(1:end-4) '_0' num2str(spkChannelNum) '.spi'],'-mat');
    
    %generate the following structure
    
    %channelData.rawSignal - vecotr of raw signal
    %                                       .fs - sampling frequency (Hz)
    %                                       .filtSignal - filteredSignal
    %                                       .lpFiltFreq - low pass freq in Hz
    %                                       .hpFiltFreq - high pass freq in Hz
    %                                       .timesToInclude - Currently unused so set to [-inf inf]
    %                                       .timesToExclude - table with times to exclude form further analysis (periodNum,[start end])
    
    channelData.rawSignal = spkData;
    channelData.filtSignal = spkData;
    channelData.fullTimeVector = linspace(0,length(spkData)/spkFs,length(spkData));
    channelData.fs = spkFs;
    channelData.lpFiltFreq = NaN;
    channelData.hpFiltFreq = NaN;
    channelData.timesToExclude = [];
    %loop over units and build spikeTimes cell
    activeUnits = unique(S.NetId);
    activeUnits = activeUnits(activeUnits~=0);
    spikeTimes = cell(1,max(activeUnits));
    for unitNum = activeUnits'
        theseSpikes = S.Tm(S.NetId==unitNum);
        spikeTimes{unitNum} = theseSpikes;
    end
    [hSpikeGui] = launchSpikeGui(channelData,spikeTimes)
end
