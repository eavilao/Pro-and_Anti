function [hMainFig] = markStablePeriods(fileList,eyeChannels)
%eric stability fucntion
%TODO: help
%INPUTS: fileList - cell with row for each file,
%                   1st column filename (inc dir and ext)

%TODO: should load all files and spike sorting results
%TODO cycle through units and ask to mark stable section ends (only 1 for now)

%INITIALISE FIGURE
hMainFig = figure('units','normalized',...
    'position',[0.1 0.1 0.8 0.8],'CloseRequestFcn',{@closeMainWindow},...
    'toolbar','figure','tag','tMainFig');

setappdata(0,'hMainFig',hMainFig)

%main axis for signal
hMainPlot = axes('Parent',hMainFig,'units','normalized',...
    'position',[0.05 0.55 0.8 0.4],'tag','tMainPlot');

hMarkerPlot = axes('Parent',hMainFig,'units','normalized',...
    'position',[0.05 0.05 0.8 0.4],'tag','tMarkerPlot');

%hEyePlot = axes('Parent',hMainFig,'units','normalized',...
%    'position',[0.05 0.05 0.8 0.1],'tag','tEyePlot');
%TODO

set(get(hMarkerPlot,'xLabel'),'string','time (s)')
set(get(hMarkerPlot,'yLabel'),'string','Fr Modulation (z-stat)')

linkaxes([hMainPlot hMarkerPlot],'x')

numFiles = size(fileList,1);
structInitializer = cell(numFiles,1);
dataStruct = struct('spkData',structInitializer,'spkFs',...
    structInitializer,'spkDataStartTime',structInitializer,...
    'spkTimeVec',structInitializer,'alignedSpikes',structInitializer,...
    'spikeDensities',structInitializer);

cc = lines(12); %colour order for units and files

%DATA

%Create 1 long plot of all data using file starting times to check
%continous, colour each file differently

%loop over files and load data
for fileNum = 1:numFiles
    spkChannelNum = fileList(fileNum).spkChannelNumber;
    
    %build variable strings
    
    %because some files use 'Spk' and some use 'SPK' let's check first
    spkCheck = whos('-file',fileList(fileNum).fileName,'-regexp','Spk');
    if numel(spkCheck)<1
        channelType = 'SPK';
    else
        channelType = 'Spk';
    end
    dataString = ['C' channelType '_' sprintf('%03d',spkChannelNum)];
    dataFsString = [dataString '_KHz'];
    dataStartTimeString = [dataString '_TimeBegin'];
    
    %build variable list to load
    %TODO need to use eyeChannels here to build correct variables to load
    
    %     varList = {dataString,dataFsString,dataStartTimeString,...
    %         'CAI_017_TimeBegin','CAI_017_KHz','CAI_017','CAI_018'};
    
    eyeChannelVarNameX = ['CAI_' sprintf('%03d',eyeChannels(1))];
    eyeChannelVarNameY = ['CAI_' sprintf('%03d',eyeChannels(2))];
    
    varList = {dataString,dataFsString,dataStartTimeString,...
        eyeChannelVarNameX, eyeChannelVarNameY, [eyeChannelVarNameX '_KHz'],...
        [eyeChannelVarNameX '_TimeBegin']};
    %load data
    
    load(fileList(fileNum).fileName,varList{:})
    
    %TODO replace with dynamic field names?
    %some file use upper case letters so try first time and then switch to
    %upper
    %     try
    %     eval(['spkData = double(' dataString ');']);
    %     catch
    %         dataString = upper(dataString);
    %         dataFsString = [dataString '_KHz'];
    %         dataStartTimeString = [dataString '_TimeBegin'];
    %          varList = {dataString,dataFsString,dataStartTimeString,...
    %         'CAI_017_TimeBegin','CAI_017_KHz','CAI_017','CAI_018'};
    %         load(fileList(fileNum).fileName,varList{:})
    %         eval(['spkData = double(' dataString ');']);
    %     end
    eval(['spkData = double(' dataString ');']);
    eval(['spkFs = ' dataFsString '*1000;']);
    eval(['spkDataStartTime = ' dataStartTimeString ';']);
    dataStruct(fileNum).spkData = spkData;
    dataStruct(fileNum).spkFs = spkFs;
    dataStruct(fileNum).spkDataStartTime = spkDataStartTime; %time in s since system on
    
    
    
    
    
    %load eye data and flip upside down
     eval(['eyeFs = double(' [eyeChannelVarNameX '_KHz'] ');']);
      eval(['eyeStartTime = double(' [eyeChannelVarNameX '_TimeBegin'] ');']);
       eval(['eyeX = double(' eyeChannelVarNameX ');']);
        eval(['eyeY = double(' eyeChannelVarNameY ');']);
    dataStruct(fileNum).eyeFs = eyeFs*1000;
    dataStruct(fileNum).eyeX = -eyeX;
    dataStruct(fileNum).eyeY = -eyeY;
    dataStruct(fileNum).eyeT = linspace(eyeStartTime,eyeStartTime+length(dataStruct(fileNum).eyeX)/dataStruct(fileNum).eyeFs,length(dataStruct(fileNum).eyeX));
    
    
    
    
    
    endTime = spkDataStartTime+length(spkData)/spkFs;
    tAxis = linspace(spkDataStartTime,endTime,length(spkData));
    
    %plot on grapth
    line(tAxis,spkData,'parent',hMainPlot,'tag',['tSignal ' num2str(fileNum)],'color',cc(fileNum,:))
    dataStruct(fileNum).spkTimeVec = tAxis;
    
    
    %load spikes
    %need to determine if to load the original spikes (output of Beerend's
    %sorter) or if teh resortFlag for this file is 1 then should load the
    %resorted directory
    unitNumbers = fileList(fileNum).unitNumbers; %TODO this will ignore additional units added during resorting
    alignedSpikes = cell(1,max(unitNumbers));
    if fileList(fileNum).resortFlag
        
        %build the file name and load
        resortName = [fileList(fileNum).fileName(1:end-4) '-editedSpikeTimes.mat'];
        inStruct = load(resortName);
        for unitNum = unitNumbers
            unitSpikes = inStruct.spikeTimes{unitNum};
            unitSpikes = unitSpikes(:);
            alignedSpikes{unitNum} =unitSpikes+spkDataStartTime;
        end
        clear inStruct
    else
        inStruct = load([fileList(fileNum).fileName(1:end-4) '_0' num2str(spkChannelNum) '.spi'],'-mat');
        
        
        
        %alignedSpikes = cell(1,length(unitNumbers)); %removed to try and deal
        %with units 4 and 5
        
        for unitNum = unitNumbers
            
            unitSpikes = inStruct.S.Tm(inStruct.S.NetId==unitNum,6);
            
            %take into account start of file
            unitSpikes = unitSpikes+spkDataStartTime; %TODO check if eric always used whole file or if teh infromation is in the output structure
            alignedSpikes{unitNum} = unitSpikes;
            
        end
        %clear F L Q R S
        clear inStruct
    end
    dataStruct(fileNum).alignedSpikes = alignedSpikes;
    
    
    clear(varList{:})
end

timeLims = [min([dataStruct.spkTimeVec]),max([dataStruct.spkTimeVec])];
set(hMainPlot,'xlim',timeLims)

%SORTED SPIKES
%Calculate sliding window firing rate using ksdensity for each unit and
%draw lines on graph
%
binSize = 1; %for calculating spike density
%collect all spike times and reshape
allSpikes = {dataStruct.alignedSpikes};
allSpikes = vertcat(allSpikes{:});

unitNumbers = unique([fileList.unitNumbers]);%unit numbers used across al files
allSpikeDensities = cell(1,max(unitNumbers));
for unitNum = unitNumbers
    allUnitSpikes = vertcat(allSpikes{:,unitNum});
    if ~isempty(allUnitSpikes)
        %calculate spike density
        spDense = ksdensity(allUnitSpikes,[timeLims(1):1:timeLims(2)],'width',binSize);
        allSpikeDensities{unitNum} = spDense*length(allUnitSpikes); %convert to Hz
        allSpikeDensities{unitNum} = (allSpikeDensities{unitNum}-nanmean(allSpikeDensities{unitNum}(:)))./nanstd(allSpikeDensities{unitNum}); %try to make mean 0 to see both traces
        %plot line on marker axis
        line([timeLims(1):1:timeLims(2)],allSpikeDensities{unitNum},'color',cc(unitNum,:),'parent',hMarkerPlot)
        
        
        %line for spike markers at mean spike density
        lineHeight = mean(allSpikeDensities{unitNum});
        hSpikeMarkers = line(allUnitSpikes,lineHeight*ones(length(allUnitSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
            'linestyle','none','marker','s','MarkerFaceColor',cc(unitNum,:),'tag',['tSpikeMarkers ' num2str(unitNum)]);
        
    end
end


dataStruct(fileNum).spikeDensities = allSpikeDensities; %TODO only held for final file at the moment






%mark in graph

%BUTTONS:
%mark period
hMarkLaunch = uicontrol('style','push',...
    'Parent',hMainFig,'units','normalized',...
    'position',[0.875 0.125 0.11 0.05],...
    'string','Mark Stable Period','Callback',{@markLaunch,1,timeLims});

%delete last section
hDeleteMarks = uicontrol('style','push',...
    'Parent',hMainFig,'units','normalized',...
    'position',[0.875 0.075 0.11 0.05],...
    'string','Delete Marks','Callback',{@deleteMarks});

%save and exit
hSaveAndExit = uicontrol('style','push',...
    'Parent',hMainFig,'units','normalized',...
    'position',[0.875 0.025 0.11 0.05],...
    'string','Save and Exit','Callback',{@saveAndExit,dataStruct,fileList});
%mark all
hMarkAll = uicontrol('style','push',...
    'Parent',hMainFig,'units','normalized',...
    'position',[0.875 0.175 0.11 0.05],...
    'string','Mark All','Callback',{@markLaunch,2,timeLims});
end


function [] = markLaunch(src,eventData,optionArg,optionLims)
disp('markLaunch')

%create lines at points marked in graph
%TODO should be draggable
%mark artefact period using two impoints, add to table and draw patch
hMainPlot = findobj('tag','tMainPlot');
hMarkerPlot = findobj('tag','tMarkerPlot');
%create two impoints to select artefact periods and get positions
switch optionArg
    case 1
        hFirstSel = impoint(hMainPlot);
        firstSelPos = getPosition(hFirstSel)
        hSecondSel = impoint(hMainPlot);
        secondSelPos = getPosition(hSecondSel)
        delete([hFirstSel hSecondSel])
        timeLimits = sort([firstSelPos(1) secondSelPos(1)]);
    case 2
        
        timeLimits = optionLims;
end
%get time limits


%TODO error checking to check within file range

%create patch in both axes to display marked points
curYLims = get(hMainPlot,'Ylim');
%draw patch around points
patch([timeLimits(1) timeLimits(1) timeLimits(2) timeLimits(2)],[curYLims(1) curYLims(2) curYLims(2) curYLims(1)],...
    [1 0 0],'FaceAlpha',0.5,'Parent',hMainPlot,'FaceColor','g',...
    'tag','tStablePeriod');
%now marker plot
curYLims = get(hMarkerPlot,'Ylim');
patch([timeLimits(1) timeLimits(1) timeLimits(2) timeLimits(2)],[curYLims(1) curYLims(2) curYLims(2) curYLims(1)],...
    [1 0 0],'FaceAlpha',0.5,'Parent',hMarkerPlot,'FaceColor','g',...
    'tag','tStablePeriod');


end

function [] = deleteMarks(varargin)
disp('deleteMarks')
delete(findobj('tag','tStablePeriod'))
end

function [] = saveAndExit(buttnHandle,dummy,dataStructure,fileList)
disp('saveAndExit')
%TODO warning
hMainPlot = findobj('tag','tMainPlot');
hPatch = findobj(hMainPlot,'tag','tStablePeriod');

if ~isempty(hPatch)
    patchLimits = get(hPatch,'XData');
    timeLimits = [patchLimits(1) patchLimits(3)];
    %startOfFirstFile = dataStructure(1).spkTimeVec(1);
    %relativeTimeLimits = timeLimits - startOfFirstFile
    %fill in start and stop times across file list
    %Codes: [] - not analysed
    %       NaN - none to be used
    %       Inf from beginning/to end of file
    %       any number is time in s to be used
    numFiles = size(fileList,1);
    %editedFileList = fileList;
    for fileNum = 1:numFiles
        fileStart = dataStructure(fileNum).spkTimeVec(1);
        fileEnd = dataStructure(fileNum).spkTimeVec(end);
        
        if fileEnd < timeLimits(1) || fileStart > timeLimits(2)%file before stable period or beyond end
            fileList(fileNum).startTime = NaN;
            fileList(fileNum).stopTime = NaN;
        elseif fileEnd <= timeLimits(2) && fileStart >= timeLimits(1)%file totally within period
            fileList(fileNum).startTime = inf;
            fileList(fileNum).stopTime = inf;
        elseif fileEnd < timeLimits(2) && fileStart <= timeLimits(1)%start midway through file
            fileList(fileNum).startTime = timeLimits(1)-fileStart;
            fileList(fileNum).stopTime = inf;
        elseif fileEnd >= timeLimits(2) && fileStart > timeLimits(1)%end midway through file
            fileList(fileNum).startTime = inf ;
            fileList(fileNum).stopTime = timeLimits(2)-fileStart;
        elseif fileEnd >= timeLimits(2) && fileStart <= timeLimits(1)%start and end within file
            fileList(fileNum).startTime = timeLimits(1)-fileStart ;
            fileList(fileNum).stopTime = timeLimits(2)-fileStart;
        end
    end
    %assign variable to baseworkspace %TODO-needs to be caller function?
else
    numFiles = size(fileList,1);
    newVals = repmat({inf},1,numFiles);
    [fileList.startTime] = newVals{:};
    [fileList.stopTime] =  newVals{:};
    %TODO should warn that selecting all and mark [] if not agree
end
assignin('base','dataStructure',dataStructure)
assignin('base','editedFileList',fileList)

closeMainWindow
end

function [] = closeMainWindow(varargin)
disp('crf!!!!')
uiresume
delete(gcf)
end
