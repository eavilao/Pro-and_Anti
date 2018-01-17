function [hSpikeGui] = launchSpikeGui(channelData,spikeTimes,varargin)
%gui for plotting data and spike markers and enable their deletion or
%reassignment to another unit 2
% p.holland@erasmusmc.nl
%INPUTS:
%channelData.rawSignal - vecotr of raw signal
%                                       .fs - sampling frequency (Hz)
%                                       .filtSignal - filteredSignal
%                                       .lpFiltFreq - low pass freq in Hz
%                                       .hpFiltFreq - high pass freq in Hz
%                                       .timesToInclude - Currently unused so set to [-inf inf]
%                                       .timesToExclude - table with times to exclude form further analysis (periodNum,[start end])
%                                       .fullTimeVector - vector of same
%                                       size as data giving time at each
%                                       point
%   spikeTimes - vector of 1 unit of spike times in s
%       if is cell array then each cell should be vctor of spike times for
%       a different unit
%   param/value pairs
%   'otherEvents' - cell array - different event times in seconds, format
%       is the same as for spikeTimes
%OUTPUTS:
%   assigns, editedSpikeTimes - cell array in same format as input
%   handle to figure

%delete(findobj('tag','tMainFig'))

%TODO filter plot toggle doesn't switch off when plot is killed
%PH edits in 2015 for eric
%1 TODO need to be able to jmp through spikes not just time
%2  jump to a spike by clicking on it in a raster, this works (inst rast
%button) but currently the raster doesn't update itself unless you close
%and reopen
%% TAGS
%'tZoomEdit' - edit box controlling time step size
%'tFullPlot' - main plotting axis
%'tMarkerPlot' - axis controlling the spike markers
%'tSliderPlot' - axis with full signal displayed and slider
%'tUnitTargetForSel' - list box for unit of selected for viewing spike
%'tSpikeMaxSel' - editbox (diabled) to show number of spikes on this unit
%'tSpikeIndSel' - edit box (enabled) to show which spike ind to centre on

%cPalette.snrCmap = [0.2 0.8 0;0.8 0.8 0;0 0.8 0.8; 1 0 0;0 0 0]; %TODO add more colours
%settings
startDispTime = 5; %number of seconds to start displaying of signal
%% INPUTS

%TODO input parser and run silent
% if nargin>1
%     spikeTimes = varargin{1};
% else
%     spikeTimes = [];
% end


% otherEvents = [];
% if nargin>3
% otherEvents = varargin{2};
% end

p = inputParser();
p.addParamValue('multipleFileInfo', [], @(x) iscell(x)); %cell array containng namesand length of individual sections in trace
p.addParamValue('otherEvents', [], @(x) iscell(x)); %marked as patches
p.addParamValue('colourScheme', 1, @(x) isnumeric(x)); %to swtich colour schemes
p.parse(varargin{:});
multipleFileInfo = p.Results.multipleFileInfo;
otherEvents = p.Results.otherEvents;
colourScheme = p.Results.colourScheme;

if colourScheme==1
    %colours:
    cPalette.subPanelbg = [0.8 0.8 0.8];
    cPalette.bg = [0.9 0.9 0.9]; %main background colour
    cPalette.axBg = [1 1 1];
    cPalette.spikeColours = lines(12); %TODO SNR colour scheme?
    cPalette.sigCol = [0 0 0];
else
    cPalette.subPanelbg = [0.5 0.5 0.5];
    cPalette.bg = [0.2 0.2 0.2]; %main background colour
    cPalette.axBg = [0.3 0.3 0.3];
    cPalette.spikeColours = lines(12); %TODO SNR colour scheme?
    cPalette.sigCol = [0.1 0.8 0.8];
end

%use filtered if available
if ~isempty(channelData.filtSignal)
    data = channelData.filtSignal(:); %force column
else
    data = channelData.rawSignal(:); %force column
end

fs = channelData.fs;

%now nan all the artefact periods
if ~isempty(channelData.timesToExclude)
    numPeriods = size(channelData.timesToExclude,1);
    for periodNum = 1:numPeriods
        
        thisPeriodStartSample = floor(channelData.timesToExclude(periodNum,1)*fs);
        thisPeriodStopSample = ceil(channelData.timesToExclude(periodNum,2)*fs);
        data(thisPeriodStartSample:thisPeriodStopSample) = nan;
    end
end

processedData.data = data;
processedData.fs = fs;
channelData.multipleFileInfo = multipleFileInfo; %for use later when saving
%% INITIALISE FIGURE
hSpikeGui = figure('units','normalized',...
    'position',[0.1 0.1 0.8 0.8],'CloseRequestFcn',{@closeMainWindow},...
    'toolbar','figure','tag','tMainFig','KeyPressFcn', {@keyPress},'color',cPalette.bg); %,...
%TODO close request function asking if to assignout or not
setappdata(0,'hSpikeGui',hSpikeGui)

set(gcf,'renderer','zBuffer')
%main axis for signal
hFullPlot = axes('Parent',hSpikeGui,'units','normalized',...
    'position',[0.05 0.2 0.795 0.65],'tag','tFullPlot','color',cPalette.axBg);
%main axis for signal
hMarkerPlot = axes('Parent',hSpikeGui,'units','normalized',...
    'position',[0.05 0.85 0.795 0.1],'tag','tMarkerPlot',...
    'ButtonDownFcn',{@clickMarker},'color',cPalette.axBg); %
%full signal to enable zoom control of main plot
hSliderPlot = axes('Parent',hSpikeGui,'units','normalized',...
    'position',[0.05 0.05 0.8 0.1],'tag','tSliderPlot','color',cPalette.axBg);


set(get(hFullPlot,'xLabel'),'string','Time (s)')
set(get(hFullPlot,'yLabel'),'string','Amplitude (\muV)')
set(get(hMarkerPlot,'yLabel'),'string','Unit')
set(get(hMarkerPlot,'xLabel'),'string','')
set(hMarkerPlot,'xTick',[])

%plot signal
%timeVector = linspace(0,length(data)/fs,length(data));
%now get logical to show which times to display
timeToDisplay = channelData.fullTimeVector>0 & channelData.fullTimeVector<startDispTime;
hSignal = line(channelData.fullTimeVector(timeToDisplay),...
    data(timeToDisplay),'color',cPalette.sigCol,'tag','tSignal','Parent',hFullPlot);
initialYLims = 1.5*[nanmin(data) nanmax(data)]; %Y axis limits for main plot
set(hFullPlot,'ylim',initialYLims);

%plot downsampled version of signal in the slider plot
dsFactor = 20;
dsData = data(1:dsFactor:end);
dsTimeVetor = linspace(0,length(data)/fs,length(dsData));

hDsSig = line(dsTimeVetor,...
    dsData,'color',cPalette.sigCol,'Parent',hSliderPlot);
set(hSliderPlot,'xlim',[0 length(data)/fs])
set(hSliderPlot,'ylim',initialYLims)
% now create 2 lines and a patch to act as the slider %TODO draggable
% slider edges
% hVcursor1 = line([0 0],2*[min(data) max(data)],'Color',[0 1 0],...
%     'LineWidth',3,'Parent',hSliderPlot,'yliminclude','off');
% hVcursor2 = line([startDispTime startDispTime],2*[min(data) max(data)],'Color',[1 1 0],...
%     'LineWidth',3,'Parent',hSliderPlot,'yliminclude','off');

%draw patch for slider main button between vertical cursors
hTimeSlider = patch([startDispTime startDispTime 0 0],3*[initialYLims(1) initialYLims(2) initialYLims(2) initialYLims(1)],...
    [0.8 0.8 0.8],'FaceAlpha',0.5,'tag','tTimeSliderBar','hittest','on',...
    'yliminclude','off','parent',hSliderPlot);

set(hTimeSlider,'ButtonDownFcn',{@clickTimeSlider,hSpikeGui})
if ~isempty(spikeTimes)
    %process spikes, switch for cell or vector input
    if iscell(spikeTimes)
        % disp('currently only vector input')
        numUnits = numel(spikeTimes); %TODO check if correct and consistent
        
    elseif isvector(spikeTimes)
        
        spikeTimes = {spikeTimes(:)}; %force column
        numUnits = 1;
        
    end
    cc = cPalette.spikeColours;
    %plot line
    for unitNum = 1:numUnits
        lineHeight = unitNum;
        unitSpikes = spikeTimes{unitNum};
        line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
            'linestyle','none','marker','s','MarkerFaceColor',cPalette.spikeColours(unitNum,:),...
            'tag',['tUnit ' num2str(unitNum)]);
    end
    set(hMarkerPlot,'ylim',[0 numUnits+1]);
    set(hMarkerPlot,'yTick',[1:numUnits])
else
    set(hMarkerPlot,'ylim',[0 1],'yTick',[0 1]);
    numUnits = 0;
end
linkaxes([hFullPlot hMarkerPlot],'x')

set(hFullPlot,'xlim',[0 startDispTime]) %TODO shoudl be replaced with setting xdata

%plot other evnts if they exist
if ~isempty(otherEvents)
    numEventTypes = size(otherEvents,2); %TODO should rearrange to find max size
    for eventTypeNum = 1:numEventTypes
        lineHeight = eventTypeNum+0.5;
        eventTimes = otherEvents{eventTypeNum};
        line(eventTimes,lineHeight*ones(length(eventTimes),1),'parent',hMarkerPlot,'color',cc(12-eventTypeNum,:),...
            'linestyle','none','marker','d','MarkerFaceColor',cPalette.spikeColours(12-eventTypeNum,:),...
            'tag',['tEventLine ' num2str(eventTypeNum)]);
    end
end
%% BUTTONS
hAnalysisButPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.725 0.11 0.15],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Analysis')
%button launch interactive waveform pciker
hFir = uicontrol('style','push',...
    'Parent',hAnalysisButPan,'units','normalized',...
    'position',[0.05 0.85 0.9 0.2],...
    'string','Int Wave','Callback',{@launchInteractiveWaveView,0});
%button launch interactive raster
hFir = uicontrol('style','push',...
    'Parent',hAnalysisButPan,'units','normalized',...
    'position',[0.05 0.65 0.9 0.2],...
    'string','Int Rast','Callback',{@launchInteractiveRaster});
%button to launch firing rate figure
hFr = uicontrol('style','push',...
    'Parent',hAnalysisButPan,'units','normalized',...
    'position',[0.05 0.45 0.9 0.2],...
    'string','Firing Rates','Callback',{@calculateFiringRates});
%button to launch isi figure
hISI = uicontrol('style','push',...
    'Parent',hAnalysisButPan,'units','normalized',...
    'position',[0.05 0.25 0.45 0.2],...
    'string','ISIs','Callback',{@createIsiFig});
%button to launch spike summary figure
hSpikeSum = uicontrol('style','push',...
    'Parent',hAnalysisButPan,'units','normalized',...
    'position',[0.5 0.25 0.45 0.2],...
    'string','Spike Summary','Callback',{@createSpikeSumFig});
%button to launch spike summary figure
hSpikeAlign = uicontrol('style','push',...
    'Parent',hAnalysisButPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.2],...
    'string','Spike Align','Callback',{@redoSpikeAligment});
%period to analyse selection + tick box for all
stopTimesS = 30;
startTimesS = 0;
hAnSelButPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.675 0.11 0.05],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Analysis Period');
%auto scale tick box
hAnWholeFile = uicontrol('style','checkbox',...
    'Parent',hAnSelButPan,'units','normalized',...
    'position',[0.35 0.05 0.3 0.9],...
    'tag','tAnalyseWholeFileTb',...
    'tooltipstring','Whole File','value',1,... %TODO make starting default off?
    'callback',{@anWholeFile});

%
%Start Period edit
hAnTimes.Start = uicontrol('style','edit',...
    'Parent',hAnSelButPan,'units','normalized',...
    'position',[0.05 0.05 0.3 0.9],...
    'tag','tAnStartTime','String',num2str(startTimesS),...
    'TooltipString','Analysis Start (s)');
%Stop period edit
hAnTimes.Stop = uicontrol('style','edit',...
    'Parent',hAnSelButPan,'units','normalized',...
    'position',[0.65 0.05 0.3 0.9],...
    'tag','tAnStopTime','String',num2str(stopTimesS),...
    'TooltipString','Analysis Stop (s)');
%'callback',{@anManTime,hFullPlot,hAnTimes},...


set(hAnTimes.Start,'callback',{@anManTime,hAnWholeFile,hAnTimes})
set(hAnTimes.Stop,'callback',{@anManTime,hAnWholeFile,hAnTimes})

%now draw patch on main plot and slider plot
hAnPatch.Full =  patch([stopTimesS stopTimesS startTimesS startTimesS],3*[initialYLims(1) initialYLims(2) initialYLims(2) initialYLims(1)],...
    [0.2 0.8 0.2],'FaceAlpha',0.5,'tag','tAnPatch',... %'hittest','on',...
    'yliminclude','off','parent',hFullPlot);

%draw patch for slider
%now draw analysis period patch on main plot and slider plot
hAnPatch.Slid =  patch([stopTimesS stopTimesS startTimesS startTimesS],3*[initialYLims(1) initialYLims(2) initialYLims(2) initialYLims(1)],...
    [-1 -1 -1 -1],...     %3d is to ensure is displayed at back and therefor slider remains clickable
    [0.2 0.8 0.2],'FaceAlpha',0.5,'tag','tAnPatchSlid',... %'hittest','on',...
    'yliminclude','off','parent',hSliderPlot);

setappdata(hSpikeGui,'hAnPatch',hAnPatch) %TODO start of repalcement of findobj commands with ahndles in appdata

if get(hAnWholeFile,'value')
    
    set(hAnPatch.Slid,'visible','off')
    set(hAnPatch.Full,'visible','off')
end
%button for selecting spike as number
hSpikeSelectedPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.525 0.11 0.1],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Spike Selected');
%unit selector
%list to select which unit to move to
if numUnits>0
    a = 1:numUnits;
    listStrings = cellstr(num2str(a(:)));
    boxStringSel = num2str(1); %initial number in the edit boxes
    boxStringMax = length(spikeTimes{1});
else
    %TODO the calbacks to these functions have no error checking for no
    %spikes
    listStrings = {''};
    boxStringSel = '';
    boxStringMax = '';
end
hUnitSel = uicontrol('style','popupmenu',...
    'Parent',hSpikeSelectedPan,'units','normalized',...
    'position',[0.35 0.75 0.3 0.2],...
    'string',listStrings,'value',1,... %{'1','2','3','4','+'}
    'TooltipString','Select Unit, Right click to update max spikes',...
    'tag','tUnitTargetForSel');
%select spike by manual text entry box
hSpikeIndSel = uicontrol('style','edit',...
    'Parent',hSpikeSelectedPan,'units','normalized',...
    'position',[0.15 0.05 0.3 0.4],...
    'tag','tSpikeIndSel','String',boxStringSel,...
    'callback',{@centreOnSpike,hUnitSel},...
    'TooltipString','Spike Ind');
% box to say maximum number of spikes for this unit %TODO should update
%when each unit selected
hSpikeMaxSel = uicontrol('style','edit',...
    'Parent',hSpikeSelectedPan,'units','normalized',...
    'position',[0.55 0.05 0.3 0.4],...
    'enable','off',...
    'tag','tSpikeMaxSel','String',boxStringMax,...
    'TooltipString','Max Spike Ind');
%now set the callback to the unit selector and the edit box so that they
%have the same input order
set(hUnitSel,'callback',{@centreOnSpike,hUnitSel,hSpikeIndSel,0},...
    'ButtonDownFcn',{@updateMaxSpikeInd,hSpikeMaxSel})
set(hSpikeIndSel,'callback',{@centreOnSpike,hUnitSel,hSpikeIndSel,0})
%now we have the handles to the selectors we can now pass them to the left
%and right buttons
hSpikeStepForward = uicontrol(hSpikeSelectedPan,'style','push',...
    'units','normalized','position',[0.65 0.65 0.3 0.3],...
    'string','->','tag','zoomUpSmall','callback',{@centreOnSpike,hUnitSel,hSpikeIndSel,1});
hSpikeStepBackward = uicontrol(hSpikeSelectedPan,'style','push',...
    'units','normalized','position',[0.05 0.65 0.3 0.3],...
    'string','<-','tag','zoomDownSmall','callback',{@centreOnSpike,hUnitSel,hSpikeIndSel,-1});

%buttons for editing spikes
hSpikeButPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.325 0.11 0.2],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Spike Control');
%button to launch spike marker creation
hNewSpike = uicontrol('style','push',...
    'Parent',hSpikeButPan,'units','normalized',...
    'position',[0.05 0.75 0.9 0.2],...
    'string','New Spike','Callback',{@markNewSpike});
%button to delete spikes
hDeleteSpikes = uicontrol('style','push',...
    'Parent',hSpikeButPan,'units','normalized',...
    'position',[0.05 0.25 0.9 0.2],...
    'string','Delete Spikes','Callback',{@deleteSpikes});


%button to move to otehr unit
hMoveSpikes = uicontrol('style','push',...
    'Parent',hSpikeButPan,'units','normalized',...
    'position',[0.05 0.45 0.9 0.2],...
    'string','Move Spikes','Callback',{@moveSpikes});
%list to select which unit to move to
if numUnits>0
    a = 1:numUnits;
    listStrings = cellstr(num2str(a(:)));
    listStrings = vertcat(listStrings,'+');
else
    listStrings = {'+'};
end
hMoveSpikesTo = uicontrol('style','popupmenu',...
    'Parent',hSpikeButPan,'units','normalized',...
    'position',[0.05 0.65 0.9 0.1],...
    'string',listStrings,'value',1,... %{'1','2','3','4','+'}
    'tag','tUnitTarget');
%clear selection button
hClSpikes = uicontrol('style','push',...
    'Parent',hSpikeButPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.2],...
    'string','Clear Selection','Callback',{@clSelSpikes});

%buttons for controlling automatic spike sorting
hSortSpikesButPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.875 0.11 0.1],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Auto Sorting');
%launch auto sorting
hStartSorting = uicontrol('style','push',...
    'Parent',hSortSpikesButPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.45],...
    'string','Auto-Sort','Callback',{@launchSorter});
%launch template match
hStartTempMatch = uicontrol('style','push',...
    'Parent',hSortSpikesButPan,'units','normalized',...
    'position',[0.05 0.5 0.9 0.45],...
    'string','Template Match','Callback',{@launchTemplateMatch});

%buttons for importing or exporting spikes
hImpExpButPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.025 0.11 0.1],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Import/Export');
%save and exit
hSaveAndExit = uicontrol('style','push',...
    'Parent',hImpExpButPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.45],...
    'string','Save and Exit','Callback',{@saveAndExit});
%save and exit
hImportSpikes = uicontrol('style','push',...
    'Parent',hImpExpButPan,'units','normalized',...
    'position',[0.05 0.5 0.9 0.45],...
    'string','Import','Callback',{@importSpikes});


% amplitude axis scale control
hYaxButPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.225 0.11 0.05],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Amp Axis Control')
%auto scale tick box
hYautoScale = uicontrol('style','checkbox',...
    'Parent',hYaxButPan,'units','normalized',...
    'position',[0.35 0.05 0.3 0.9],...
    'tag','tAutoYScaleTb',...
    'string','Auto Scale','value',0,...
    'callback',{@ampAxisControl,hFullPlot});


%ymin edit box
hYaxLimLow = uicontrol('style','edit',...
    'Parent',hYaxButPan,'units','normalized',...
    'position',[0.05 0.05 0.3 0.9],...
    'tag','tYaxLimLow','String',num2str(initialYLims(1)),...
    'callback',{@ampAxisMan,hFullPlot},...
    'TooltipString','Y lim min');
%ymax edit box
hYaxLimHigh = uicontrol('style','edit',...
    'Parent',hYaxButPan,'units','normalized',...
    'position',[0.65 0.05 0.3 0.9],...
    'tag','tYaxLimHigh','String',num2str(initialYLims(2)),...
    'callback',{@ampAxisMan,hFullPlot},...
    'TooltipString','Y lim max');
set(hFullPlot,'ylimmode','manual')


% toggle button to control the plotting of a differently filtered signal
hFiltPlotButPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.275 0.11 0.05],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Filter Plot')
%auto scale tick box
hShowFiltPlot = uicontrol('style','checkbox',...
    'Parent',hFiltPlotButPan,'units','normalized',...
    'position',[0.35 0.05 0.3 0.9],...
    'tag','tLaunchFiltPlot',...
    'string','Auto Scale','value',0,...
    'callback',{@launchFiltPlot});

%threshold based detection, %TODO set manually the detection threshold
hThreshButPan = uipanel('parent',hSpikeGui,'units','normalized',...
    'position',[0.875 0.625 0.11 0.05],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Threshold Crossings');
%tick box for positive or negative crossing  thresh or both or
%either
% hPosNegTheshBg = uibuttongroup('parent',hThreshButPan,'units','normalized',...
%     'position',[0.05 0.05 0.9 0.45],'tag','tPosNegTheshBg');
%for no the lines are linked in the sense that moving one line moves the
%other (even if it is invisible)
startingThreshLevel = [3*nanstd(data) 3*nanstd(data)];
%horizontal cursor
hHcursorPos = line([channelData.fullTimeVector(1) channelData.fullTimeVector(end)],startingThreshLevel,...
    'Color',[1 0 0],'LineWidth',2,'tag','tHorPosLine','Parent',hFullPlot);

%horizontal cursor
hHcursorNeg = line([channelData.fullTimeVector(1) channelData.fullTimeVector(end)],-startingThreshLevel,...
    'Color',[1 0 0],'LineWidth',2,'tag','tHorNegLine','Parent',hFullPlot);
set(hHcursorPos,'ButtonDownFcn',{@clickLine,hSpikeGui,hHcursorNeg})
set(hHcursorNeg,'ButtonDownFcn',{@clickLine,hSpikeGui,hHcursorPos})



hPosThreshTb = uicontrol('parent',hThreshButPan,'style','check',...
    'units','normalized',...
    'position',[0.05 0.1 0.45 0.6],'string','pos',...
    'value',1,'callback',{@turnOnOffLine,hHcursorPos});
hNegThreshTb = uicontrol('parent',hThreshButPan,'style','check',...
    'units','normalized',...
    'position',[0.55 0.1 0.45 0.6],'string','neg',...
    'value',1,'callback',{@turnOnOffLine,hHcursorNeg});
%set(hPosNegTheshBg,'selectedobject',hPosThreshRb)




%zoom buttons and time step buttons
hMainPlotPanel.zoomBg = uibuttongroup(hSpikeGui,'units','normalized',...
    'BackGroundColor',cPalette.subPanelbg,...
    'position',[0.875 0.125 0.11 0.1],'title','time axis');
hMainPlotPanel.zoomIn.pb = uicontrol(hMainPlotPanel.zoomBg,'style','push',...
    'units','normalized','position',[0.2 0.675 0.6 0.3],...
    'string','+','tag','zoomUpBig','callback',{@timeAxisControl,1});
hMainPlotPanel.zoomOut.pb = uicontrol(hMainPlotPanel.zoomBg,'style','push',...
    'units','normalized','position',[0.2 0.025 0.6 0.3],...
    'string','-','tag','zoomDownBig','callback',{@timeAxisControl,2});
hMainPlotPanel.stepForward.pb = uicontrol(hMainPlotPanel.zoomBg,'style','push',...
    'units','normalized','position',[0.675 0.35 0.3 0.3],...
    'string','->','tag','zoomUpSmall','callback',{@timeAxisControl,3});
hMainPlotPanel.stepBackward.pb = uicontrol(hMainPlotPanel.zoomBg,'style','push',...
    'units','normalized','position',[0.025 0.35 0.3 0.3],...
    'string','<-','tag','zoomDownSmall','callback',{@timeAxisControl,4});
hMainPlotPanel.zoomDownAmount.ed = uicontrol(hMainPlotPanel.zoomBg,'style','edit',...
    'units','normalized','position',[0.35 0.35 0.3 0.3],...
    'string','1','tag','tZoomEdit');
%% FINISH PROCESSING
selectedSpikeInds = cell(1,numUnits);
setappdata(hSpikeGui,'spikeTimes',spikeTimes)
setappdata(hSpikeGui,'selectedSpikeInds',selectedSpikeInds)
setappdata(hSpikeGui,'channelData',channelData);
setappdata(hSpikeGui,'processedData',processedData);
setappdata(hSpikeGui,'cPalette',cPalette)


%TODO add all ahndles to this tructrue to stop need for findobj calls

guiHands.hFullPlot = hFullPlot; %handle for the full signal displaying axis
guiHands.hMarkerPlot = hMarkerPlot; %handle for the spike marker axis
guiHands.hAnPatch = hAnPatch; %handles to slider and full p[lot patches representing analysi time
guiHands.hAnWholeFile = hAnWholeFile; %handle to tickbox for anlysing whole file
guiHands.hAnTimes = hAnTimes; %edit boxes with analysis times
guiHands.hPosThreshTb = hPosThreshTb; %threshold tickboxs
guiHands.hNegThreshTb = hNegThreshTb;
guiHands.hHcursorPos = hHcursorPos;
guiHands.hHcursorNeg = hHcursorNeg;
guiHands.hMoveSpikesTo = hMoveSpikesTo; %handle for listbox displaying available units
%spike selection handles
guiHands.hUnitSel = hUnitSel;
guiHands.hSpikeMaxSel = hSpikeMaxSel;
guiHands.hSpikeIndSel = hSpikeIndSel;
setappdata(hSpikeGui,'guiHands',guiHands)
end


%% ANALYSIS PERIOD FUNCTIONS

function redoSpikeAligment(src,eventData)
%callback from spike alignemnt button should get all spikes and then
%realign them based on some options. Needs to open separate figure to
%display results with ana ccept button somewhere.
%TODO, think about conflict handling
bigTimeWin = [1 4]; %how many ms to take in total
spikeTimeWindow = bigTimeWin./2; %[0.5 2]; %where to look for alignment
%I think it has to be exactly double
hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
spikeTimes = getappdata(hSpikeGui,'spikeTimes');
channelData = getappdata(hSpikeGui,'channelData');
fs = channelData.fs; %TODO switch to processed data if below is ok with nans
data = channelData.filtSignal;
%get active units and then get spike shapes for each
activeUnits = find(~cellfun('isempty',spikeTimes));
spikeShapes = cell(1,max(activeUnits));
removedSpikeInds = cell(1,max(activeUnits));
for unitNum = activeUnits
    
    unitSpikeTimes = spikeTimes{unitNum};
    [spikeShapes{unitNum} timeVector removedSpikeInds{unitNum}] = extractSpikeShapes(data,unitSpikeTimes,bigTimeWin,fs);
    
    if ~isempty(removedSpikeInds{unitNum})
        %if removed then also remove from times
        spikeTimes{unitNum}(removedSpikeInds{unitNum}) = [];
    end
end



hAlignGui.hF = figure;
hAlignGui.hA = axes('parent',hAlignGui.hF,'units','normalized','position',[0.05 0.05 0.6 0.9]);

%draw spikes
timeVector = timeVector*1000;
%now plot each unit in the correct colour
hL = cell(max(activeUnits),1); %handles for lines
for unitNum = activeUnits
    if ~isempty(spikeShapes{unitNum})
        hL{unitNum} = line(timeVector(:),spikeShapes{unitNum},'parent',hAlignGui.hA,'color',cc(unitNum,:));
    end
end

hAlignGui.butPan.pan = uipanel('parent',hAlignGui.hF,'units','normalized',...
    'position',[0.7 0.1 0.25 0.5],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Alignment');
%TODO lisbox of aligment options
listStrings = {'MAX','MIN','MAX DIFF','MAX ABS DIFF'};
hAlignGui.butPan.ls =  uicontrol('style','popupmenu',...
    'Parent',hAlignGui.butPan.pan,'units','normalized',...
    'position',[0.05 0.75 0.9 0.2],...
    'string',listStrings,'value',3,... %{'1','2','3','4','+'}
    'TooltipString','Select Alignment method');
%time window option

%button for go
hAlignGui.butPan.goPb = uicontrol('style','push',...
    'Parent',hAlignGui.butPan.pan,'units','normalized',...
    'position',[0.05 0.45 0.9 0.2],...
    'string','Align Spikes','Callback',{@doAlignment,spikeTimeWindow}); %TODO should be handle to edit boxes for time window
%button for accept
hAlignGui.butPan.yesPb = uicontrol('style','push',...
    'Parent',hAlignGui.butPan.pan,'units','normalized',...
    'position',[0.05 0.25 0.9 0.2],...
    'string','Accept Spikes','Callback',{@acceptAlignment});





%add lines for where alignment will happen
%TODO make draggable
hAlignGui.hLvert.low = line([-spikeTimeWindow(1) -spikeTimeWindow(1)],[-1000 1000],'yliminclude','off','parent',hAlignGui.hA,'color','k');
hAlignGui.hLvert.high = line([spikeTimeWindow(2) spikeTimeWindow(2)],[-1000 1000],'yliminclude','off','parent',hAlignGui.hA,'color','k');
hAlignGui.origSpikeShapes = spikeShapes;
hAlignGui.origSpikeTimes = spikeTimes; %make a copy of original spike times
hAlignGui.bigTimeWin = bigTimeWin;
hAlignGui.spikeTimeWindow = spikeTimeWindow;
hAlignGui.hL = hL;
hAlignGui.fs = fs;
setappdata(0,'hAlignGui',hAlignGui);
end

function [] = acceptAlignment(src,eventData)
%function clled from 'accept' pushbutton in alignment gui, should take the
%new spike times and overwrite the main guis spikeTimes apppdata.
hAlignGui = getappdata(0,'hAlignGui');
hSpikeGui = getappdata(0,'hSpikeGui');
guiHands = getappdata(hSpikeGui,'guiHands'); %main gui handles
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
spikeTimes = hAlignGui.origSpikeTimes;

allSelectedSpikeInds = getappdata(hSpikeGui,'selectedSpikeInds'); %need to delete any elected to stop confusion
delete(findobj(guiHands.hMarkerPlot,'-regexp','tag','tSelUnit ')); %delte marks
delete(findobj(guiHands.hMarkerPlot,'-regexp','tag','tUnit '));
activeUnits = find(~cellfun('isempty',spikeTimes));
%now we need to redraw the spike markers
for unitNum = activeUnits
    
    
    
    
    unitSpikes = spikeTimes{unitNum};
    
    
    %delte existing line and redraw
    delete(findobj('tag',['tUnit ' num2str(unitNum)]))
    lineHeight = unitNum;
    line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',guiHands.hMarkerPlot,'color',cc(unitNum,:),...
        'linestyle','none','marker','s','MarkerFaceColor',cc(unitNum,:),...
        'tag',['tUnit ' num2str(unitNum)]);
    
    %delete selected marker lines and blank selectedSpikeInds
    
    allSelectedSpikeInds{unitNum} = [];
end


setappdata(hSpikeGui,'selectedSpikeInds',allSelectedSpikeInds)
setappdata(hSpikeGui,'spikeTimes',spikeTimes);
delete(hAlignGui.hF);
end

function [] = doAlignment(src,eventData,timeWin)
% realign the spikes based on method slected in  lshandle
hAlignGui = getappdata(0,'hAlignGui');
hSpikeGui = getappdata(0,'hSpikeGui');
spikeTimes = hAlignGui.origSpikeTimes;
activeUnits = find(~cellfun('isempty',spikeTimes));
%channelData = getappdata(hSpikeGui,'channelData');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
fs = hAlignGui.fs; %TODO switch to processed data if below is ok with nans
% data = channelData.filtSignal;
% bigTimeWin = hAlignGui.bigTimeWin;

spikeShapes = hAlignGui.origSpikeShapes;

% spikeShapes = cell(1,max(activeUnits));
% for unitNum = activeUnits
%
%     unitSpikeTimes = spikeTimes{unitNum};
% [spikeShapes{unitNum} timeVector] = extractSpikeShapes(data,unitSpikeTimes,bigTimeWin,fs);
%
% end
%temporary figure
% hF = figure;
% hA = axes('parent',hF);
% for unitNum = activeUnits
%     if ~isempty(spikeShapes{unitNum})
%     hL{unitNum} = line(timeVector(:),spikeShapes{unitNum},'parent',hA,'color',cc(unitNum,:));
%     end
% end


spikeSampleTimes = cellfun(@(x) round(x*fs),spikeTimes,'uniformoutput',false);
spikeTimeWindow = hAlignGui.spikeTimeWindow;
settings.postPeakLength =  floor(spikeTimeWindow(2)*fs/1000);
settings.prePeakLength = floor(spikeTimeWindow(1)*fs/1000); %needs to be exactly half of big window
settings.interpolationFactor = 1;
settings.alignMethod = get(hAlignGui.butPan.ls,'value'); %0 none, 1 min, 2 max , 3 mx diff, 4 max abs diff
epochsAligned = cell(1,max(activeUnits));
spikeIndsAligned = cell(1,max(activeUnits));
for unitNum = activeUnits
    thisUnitSpikes = spikeSampleTimes{unitNum};
    thisUnitSpikes = thisUnitSpikes(:)';
    [epochsAligned{unitNum} spikeIndsAligned{unitNum}] = alignSpikesAO(spikeShapes{unitNum}',thisUnitSpikes,settings)
    
end
newTimeVector = 1000*linspace(-settings.prePeakLength/fs,settings.postPeakLength/fs,size(epochsAligned{unitNum},2));
%temporary figure
hF = figure;
hA = axes('parent',hF);
for unitNum = activeUnits
    if ~isempty(spikeShapes{unitNum})
        hL{unitNum} = line(newTimeVector,epochsAligned{unitNum},'parent',hA,'color',cc(unitNum,:));
    end
end

handlesToLines = hAlignGui.hL;
hAx = hAlignGui.hA; %TODO could switch to appdata/handles struct method if I like this tool


delete(vertcat(handlesToLines{:}))

for unitNum = activeUnits
    if ~isempty(spikeShapes{unitNum})
        hL{unitNum} = line(newTimeVector,epochsAligned{unitNum},'parent',hAx,'color',cc(unitNum,:));
    end
end
delete(hF)
%now need to take into account indices corrections and turn into new spike
%times
spikeIndsAligned{1}(1)
newSpikeTimes = cellfun(@(x) x./fs,spikeIndsAligned,'uniformoutput',false);
%newSpikeTimes = cellfun(@(x) x-(settings.prePeakLength/fs),newSpikeTimes,'uniformoutput',false);
%now add new line handles and spike times back to appdata

%also diabel the button to do alignment as it won't wor a second time
set(hAlignGui.butPan.goPb,'enable','off')
hAlignGui.hL = hL;
hAlignGui.origSpikeTimes = newSpikeTimes;
hAlignGui.origSpikeShapes = epochsAligned;
setappdata(0,'hAlignGui',hAlignGui);
end

function [epochsAligned spikeIndsAligned] = alignSpikesAO(epochs,SpikeInds,settings)
%function that performs alignment of spikes using various options
%contained in settings structure.
%input is a matrix with a row for each spike and a column for each dp
%
%takes spikes that are twice as many samples as requested
%
%update of PJHAlignOffline2
%taken staight from cleanAO folder and made nested here so I can edit it
postSamples = settings.postPeakLength;
epochSamples = settings.postPeakLength + settings.prePeakLength+1;
preSamples = settings.prePeakLength;
numSpikes = length(SpikeInds);

%if interpolation factor is greater than 1 then interpolate
if settings.interpolationFactor == 0
    settings.interpolationFactor = 1;
elseif settings.interpolationFactor > 1
    epochs = reshape(([zeros(1,settings.interpolationFactor-1) interp1(settings.interpolationFactor*[1:prod(size(epochs))],reshape(epochs',1,prod(size(epochs))),[settings.interpolationFactor:settings.interpolationFactor*prod(size(epochs))])]),size(epochs,2)*settings.interpolationFactor,size(epochs,1))';
else
end

%create blank matrix
epochsAligned = nan(numSpikes,epochSamples+1);

%find alignment indexs for each detected spike (finds start)
if settings.alignMethod == 0 %no alignment
    Indexs = settings.interpolationFactor*settings.prePeakLength*ones(numSpikes,1);
elseif settings.alignMethod == 1 % simple minimum alignment
    [C Indexs] = min(epochs(:,settings.interpolationFactor*settings.prePeakLength+1:settings.interpolationFactor*(2*epochSamples-postSamples)-1),[],2);
elseif settings.alignMethod == 2 %max
    [C Indexs] = max(epochs(:,settings.interpolationFactor*settings.prePeakLength+1:settings.interpolationFactor*(2*epochSamples-postSamples)-1),[],2);
elseif settings.alignMethod == 3 %max slope
    [C Indexs] = max(diff(epochs(:,settings.interpolationFactor*settings.prePeakLength+1:settings.interpolationFactor*(2*epochSamples-postSamples)-1),[],2),[],2);
elseif settings.alignMethod == 4 %max abs slope
    [C Indexs] = max(abs(diff(epochs(:,settings.interpolationFactor*settings.prePeakLength+1:settings.interpolationFactor*(2*epochSamples-postSamples)-1),[],2)),[],2);
else
    disp('Warning: Invalid alignment method choice, using no alignement')
    Indexs = settings.interpolationFactor*settings.prePeakLength*ones(numSpikes,1);
end


%build aligned epochs

for spikeNum = 1:numSpikes
    epochsAligned(spikeNum,:) = epochs(spikeNum,Indexs(spikeNum):settings.interpolationFactor:Indexs(spikeNum)+settings.interpolationFactor*epochSamples);
end

%update aligned spike indexs
%spikeIndsAligned = SpikeInds+round((Indexs'-settings.prePeakLength)/(settings.interpolationFactor));
spikeIndsAligned = SpikeInds+round(Indexs'/settings.interpolationFactor)-(settings.prePeakLength);

%remove spikes that have been aligned exactly on top of each other
%[uniqueSpikes m n] = unique(spikeIndsAligned,'first'); %removed for
%compatability with r14sp3
[uniqueSpikes m n] = unique(spikeIndsAligned);
epochsAligned = epochsAligned(m,:);
spikeIndsAligned = uniqueSpikes;
end
%%
function [] = anWholeFile(hTb,eventData)
%callback from tick box to toggle analysis of whole file, also should turn
%on or off patches
hSpikeGui = getappdata(0,'hSpikeGui');
guiHands = getappdata(hSpikeGui,'guiHands');
hAnPatch = guiHands.hAnPatch;



switch get(hTb,'value')
    case 1
        set(hAnPatch.Full,'visible','off')
        set(hAnPatch.Slid,'visible','off')
    case 0
        set(hAnPatch.Full,'visible','on')
        set(hAnPatch.Slid,'visible','on')
end
end

function [] = anManTime(hEdit,eventData,hAnWholeFile,hAnTimes)
%callback when the edit boxes specifying start and end of analysis period are changed
hSpikeGui = getappdata(0,'hSpikeGui');
%get new values form strings

guiHands = getappdata(hSpikeGui,'guiHands');
hAnPatch = guiHands.hAnPatch;

stopTimesStr = get(hAnTimes.Stop,'string');
startTimesStr = get(hAnTimes.Start,'string');
%TODO need error checking
startTime = str2num(startTimesStr)
stopTime = str2num(stopTimesStr)

%change pathces to new values
set(hAnPatch.Full,'xdata',[stopTime stopTime startTime startTime])
set(hAnPatch.Slid,'xdata',[stopTime stopTime startTime startTime])

set(hAnPatch.Full,'visible','on')
set(hAnPatch.Slid,'visible','on')
% set tick box to manual
set(hAnWholeFile,'value',0)

end
%% TEMPLATE MATCH FUNCTIONS
function [] = launchTemplateMatch(src,eventData)
disp('launch tempalte match using existing spikes as templates')
hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
spikeTimes = getappdata(hSpikeGui,'spikeTimes');
channelData = getappdata(hSpikeGui,'channelData');
fs = channelData.fs; %TODO switch to processed data if below is ok with nans
data = channelData.filtSignal;


if isempty(spikeTimes)
    disp('currently only possible with spikes')
    
    return;
end

%get active units and then get spike shapes for each
activeUnits = find(~cellfun('isempty',spikeTimes));
spikeShapesCell = cell(1,max(activeUnits));
for unitNum = activeUnits
    
    unitSpikeTimes = spikeTimes{unitNum};
    [spikeShapesCell{unitNum} t] = extractSpikeShapes(data,unitSpikeTimes,[0.8 2],fs);
    
end

%TODO move this spike plot to analyssi button for now is debuggin here
hF = figure;
hA1 = subplot(1,2,1,'parent',hF); %for shapes
hA2 = subplot(1,2,2,'parent',hF); %for means and templates
set([hA1 hA2],'xlim',[t(1) t(end)])
for unitNum = activeUnits
    line(t,spikeShapesCell{unitNum},'color',cc(unitNum,:),'parent',hA1);
    line(t,mean(spikeShapesCell{unitNum},2),'color',cc(unitNum,:),'parent',hA2,...
        'linewidth',2);
end
%TODO get templates from spike shapes
%clusterStructure = generateTemplatesAO(epochsClustered, distances, centers)
%for now just use mean
templatesForMatching = cell(max(activeUnits),1); %to store the template shapes
templatesForMatchingTimeVec = cell(max(activeUnits),1); %to store aligned time vectors
templatePoints = [-5 0 5 10 15 20 25 30]; %indices relative to peak
indsToTest = [-templatePoints(1):length(mean(spikeShapesCell{unitNum},2))-templatePoints(end)-1]; %currently second poitn along at peak
for unitNum = activeUnits
    theseShapes = spikeShapesCell{unitNum};
    muShape = mean(theseShapes,2);
    
    %find max and minimum and which is bigger muShape
    [posPeak posPeakInd] = max(muShape(indsToTest));
    [negPeak negPeakInd]= min(muShape(indsToTest));
    
    
    if abs(negPeak)>posPeak
        %negative spike
        alignMentInd = negPeakInd+indsToTest(1)-1;
    else
        alignMentInd = posPeakInd+indsToTest(1)-1;
    end
    
    
    
    %get mean shape at those inds
    templatesForMatching{unitNum} = muShape(templatePoints+alignMentInd);
    templatesForMatchingTimeVec{unitNum} = t(templatePoints+alignMentInd);
    %draw as line on second axis
    line(templatesForMatchingTimeVec{unitNum},templatesForMatching{unitNum},....
        'color',cc(unitNum,:),'parent',hA2,'linestyle','none',...
        'marker','s','markersize',10,'linewidth',2)
    
end

%TODO run the template match and display the matched spikes and ssq
%histogram
%regenereate clusterStruc
numUnits = max(activeUnits);
for templ8Num = 1:numUnits
    
    cluster_struc(templ8Num).templ8 = templatesForMatching{templ8Num}';
end


%find which period of data to test
guiHands = getappdata(hSpikeGui,'guiHands');
hAnWholeFile = guiHands.hAnWholeFile; %handle to tick box controlling whole file analysis or sleected period

if get(hAnWholeFile,'value')
    startTimeS = 0;
    stopTimeS = length(processedData.data)/processedData.fs;
else
    hAnTimes = guiHands.hAnTimes; %handles to edit boxes with start and stop times for analysis
    startTimeS = str2num(get(hAnTimes.Start,'string')) %TODO error checking must be done at callback from edit boxes
    stopTimeS = str2num(get(hAnTimes.Stop,'string'))
    
end
settings.rate = fs;
stopSample = stopTimeS*settings.rate;
startSample =  startTimeS*settings.rate+1;

%get times that matter by trimming table of aretefact times to exclude to be within
%rate

releventTimesToExclude = channelData.timesToExclude;
if ~isempty(releventTimesToExclude)
    releventTimesToExclude = releventTimesToExclude(releventTimesToExclude(:,1)>startTimeS & releventTimesToExclude(:,2)<stopTimeS,:);
end




%basic settings
settings.ddthresh = 80;
settings.mergeClusters = 1;
settings.run_silent=1;
tic
[itm merged ddetect_all ssqFunction] = templateMatchAOdev(cluster_struc,numUnits,data(startSample:stopSample),settings);
toc



%TODO switch to using all itm.min not tightthmins
%TODO remove arefact times
% if ~isempty(timesToExclude)
%     %first convert table from s into samples, floor first colum, ceil
%     %second.
%     %loop over each perido and remove spikes
%     timesToExcludeSamples = timesToExclude*settings.rate;
%     numPeriods = size(timesToExclude,1);
%     for periodNum = 1:numPeriods
%
%        artSpikes = spikeIndexsAligned>timesToExcludeSamples(periodNum,1) & spikeIndexsAligned<timesToExcludeSamples(periodNum,2);
%
%        spikeIndexsAligned(artSpikes) = [];
%        spikesAligned(artSpikes,:) = [];
%
%     end
%
% end


%get times and spike shapes for all
templateMatchTimes = cell(1,numUnits);
templateMatchShapes = cell(1,numUnits);
for unitNum = 1:numUnits
    templateMatchTimes{unitNum} = startTimeS+sort(itm(unitNum).tightthmins(2,:)/fs); %TODO currently run whoe file
    
    
    [templateMatchShapes{unitNum} t] = extractSpikeShapes(data,templateMatchTimes{unitNum},[0.8 2],fs);
end

%create figure for displaying template matching results
hTempMatchResFig = figure('units','normalized',...
    'position',[0.1 0.1 0.8 0.8],...
    'toolbar','figure','tag','tTmResFig'); %,'CloseRequestFcn',{@closeMainWindow}
%1 axes for shapes, 1 for ssq histogram and a listbox for unit selection
hTmSsqHistAx = axes('parent',hTempMatchResFig,'units','normalized',...
    'position',[0.05 0.05 0.8 0.45])
hTmShapesAx = axes('parent',hTempMatchResFig,'units','normalized',...
    'position',[0.05 0.55 0.8 0.45])

%TODO swithc this plot to use changeTMUnit

%plot spike shapes for unit 1
unitNum =1;
hTempMatchesLine = line(t,templateMatchShapes{unitNum},'parent',hTmShapesAx,'color',cc(unitNum,:));
%also show template
%draw as line on second axis
templatesForMatchingTimeVec{unitNum} = templatesForMatchingTimeVec{unitNum}-templatesForMatchingTimeVec{unitNum}(1);
hTempLine = line(templatesForMatchingTimeVec{unitNum},templatesForMatching{unitNum},....
    'color','k','parent',hTmShapesAx,'linestyle','--',...
    'marker','s','markersize',10,'linewidth',2,'markerfacecolor','k')
%also show noise or background units %TODO find the level line crosses not
%classified
currentTmResults.templateMatchShapes = templateMatchShapes; %results structure to pass to otehr functions
currentTmResults.hTmShapesAx = hTmShapesAx; %spike shpes axis
currentTmResults.hTempLine = hTempLine; %template line handle
currentTmResults.hTempMatchesLine  = hTempMatchesLine;%handle to spike shapes
currentTmResults.hTmSsqHistAx = hTmSsqHistAx; %axis handle to ssq hist axis
currentTmResults.itm = itm;
currentTmResults.templatesForMatchingTimeVec = templatesForMatchingTimeVec;
currentTmResults.templatesForMatching = templatesForMatching;
setappdata(hSpikeGui,'currentTmResults',currentTmResults) %currently sotred in main gu appdata but could switch to separate one for this figure

%make sure tempalte line is ont top
uistack(hTempLine,'top')

%panel to hold drop down menu to selct which unit to view
hUnitSelPan = uipanel('parent',hTempMatchResFig,...
    'units','normalized','position',[0.85 0.7 0.1 0.2],...
    'title','Select Unit');
%list to select which unit to display and all option
a = 1:numUnits;
listStrings = cellstr(num2str(a(:)));
listStrings = vertcat(listStrings,'All');

hSelectedUnit = uicontrol('style','popupmenu',...
    'Parent',hUnitSelPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.9],...
    'string',listStrings,'value',1,... %{'1','2','3','4','+'}
    'tag','tUnitTarget','callback',{@changeTMUnit});


%plot histogram of ssq's
%should have a cutoff threshold set at the median
%this limits the maximum value of the ssq minima to bother plotting
%inhistogram
allMinima = itm(unitNum).mins(1,:);
matchMinima = itm(unitNum).tightthmins(1,:);
dispCutOff = itm(unitNum).noiseThreshold;


%panel to hold edit boxes ofor displaying and controlling thresolds
hSsqPan = uipanel('parent',hTempMatchResFig,...
    'units','normalized','position',[0.85 0.2 0.1 0.4],...
    'title','Hist Controls');
hDispCutOffEdit = uicontrol('parent',hSsqPan,'style','edit',...
    'units','normalized','position',[0.05 0.05 0.9 0.35],...
    'tooltipstring','display max','string',num2str(dispCutOff)); %TODO callback should be to change maxima and recalculate
hTightThreshEdit = uicontrol('parent',hSsqPan,'style','edit',...
    'units','normalized','position',[0.05 0.45 0.9 0.35],...
    'tooltipstring','threshold','string',num2str(itm(unitNum).detectionThreshold));



dispMinima = allMinima(allMinima<dispCutOff);

[minimaCount minimaBins] = hist(dispMinima,1000);
bar(minimaBins,minimaCount,'parent',hTmSsqHistAx);
set(hTmSsqHistAx,'nextplot','add')
tightMinimaCount = histc(matchMinima,minimaBins);

hB = bar(minimaBins,tightMinimaCount,'parent',hTmSsqHistAx,'facecolor','r');


uistack(hB,'top')




%TODO ask to accept results?
%copy section for adding to marker plot but with different marker style.

%TODO shift this updating of marker plot to separate function
%and then add to marker trace
hMarkerPlot = findobj(hSpikeGui,'tag','tMarkerPlot');
for unitNum = activeUnits
    %delte existing line and redraw
    %delete(findobj(hMarkerPlot,'tag',['tUnit ' num2str(unitNum)]))
    lineHeight = unitNum;
    unitSpikes = templateMatchTimes{unitNum};
    line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
        'linestyle','none','marker','d','MarkerFaceColor','k',...%cc(unitNum,:),...
        'tag',['tTmUnit ' num2str(unitNum)]);
    
end


end


function  [] = changeTMUnit(hDrop,eventData,hTmShapesAx,hTmSsqHistAx)
%callback from drop down menu to change teh displayed unit
hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;

currentTmResults = getappdata(hSpikeGui,'currentTmResults')
%get chosen value and check agains
selString = get(hDrop,'string');
selVal = get(hDrop,'value');
if strcmpi('all',selString{selVal})
    
else
    unitNum = str2num(selString{selVal});
    hTempMatchesLine = currentTmResults.hTempMatchesLine;
    theseShapes = currentTmResults.templateMatchShapes{unitNum};
    oldX = get(hTempMatchesLine,'xdata');
    
    %TODO this should be replaced with a redraw method
    delete(hTempMatchesLine)
    hTempMatchesLine = line(repmat(oldX{1}(:),1,size(theseShapes,2)),theseShapes,...
        'parent',currentTmResults.hTmShapesAx,'color',cc(unitNum,:))
    
    
    %also replace template
    hTempLine = currentTmResults.hTempLine; %template line handle
    templatesForMatchingTimeVec{unitNum} = currentTmResults.templatesForMatchingTimeVec{unitNum}-currentTmResults.templatesForMatchingTimeVec{unitNum}(1);
    set(hTempLine,'xdata',templatesForMatchingTimeVec{unitNum},...
        'ydata',currentTmResults.templatesForMatching{unitNum}','color','k');
    uistack(hTempLine,'top')
    currentTmResults.hTempMatchesLine =  hTempMatchesLine;
    currentTmResults.hTempLine =  hTempLine;
    
    
    
    hTmSsqHistAx = currentTmResults.hTmSsqHistAx;
    itm = currentTmResults.itm;
    %plot histogram of ssq's
    %should have a cutoff threshold set at the median
    %this limits the maximum value of the ssq minima to bother plotting
    %inhistogram
    
    delete(get(hTmSsqHistAx,'children')) %TODO replace with redraw method
    
    
    allMinima = itm(unitNum).mins(1,:);
    matchMinima = itm(unitNum).tightthmins(1,:);
    dispCutOff = itm(unitNum).noiseThreshold;
    
    dispMinima = allMinima(allMinima<dispCutOff);
    
    [minimaCount minimaBins] = hist(dispMinima,1000);
    bar(minimaBins,minimaCount,'parent',hTmSsqHistAx);
    set(hTmSsqHistAx,'nextplot','add')
    tightMinimaCount = histc(matchMinima,minimaBins);
    
    bar(minimaBins,tightMinimaCount,'parent',hTmSsqHistAx,'facecolor',cc(unitNum,:),'edgecolor',cc(unitNum,:));
    
end


setappdata(hSpikeGui,'currentTmResults',currentTmResults)
end
%% SPIKE SORTING FUNCTIONS

function [spikeShapes timeVector removedSpikeInds] = extractSpikeShapes(data,spikeTimes,eventWindow,fs)
%extract the spike shape for each of the given times
%


data = data(:); %ensure column

eventWindowSamp = [floor(eventWindow(1)*fs/1000) ceil(eventWindow(2)*fs/1000)]; %convert window to samples
eventTimesSamp = round(spikeTimes*fs); %convert triggerTimes to samples

%remove those too close to end or beginning
tooLate = eventTimesSamp>length(data)-eventWindowSamp(2)-1;
tooEarly = eventTimesSamp<eventWindowSamp(1)+1;
toRemove = sum([tooEarly tooLate],2);
removedSpikeInds = find(toRemove);
eventTimesSamp(logical(toRemove)) = [];



numEvents = length(eventTimesSamp);
eventLength = sum(eventWindowSamp)+1;
spikeShapes = zeros(eventLength,numEvents);
for eventNum = 1:numEvents
    spikeShapes(:,eventNum) = data(eventTimesSamp(eventNum)-eventWindowSamp(1):eventTimesSamp(eventNum)+eventWindowSamp(2));
end

timeVector = linspace(-eventWindowSamp(1)/fs,eventWindowSamp(2)/fs,eventLength);

end



function [] = launchSorter(src,eventData)
%function to launch the autmoatic sorting, should do just detection as well
addpath(genpath('C:\Users\Pete\CloudStation\cleanAO'))
disp('launch automatic sorting')
hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;

%TODO need noise or unsorted class
channelData = getappdata(hSpikeGui,'channelData');
processedData.fs = channelData.fs;
if isempty(channelData.filtSignal)
    processedData.data = channelData.rawSignal;
else
    
    processedData.data = channelData.filtSignal;
end
%TODO processedData = getappdata(hSpikeGui,'processedData');
%for now use non nan data to test


%TODO replace with getting from main gui (loading from file for now)
%settings = loadedSettings(selSet).settings;

load('C:\Users\Pete\CloudStation\cleanAO\settings files\sumanSettings.mat')
settingz = offlineHiddenSettings(settingz);
settingz.rate = processedData.fs;

%check if the manual threshold buttons are checked and if so the get teh
%threshold level(s)
guiHands = getappdata(hSpikeGui,'guiHands');
% hPosThreshTb; %threshold tickboxs
% guiHands.hNegThreshTb = hNegThreshTb;
% guiHands.hHcursorPos = hHcursorPos;
% guiHands.hHcursorNeg = hHcursorNeg;


%if either are not empty then switch to manual threshold level
if get(guiHands.hPosThreshTb,'value') && get(guiHands.hNegThreshTb,'value')
    %double threshold
    settingz.detection = 0;%manual mode
    posThresh = get(guiHands.hHcursorPos,'ydata');
    posThresh = posThresh(1)
    settingz.detectionSdMin = posThresh;
    settingz.threshDirection = 'both';
    settingz.detectionSdMax = 10*posThresh; %TODO this should be replaced with another line for artefacts
elseif get(guiHands.hPosThreshTb,'value') && get(guiHands.hNegThreshTb,'value')==0
    %double threshold
    settingz.detection = 0;%manual mode
    posThresh = get(guiHands.hHcursorPos,'ydata');
    posThresh = posThresh(1)
    settingz.detectionSdMin = posThresh;
    settingz.threshDirection = 'pos';
    settingz.detectionSdMax = 10*posThresh; %TODO this should be replaced with another line for artefacts
elseif get(guiHands.hPosThreshTb,'value')==0 && get(guiHands.hNegThreshTb,'value')
    %double threshold
    settingz.detection = 0;%manual mode
    posThresh = get(guiHands.hHcursorPos,'ydata');
    posThresh = posThresh(1)
    settingz.detectionSdMin = posThresh;
    settingz.detectionSdMax = 10*posThresh; %TODO this should be replaced with another line for artefacts
    settingz.threshDirection = 'neg';
end



% if ~isempty(negThresh) | ~isempty(posThresh)
% %currently both threshold must be equal so set to smallest %TODO fix
%     settings.detection = 0;%manual mode
%     settings.detectionSdMin = min(abs([negThresh]));
%     settings.threshDirection = 'both';
% end

%return; %TODO delte after testing

settingz.postPeakLengthMs = 1; %TODO write to settings file

settingz.prePeakLength = round(settingz.prePeakLengthMs*settingz.rate/1000);
%settings.postPeakLength = round(settings.epochLengthMs*settings.rate/1000);
settingz.postPeakLength = round(settingz.postPeakLengthMs*settingz.rate/1000);

settingz.epochLength = settingz.prePeakLength+settingz.postPeakLength+1;%round(settings.epochLengthMs*settings.rate/1000);
if settingz.epochLength<36
    disp('epoch too short for generating template')
    disp('length is 36 extending')
    settingz.epochLength = 37;
    settingz.postPeakLength = settingz.epochLength - settingz.prePeakLength ;
end

settingz.deadTime = round(settingz.deadTimeMs*settingz.rate/1000);
if settingz.deadTime<4
    settingz.deadTime = 4;
end

data = processedData.data(:)'; %ensure row vector %TODO fix in code
%TODO get size to analyse from patch. Limit to 1min for now for stability
%and time
% stopTimesS = 30;
% startTimesS = 0;

guiHands = getappdata(hSpikeGui,'guiHands');
%hAnPatch = guiHands.hAnPatch; %handle to patchs for analysis period
hAnWholeFile = guiHands.hAnWholeFile; %handle to tick box controlling whole file analysis or sleected period

if get(hAnWholeFile,'value')
    startTimeS = 0;
    stopTimeS = length(processedData.data)/processedData.fs;
else
    hAnTimes = guiHands.hAnTimes; %handles to edit boxes with start and stop times for analysis
    startTimeS = str2num(get(hAnTimes.Start,'string')) %TODO error checking must be done at callback from edit boxes
    stopTimeS = str2num(get(hAnTimes.Stop,'string'))
    
end

stopSample = stopTimeS*settingz.rate;
startSample =  startTimeS*settingz.rate+1;
%get times that matter by trimming table of times to exclude to be within
%rate

releventTimesToExclude = channelData.timesToExclude;
if ~isempty(releventTimesToExclude)
    releventTimesToExclude = releventTimesToExclude(releventTimesToExclude(:,1)>startTimeS & releventTimesToExclude(:,2)<stopTimeS,:);
end
[results] = sortSpikesAOdev(data(startSample:stopSample),settingz,[],releventTimesToExclude);

disp('finished sorting')
%turn results into spikeTimesCell
results.classification = results.clusteredSpikesFinal(:,1);
results.spikeIndexs = results.clusteredSpikesFinal(:,2);
results.spikeShapes = results.clusteredSpikesFinal(:,3:end);

results.spikeIndexs = results.spikeIndexs+startSample;

allAnalysisRuns = getappdata(hSpikeGui,'allAnalysisRuns') %struct(1,numRuns).releventTimesToExclude
%       .spikeTimes (1,numSpikes) (s)
%       .spikeClassification classification given (unmatched to previous or
%       anything?)
%       .stopTimeS
%       .startTimeS
currentNumRuns = size(allAnalysisRuns,2);
%results need tagging with time analysed
allAnalysisRuns(1,currentNumRuns+1).spikeTimes = results.spikeIndexs/settingz.rate;
allAnalysisRuns(1,currentNumRuns+1).spikeClassification = results.classification;
allAnalysisRuns(1,currentNumRuns+1).stopTimeS = stopTimeS;
allAnalysisRuns(1,currentNumRuns+1).startTimeS = startTimeS;
allAnalysisRuns(1,currentNumRuns+1).spikeShapes = results.spikeShapes;

setappdata(hSpikeGui,'allAnalysisRuns',allAnalysisRuns)

if currentNumRuns>=1
    %TODO ask if to add as is or to try to match.
    
    [res] = launchUnitMatchGui(allAnalysisRuns(1,currentNumRuns),allAnalysisRuns(1,currentNumRuns+1),settingz);
else
    
    activeUnits = unique(results.classification);
    spikeTimesCell = cell(1,max(activeUnits));
    for unitNum = activeUnits'
        spikeTimesCell{unitNum} = results.spikeIndexs(results.classification==unitNum)/settingz.rate;
        
        
    end
    
    %and then add to marker trace
    hMarkerPlot = findobj(hSpikeGui,'tag','tMarkerPlot');
    for unitNum = activeUnits'
        %delte existing line and redraw
        delete(findobj(hMarkerPlot,'tag',['tUnit ' num2str(unitNum)]))
        lineHeight = unitNum;
        unitSpikes = spikeTimesCell{unitNum};
        line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
            'linestyle','none','marker','s','MarkerFaceColor',cc(unitNum,:),...
            'tag',['tUnit ' num2str(unitNum)]);
        
    end
    
    numUnits = max(activeUnits); %will leave space for empty units
    set(hMarkerPlot,'ylim',[0 max(activeUnits)+1])
    
    
    %update list for adding units
    a = 1:numUnits+1;
    listStrings = cellstr(num2str(a(:)));
    listStrings = vertcat(listStrings,'+');
    set(findobj('tag','tUnitTarget'),'string',listStrings);
    
    selectedSpikeInds = cell(1,numUnits);
    setappdata(hSpikeGui,'spikeTimes',spikeTimesCell)
    setappdata(hSpikeGui,'selectedSpikeInds',selectedSpikeInds)
    sortingResults.results = results;
    sortingResults.settings = settingz;
    setappdata(hSpikeGui,'sortingResults',sortingResults); %TODO check how to lable for export in AO
end
end


function [res] = launchUnitMatchGui(existingSpikes,newSpikes,settings)
%called when trying to match newly sorted units to those that already exist
%will launch a figure which shows the new shapes and old shapes + common
%PCA in middle with different markers
hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cPalette.bg = [1 1 1]; %TODO remove after testinga s is set in initiliasaiton
cc = cPalette.spikeColours;


%TODO if doesn't exist then create
hUnitMatchFig = figure('menubar','none',...
    'name','Match Units',...
    'numbertitle','off',...
    'tag','tUnitMatchFig',...
    'units','normalized','position',[0.2 0.2 0.6 0.7],...
    'resize','on','color',cPalette.bg); %TODO 'closerequestfcn',{@detailFigClose}

setappdata(0,'hUnitMatchFig',hUnitMatchFig)

%need 3 panels

hOldSpikes.pan = uipanel('parent',hUnitMatchFig,...
    'units','normalized','pos',[0.005 0.005 0.325 0.8],...
    'title','Old Spikes');


hNewSpikes.pan = uipanel('parent',hUnitMatchFig,...
    'units','normalized','pos',[0.335 0.005 0.325 0.8],...
    'title','New Spikes');

hSpikeScat.pan = uipanel('parent',hUnitMatchFig,...
    'units','normalized','pos',[0.665 0.005 0.325 0.8],...
    'title','Spike Scatter');

%button to accept match
hAcceptPb = uicontrol('parent',hSpikeScat.pan,...
    'units','normalized','position',[0.05 0.05 0.9 0.2],...
    'style','push','string','Accept Labels','callback',{@acceptUnits});

%in old pikes panel split into 4 units (%TODO enable expansion of this later)
hOldSpikes.panelHandles = nan(1,4);
hOldSpikes.panelPositions = nan(4,4);
hOldSpikes.axHand = nan(1,4);
for panelNum = 1:4
    hOldSpikes.panelPositions(panelNum,:) = [0.05 1-0.225*(panelNum) 0.9 0.2];
    hOldSpikes.panelHandles(panelNum) = uipanel('parent',hOldSpikes.pan,...
        'units','normalized','pos',hOldSpikes.panelPositions(panelNum,:),...
        'title',['Unit Num ' num2str(panelNum)]);
    hOldSpikes.axHand(panelNum) = axes('parent',hOldSpikes.panelHandles(panelNum),...
        'units','normalized','position',[0.05 0.05 0.9 0.9]);
end


%do the same for new spikes panel but for all new units add a
%listbox
hNewSpikes.panelHandles = nan(1,4);
hNewSpikes.panelPositions = nan(4,4);
hNewSpikes.listBoxHandles = nan(1,4);
hNewSpikes.axHand = nan(1,4);
hNewSpikes.panelNum = nan(1,4);
for panelNum = 1:4
    hNewSpikes.panelPositions(panelNum,:) = [0.05 1-0.225*(panelNum) 0.9 0.2];
    hNewSpikes.panelHandles(panelNum) = uipanel('parent',hNewSpikes.pan,...
        'units','normalized','pos',hNewSpikes.panelPositions(panelNum,:),...
        'title',['Unit Num ' num2str(panelNum)]);
    
    hNewSpikes.listBoxHandles(panelNum) = uicontrol('parent',hNewSpikes.panelHandles(panelNum),...
        'units','normalized','pos',[0.8 0.3 0.15 0.3],...
        'style','popupmenu','value',panelNum,...
        'string',{'1';'2';'3';'4'},'callback',{@reassignUnit,newSpikes})
    hNewSpikes.axHand(panelNum) = axes('parent',hNewSpikes.panelHandles(panelNum),...
        'units','normalized','position',[0.05 0.05 0.7 0.9]);
    hNewSpikes.panelNum(panelNum) = panelNum; %to keep track of original
end




%find amplitude axis limits and set the same for all axis
graphMax = max(max(max(existingSpikes.spikeShapes)),max(max(newSpikes.spikeShapes)))

graphMin = min(min(min(existingSpikes.spikeShapes)),min(min(newSpikes.spikeShapes)))

set([hNewSpikes.axHand],'ylim',[graphMin graphMax])
set([hOldSpikes.axHand],'ylim',[graphMin graphMax])


oldCMap = cc;
newCMap = [cc+0.4];

newCMap(newCMap>1)=1;
%TODO rearrange colours here when units are switched around
%plot all the spike shapes split into units for both new and old units

hOldSpikes.lineHand = cell(1,4);

oldUnitNums = unique(existingSpikes.spikeClassification);
spikeT = linspace(-settings.prePeakLengthMs,settings.postPeakLengthMs,size(existingSpikes.spikeShapes,2));
for unitNum = oldUnitNums'
    
    theseSpikesLog = existingSpikes.spikeClassification==unitNum;
    theseSpikeShapes = existingSpikes.spikeShapes(theseSpikesLog,:)';
    
    
    hOldSpikes.lineHand{unitNum} = line(spikeT,theseSpikeShapes,'parent',hOldSpikes.axHand(unitNum),'color',oldCMap(unitNum,:)) %ccToUse(2*(unitNum-1)+1,:)
end


hNewSpikes.lineHand = cell(1,4);


newUnitNums = unique(newSpikes.spikeClassification);
for unitNum = newUnitNums'
    
    theseSpikesLog = existingSpikes.spikeClassification==unitNum;
    theseSpikeShapes = existingSpikes.spikeShapes(theseSpikesLog,:)';
    
    
    hNewSpikes.lineHand{unitNum} = line(spikeT,theseSpikeShapes,'parent',hNewSpikes.axHand(unitNum),'color',newCMap(unitNum,:))
end

%now calculate PCA of all spikes together
[COEFF,SCORE] = princomp([existingSpikes.spikeShapes;newSpikes.spikeShapes]);
% selComps = 1:settings.numFeatures;
% comps = SCORE(:,selComps);
comps = SCORE(:,1:2);

%now draw a scatter showing
hScat.Ax = axes('parent',hSpikeScat.pan,...
    'units','normalized','pos',[0.05 0.3 0.9 0.65]);

numOldSpikes = size(existingSpikes.spikeShapes,1);
hScat.oldSpikes = scatter3(comps(1:numOldSpikes,1),comps(1:numOldSpikes,2),existingSpikes.spikeTimes,20,oldCMap(existingSpikes.spikeClassification,:))
hold on;
hScat.newSpikes = scatter3(comps(numOldSpikes+1:end,1),comps(numOldSpikes+1:end,2),newSpikes.spikeTimes,20,newCMap(newSpikes.spikeClassification,:),'x')
axis tight %TODO some proper way of setting limts

uMFigHands.hScatAx = hScat; %axis with sactter plots on
uMFigHands.hNewSpikes = hNewSpikes;
uMFigHands.hOldSpikes = hOldSpikes;
setappdata(hUnitMatchFig,'uMFigHands',uMFigHands)
setappdata(hUnitMatchFig,'existingSpikes',existingSpikes)
setappdata(hUnitMatchFig,'newSpikes',newSpikes)


%need some way of comparing the similarity of the spike trains in the
%overlapping section
%round all spike times to ms precision
tPrecis = 0.001;
tPrecisMult = round(1/tPrecis);


%trim to the overlapping section
timeSlices = intersect([newSpikes.startTimeS:newSpikes.stopTimeS],[existingSpikes.startTimeS:existingSpikes.stopTimeS]);
overlapLims = [min(timeSlices) max(timeSlices)];

allPerShared = nan(4,4);
for newUnitNum = 1:4
    if ismember(newUnitNum,newUnitNums)
        theseNewSpikesLog = newSpikes.spikeClassification==newUnitNum;
        theseNewSpikeTimes = newSpikes.spikeTimes(theseNewSpikesLog);
        theseNewSpikeTimes = round((theseNewSpikeTimes*tPrecisMult))/tPrecisMult;
        
        theseNewSpikeTimes = theseNewSpikeTimes(theseNewSpikeTimes>overlapLims(1) & theseNewSpikeTimes<overlapLims(2));
        
        for oldUnitNum = 1:4
            
            if ismember(oldUnitNum,oldUnitNums)
                %compare all units against each other
                theseSpikesLog = existingSpikes.spikeClassification==oldUnitNum;
                theseSpikeTimes = existingSpikes.spikeTimes(theseSpikesLog);
                
                theseSpikeTimes = round((theseSpikeTimes*tPrecisMult))/tPrecisMult;
                %trim t only overlapping spikes
                theseSpikeTimes = theseSpikeTimes(theseSpikeTimes>overlapLims(1) & theseSpikeTimes<overlapLims(2));
                
                numShared = sum(ismember(theseSpikeTimes,theseNewSpikeTimes));
                perShared = 100*(numShared/length(theseSpikeTimes));
                allPerShared(newUnitNum,oldUnitNum) = perShared
                
            end
        end
        
    end
end
%plot percentage similarity as heat plot
figure;
hI = imagesc(allPerShared);
set(hI,'alphadata',~isnan(allPerShared));
%TODO should have way of determining most efficient arrangemnt possible

res = []; %should be matching unit numbers
end

function [] = acceptUnits(src,eventData)
%callback from push to accept teh unit match labels and add teh new spikes
%to the mainWindow spikeTimes cell.
%should check for duplicates and delete.
hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
hUnitMatchFig = getappdata(0,'hUnitMatchFig');
existingSpikes = getappdata(hUnitMatchFig,'existingSpikes');
newSpikes = getappdata(hUnitMatchFig,'newSpikes');
uMFigHands = getappdata(hUnitMatchFig,'uMFigHands');
%trim to the overlapping section
timeSlices = intersect([newSpikes.startTimeS:newSpikes.stopTimeS],[existingSpikes.startTimeS:existingSpikes.stopTimeS]);
overlapLims = [min(timeSlices) max(timeSlices)];

nonOverlapLims = existingSpikes.stopTimeS; %for no only new spikes not already assigned will be used
%TODO ask if to change overlap spikes to new or old classification


spikesToKeepLog = newSpikes.spikeTimes>existingSpikes.stopTimeS;

spikeTimesToKeep = newSpikes.spikeTimes(spikesToKeepLog);
spikeClassToKeep = newSpikes.spikeClassification(spikesToKeepLog);

numUnits = 4;
newSpikeTimesCell = cell(1,numUnits); %TODO this may give an error if merging with different sized existing cell
%get the numbers that should be assigned to each unit
finalClassification = uMFigHands.hNewSpikes.panelNum;

for unitNum = 1:numUnits
    
    newSpikeTimesCell{unitNum} = spikeTimesToKeep(spikeClassToKeep==finalClassification(unitNum));
    
end


%now merge the new spike times cell with the old
oldSpikeTimes = getappdata(hSpikeGui,'spikeTimes');
%this is a quick fix to pad the spikeTimesCell to be 4 units big and mathc
%above
if size(oldSpikeTimes,2)<4
    oldSpikeTimes = [oldSpikeTimes repmat({[]},4-size(oldSpikeTimes,2),1)];
end

spikeTimes = cellfun(@(x,y) [x(:) ;y(:)],oldSpikeTimes,newSpikeTimesCell,'uniformoutput',false);

activeUnits = find(~cellfun('isempty',spikeTimes));

%and then add to marker trace
hMarkerPlot = findobj(hSpikeGui,'tag','tMarkerPlot');
for unitNum = activeUnits
    %delte existing line and redraw
    delete(findobj(hMarkerPlot,'tag',['tUnit ' num2str(unitNum)]))
    lineHeight = unitNum;
    unitSpikes = spikeTimes{unitNum};
    line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
        'linestyle','none','marker','s','MarkerFaceColor',cc(unitNum,:),...
        'tag',['tUnit ' num2str(unitNum)]);
    
end

numUnits = max(activeUnits); %will leave space for empty units
set(hMarkerPlot,'ylim',[0 max(activeUnits)+1])




% %this was taken from the sortignfunction and needs to be done here when
% accepting units
% %update list for adding units
% a = 1:numUnits+1;
% listStrings = cellstr(num2str(a(:)));
% listStrings = vertcat(listStrings,'+');
% set(findobj('tag','tUnitTarget'),'string',listStrings);
%
% selectedSpikeInds = cell(1,numUnits);
% setappdata(hSpikeGui,'spikeTimes',spikeTimesCell)
% setappdata(hSpikeGui,'selectedSpikeInds',selectedSpikeInds)
% sortingResults.results = results;
% sortingResults.settings = settings;
% setappdata(hSpikeGui,'sortingResults',sortingResults);

%TODO close the figure
delete(hUnitMatchFig)
end

function [] = reassignUnit(hPop,eventData,newSpikes)
%callback from list boxes used for reassigning nwew units to old
hUnitMatchFig = getappdata(0,'hUnitMatchFig');
uMFigHands = getappdata(hUnitMatchFig,'uMFigHands')

%get colourmap for lines
hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
oldCMap = cc;
newCMap = [cc+0.4];
newCMap(newCMap>1)=1;



%compare hPop to the list of handles to determine which unit it was then
%figure out what to do
hNewSpikes = uMFigHands.hNewSpikes;



%find which unit slected
unitNum = find(hNewSpikes.listBoxHandles==hPop);
%find which slot it should go to
targetSlot = get(hPop,'value'); %what it is cahngin to
%find which unit is currently there
unitInPlace = hNewSpikes.panelNum(targetSlot);
%find unitNum's current place
currentSlot = find(hNewSpikes.panelNum==unitNum);


%now switch where the original units where in the list
hNewSpikes.panelNum(currentSlot) = unitInPlace%oldUnitNum

hNewSpikes.panelNum(targetSlot) = unitNum  %origUnitCurPlace;%newUnitNum

%now set positions
set(hNewSpikes.panelHandles(unitNum),'position',hNewSpikes.panelPositions(targetSlot,:),...
    'parent',hNewSpikes.pan)
set(hNewSpikes.panelHandles(unitInPlace),'position',hNewSpikes.panelPositions(currentSlot,:),...
    'parent',hNewSpikes.pan)

%now listbox values
set(hNewSpikes.listBoxHandles(unitNum),'value',targetSlot)
set(hNewSpikes.listBoxHandles(unitInPlace),'value',currentSlot)

%now set to correct colours
set(hNewSpikes.lineHand{unitNum},'color',newCMap(targetSlot,:))
set(hNewSpikes.lineHand{unitInPlace},'color',newCMap(currentSlot,:))

%now also update colours in scatter


%ifn both unit numbers spikes
unitNumSpikes = find(newSpikes.spikeClassification==unitNum);
unitInPlaceSpikes = find(newSpikes.spikeClassification==unitInPlace);

curScatColours = get(uMFigHands.hScatAx.newSpikes,'Cdata');

%set colours to new
curScatColours(unitNumSpikes,:)=repmat(newCMap(targetSlot,:),length(unitNumSpikes),1);

curScatColours(unitInPlaceSpikes,:)=repmat(newCMap(currentSlot,:),length(unitInPlaceSpikes),1);


set(uMFigHands.hScatAx.newSpikes,'Cdata',curScatColours);
drawnow;


uMFigHands.hNewSpikes = hNewSpikes;
setappdata(hUnitMatchFig,'uMFigHands',uMFigHands)
end
%% DISPLAY FUNCTIONS %TODO remove varargin
function [] = updateMaxSpikeInd(hUnitSel,eventData,hSpikeMaxSel)
unitSel = get(hUnitSel,'value')
hSpikeGui = getappdata(0,'hSpikeGui');
spikeTimes = getappdata(hSpikeGui,'spikeTimes');
bigSpikeInd = length(spikeTimes{unitSel})
set(hSpikeMaxSel,'string',num2str(bigSpikeInd))
end


function [] = centreOnSpike(src,eventData,hUnitSel,hSpikeIndSel,indModifier)
%function called when manually typing in which spike to go to
%3rd input is the handle to the listbox with the number of units
%indModifier - add this to the selected spike ind, is only called when
%buttons are pressed
%TODO need to updaet number of units in thsi listbox if more are added or
%subtracted
inString = get(hSpikeIndSel,'string');
%first check that it is an integer %TODO copy this section to the amp axis
%control manual function to accept decimals as well
if all(ismember(inString, '0123456789'))
    selSpikeInd = str2double(inString)+indModifier
    selSpikeUnit = get(hUnitSel,'value')
    
    %get spiketimes
    hSpikeGui = getappdata(0,'hSpikeGui');
    spikeTimes = getappdata(hSpikeGui,'spikeTimes');
    %check if between 1 and the number of spikes on this unit
    if selSpikeInd<=length(spikeTimes{selSpikeUnit}) && selSpikeInd>0
        
        selSpikeTime = spikeTimes{selSpikeUnit}(selSpikeInd);
        guiHands = getappdata(hSpikeGui,'guiHands');
        curWidth = diff(get(guiHands.hFullPlot,'xlim'));
        %copy to raster selecting
        newPlotLims = [selSpikeTime-curWidth/2 selSpikeTime+curWidth/2]
        updateMainPlots(newPlotLims);
        %if called by a button then update the spikeIndEditBox
        if indModifier~=0
            set(hSpikeIndSel,'string',num2str(selSpikeInd))
        end
    else
        
        disp('Please enter an integer between 1 and the maximum number of spikes in this unit')
    end
    
else
    disp('Please enter an integer between 1 and the maximum number of spikes in this unit')
end
end



function [] = turnOnOffLine(hTickBox,eventData,hTargetLine)
%callback from the tickboxes controlling the presence of the horizontal
%threshold lines
switch get(hTickBox,'value')
    case 0
        set(hTargetLine,'visible','off','hittest','off')
    case 1
        set(hTargetLine,'visible','on','hittest','on')
end
end

function [] = clickLine(hLineMoved,eventData,hSpikeGui,hOtherLine)
%function called when the threshold line is clicked before being dragged
%hOtherLine is the handle to the 2nd line which should be a moved to the
%opposite polarity but same value
set(hSpikeGui,'WindowButtonMotionFcn',{@dragLine,1,hLineMoved,hOtherLine},...
    'WindowButtonDownFcn','drawnow;',...
    'WindowButtonUpFcn',{@dragLine,2,hLineMoved,hOtherLine});
uiwait(hSpikeGui)
end

function [] = dragLine(hSpikeGui,eventData,moveType,hTargetLine,hOtherLine)
%function called when threshold line is dragged, should move line and
%update plot, currently only works for horizontal lines but code is in
%place for verticals

if moveType==1 %moving cursor
    %disp('down')
    %targetLine = varargin{4};
    targetAxis = get(hTargetLine,'Parent');
    curPoint = get(targetAxis,'CurrentPoint');
    %redraw line based on current point - need switch to work vertical and
    %horizontal
    %cursorLimits = get(targetLine,'UserData');
    %cursor limits need to load other cursor positions
    %     switch get(targetLine,'tag')
    %         case 'hor'
    %             if curPoint(1,2) > cursorLimits(1) && curPoint(1,2) < cursorLimits(2)
    set(hTargetLine,'YData',curPoint(:,2));
    %
    %if moving positive threshold then malso match by moving negtive one
    if curPoint(1,2)>0
        set(hOtherLine,'YData',-curPoint(:,2));
    else
        set(hOtherLine,'YData',-curPoint(:,2));
    end
    %else
    %                 disp('too big')
    %             end
    %         case 'ver'
    %             if curPoint(1,1) > cursorLimits(1) && curPoint(1,1) < cursorLimits(2)
    %                 set(targetLine,'XData',curPoint(:,1));
    %             else
    %                 disp('too big ver')
    %             end
    %
    %         otherwise
    %             disp('??')
    %     end
    drawnow;
elseif moveType==2 %release mouse
    %disp('up')
    %reset figure functions and resume
    set(hSpikeGui,'WindowButtonMotionFcn','',...
        'WindowButtonDownFcn','',...
        'WindowButtonUpFcn','');
    uiresume(hSpikeGui)
else
    disp('?')
end
end

function [] = launchFiltPlot(src,eventData)
%callback from toggle button to launch a separate plot at the bottom which
%displays a 30,6000 Hz fitlered version of the raw signal

hSpikeGui = getappdata(0,'hSpikeGui');
channelData = getappdata(hSpikeGui,'channelData');
filtFreq = [30 6000];
%now filter
[b a] = butter(4,filtFreq/(channelData.fs/2),'bandpass');
filtData = filtfilt(b,a,channelData.rawSignal);

mainFigPos = get(hSpikeGui,'position');
newFigPos = mainFigPos - [0 0 0 0.4];
hF = figure('units','normalized','position',newFigPos); %TODO should match the mian fig in width
hAx = axes('parent',hF,'units','normalized','position',[0.05 0.05 0.795 0.95]);

%get current plot lims
hFullPlot = findobj(hSpikeGui,'tag','tFullPlot');
curLims = get(hFullPlot,'xlim');
timeToDisplay = channelData.fullTimeVector>curLims(1) & channelData.fullTimeVector<curLims(2);


line(channelData.fullTimeVector(timeToDisplay),filtData(timeToDisplay),...
    'parent',hAx,'color','k','tag','tFiltFigLine')
figure(hF)


setappdata(hSpikeGui,'filtFigData',filtData)
%TODO this plot should have axes linked to the main plot and use a similar
%method of updating axis lims
hMarkerPlot = findobj(hSpikeGui,'tag','tMarkerPlot')
linkaxes([hFullPlot hMarkerPlot hAx],'x')
end


function [] = ampAxisMan(src,eventData,hFullPlot)
%callback from edit boxes to set ylims n main plot
hSpikeGui = getappdata(0,'hSpikeGui');
set(hFullPlot,'ylimmode','manual')
loLimStr = get(findobj(hSpikeGui,'tag','tYaxLimLow'),'string');
hiLimStr = get(findobj(hSpikeGui,'tag','tYaxLimHigh'),'string');

%TODO need to check if both proper numbers and that loLim<hiLim
loLim = str2num(loLimStr);
hiLim = str2num(hiLimStr);

if loLim >= hiLim
    errordlg('Low limit must be less than High limit')
    return;
end
set(hFullPlot,'ylim',[loLim hiLim])
%set check box to 0 (manual)
set(findobj(hSpikeGui,'tag','tAutoYScaleTb'),'value',0)
end


function [] = ampAxisControl(hTb,eventData,hFullPlot)
%callback from the tickbox that switches between manual and automatic yaxis
%scaling on main plot
%hTb - handle to tick box or edit box if called from there
hSpikeGui = getappdata(0,'hSpikeGui');

switch get(hTb,'value')
    case 1
        %switch to auto limit mode
        set(hFullPlot,'ylimmode','auto')
    case 0
        set(hFullPlot,'ylimmode','manual')
        loLimStr = get(findobj(hSpikeGui,'tag','tYaxLimLow'),'string');
        hiLimStr = get(findobj(hSpikeGui,'tag','tYaxLimHigh'),'string');
        
        %TODO need to check if both proper numbers and that loLim<hiLim
        loLim = str2num(loLimStr);
        hiLim = str2num(hiLimStr);
        
        if loLim >= hiLim
            errordlg('Low limit must be less than High limit')
            return;
        end
        set(hFullPlot,'ylim',[loLim hiLim])
        
end



end


function [] = clickTimeSlider(varargin)
%function for selecting timeSlider- should be able to drag to set new time
%points
%args: 1 - handle to slider
%3- handle to mainfig
lineMoved = varargin{1};
set(varargin{3},'WindowButtonMotionFcn',{@dragTimeSlider,1,lineMoved},...
    'WindowButtonDownFcn','drawnow;',...
    'WindowButtonUpFcn',{@dragTimeSlider,2,lineMoved});
uiwait(varargin{3})

end

function [] = dragTimeSlider(varargin)
%args: 1 - main fig handle
%2 - []
%3 - case
%4 - line handle
switch varargin{3}
    case 1
        disp('enable drag')
        
        
        targetLine = varargin{4};
        targetAxis = get(varargin{4},'Parent');
        curPoint = get(targetAxis,'CurrentPoint');
        
        %cursorLimits = get(targetLine,'UserData') ; %TODO set limits at
        %min/max of signal
        
        curVert = get(targetLine,'Vertices');
        curWidth = abs(curVert(1,1)-curVert(4,1)); %width of slider
        
        newXVert = [curPoint(1,1)-curWidth/2;curPoint(1,1)-curWidth/2;curPoint(1,1)+curWidth/2;curPoint(1,1)+curWidth/2];
        newVert = [newXVert curVert(:,2)];
        set(targetLine,'Vertices',newVert);
        
        newPlotLims = [newXVert(1) newXVert(3)]
        updateMainPlots(newPlotLims);
        drawnow;
    case 2
        disp('stop drag')
        set(varargin{1},'WindowButtonMotionFcn','',...
            'WindowButtonDownFcn','',...
            'WindowButtonUpFcn','');
        uiresume(varargin{1})
end


end

function [] = updateMainPlots(newLimits)
%function that is called to update the limits of the main signal plot
hSpikeGui = getappdata(0,'hSpikeGui');
hFullPlot = findobj(hSpikeGui,'tag','tFullPlot');
channelData = getappdata(hSpikeGui,'channelData');
processedData = getappdata(hSpikeGui,'processedData');
timeToDisplay = channelData.fullTimeVector>newLimits(1) & channelData.fullTimeVector<newLimits(2);

hSignal = findobj(hSpikeGui,'tag','tSignal');
set(hSignal,'xdata',channelData.fullTimeVector(timeToDisplay),'ydata',processedData.data(timeToDisplay));


%if it exists also update the dispalyed data in the filtered plot
hFiltFigLine = findobj('tag','tFiltFigLine');
if ~isempty(hFiltFigLine)
    filtData = getappdata(hSpikeGui,'filtFigData');
    set(hFiltFigLine,'xdata',channelData.fullTimeVector(timeToDisplay),'ydata',filtData(timeToDisplay));
    
    
end

%TODO this should be replaced with setting the xdata and ydata to a
%correctly downslampled(save memory) version of the signal
set(hFullPlot,'xlim',newLimits)

%also need to update time slider position
newSliderX = [newLimits(2); newLimits(2); newLimits(1); newLimits(1)];
hSlider = findobj('tag','tTimeSliderBar');
set(hSlider,'Xdata',newSliderX)
end

function [] = timeAxisControl(src,eventData,zoomOption)
%function called to step or zoom in time, can be called either by button
%presses or by arrow key keyboard shortcuts (via keyPress)
hSpikeGui = getappdata(0,'hSpikeGui');
hFullPlot = findobj('tag','tFullPlot');

curXlim = get(hFullPlot,'xlim');
switch zoomOption
    case 1  %zoom up big %when zooming keep center in same place
        disp('zoom in')
        curRange = diff(curXlim)
        curMiddle = mean(curXlim)
        newXlim = [curMiddle-curRange/4 curMiddle+curRange/4];
        %set(hFullPlot,'xlim',newXlim)
        updateMainPlots(newXlim);
    case 2  %zoom up small
        disp('zoom out')
        curRange = diff(curXlim)
        curMiddle = mean(curXlim)
        newXlim = [curMiddle-curRange curMiddle+curRange];
        updateMainPlots(newXlim);
    case 3  %zoom down big
        
        disp('shift time plus')
        stepSize = str2num(get(findobj('tag','tZoomEdit'),'string'))
        
        if ~isnumeric(stepSize)
            disp('please enter an appropriate number into the edit box')
            return;
        end
        
        newXlim = curXlim+stepSize;
        updateMainPlots(newXlim);
    case 4  %zoom down small
        disp('shift time neg')
        stepSize = str2num(get(findobj('tag','tZoomEdit'),'string'))
        
        if ~isnumeric(stepSize)
            disp('please enter an appropriate number into the edit box')
            return;
        end
        newXlim = curXlim-stepSize;
        updateMainPlots(newXlim);
end



end

function [] = keyPress(src,e)
%callback from keybaord press in main figure window
%used for keyboard shortcuts

disp('keypress detected')
switch e.Key
    case 'uparrow'
        %zoom in
        timeAxisControl([],[],1)
    case 'downarrow'
        %zoom out
        timeAxisControl([],[],2)
    case 'rightarrow'
        %disp('shift time pos')
        timeAxisControl([],[],3)
    case 'leftarrow'
        %disp('shift time neg')
        timeAxisControl([],[],4)
    case 'a'
        disp('you pressed a')
        %move back one spike if possible
    case 'd'
        disp('you pressed d')
        %move forward one spike if possible
    case 'x'
        disp('you pressed x')
    case '1'
        disp('Move to Unit 1')
        moveSpikes([],[],1)
    case '2'
        disp('Move to Unit 2')
        moveSpikes([],[],2)
    case '3'
        disp('Move to Unit 3')
        moveSpikes([],[],3)
    case '4'
        disp('Move to Unit 4')
        moveSpikes([],[],4)
    case '5'
        disp('Move to Unit 5')
        moveSpikes([],[],5)
    case 'backspace'
        disp('DELTING SPIKES!')
        deleteSpikes
    case 'n'
        markNewSpike
    otherwise
        disp('unrecognised keyboard command')
end
end
%% ANALYSIS FUNCTIONS
function [] = createSpikeSumFig(src,eventData)
disp('launch the spike summary figure')
%TODO in AO version this will point back to detail fig
hSpikeGui = getappdata(0,'hSpikeGui');
%spikeTimes = getappdata(hSpikeGui,'spikeTimes');
sortingResults = getappdata(hSpikeGui,'sortingResults');
channelData = getappdata(hSpikeGui,'channelData');
processedData.fs = channelData.fs;
processedData.data = channelData.filtSignal;
stopSample = floor(30*channelData.fs);
startSample = 1;
%TODO currently this runs of teh spike sorting reuslts not what has been
%edited
hPreviewFig = summariseSpikeFig(sortingResults.results,sortingResults.settings,processedData.data(startSample:stopSample),processedData.fs)

end

function [] = calculateFiringRates(src,eventData)
%function to calcualte the firing rate of each of the units
%needs to take the total non-nan amount of signal
hSpikeGui = getappdata(0,'hSpikeGui');
spikeTimes = getappdata(hSpikeGui,'spikeTimes');
if isempty(spikeTimes)
    errordlg('No Spikes');
    return;
end
processedData = getappdata(hSpikeGui,'processedData');
cPalette = getappdata(hSpikeGui,'cPalette');
fs = processedData.fs;
data = processedData.data;


cleanDataLength = sum(~isnan(data))/fs;
fullDataLength =  length(data)/fs;



numUnits = size(spikeTimes,2);
% spikeDensityFunction = cell(1,numUnits)
% for unitNum = 1:numUnits
%
%     spikeDensityFunction = ksdensity(spikeTimes{unitNum},spdT,'width',binSize);
%
%
% end
%TODO work out how to use ksdensity and correct for missing parts of data
%calculate as usual then nan the same parts as missing insignal and take
%nanmean??
blockSizeS = 30; %block size in s
blockSizeSamples = floor(blockSizeS*fs);
numBlocks = floor(fullDataLength/blockSizeS);

frResults = nan(numUnits,numBlocks);
for unitNum = 1:numUnits
    thisUnitsSpikes  = spikeTimes{unitNum};
    for blockNum = 1:numBlocks
        
        thisBlockStart = floor(blockSizeSamples*(blockNum-1)+1);
        
        thisBlockEnd = floor(blockSizeSamples*(blockNum));
        
        thisSegCleanLength = sum(~isnan(data(thisBlockStart:thisBlockEnd)))/fs;
        
        thisSegNumSpikes = sum(thisUnitsSpikes>blockSizeS*(blockNum-1) & thisUnitsSpikes<blockSizeS*(blockNum));
        frResults(unitNum,blockNum) = thisSegNumSpikes/thisSegCleanLength;
        
    end
end


hFfr = figure;
hAxfr = axes('parent',hFfr);
cc = cPalette.spikeColours;
for unitNum = 1:numUnits
    for blockNum = 1:numBlocks
        
        
        line([blockSizeS*(blockNum-1) blockSizeS*(blockNum)],[frResults(unitNum,blockNum) frResults(unitNum,blockNum)],...
            'parent',hAxfr,'color',cc(unitNum,:),'linewidth',2)
        
    end
    
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')


end

function [] = createIsiFig(src,eventData)
%callback from the ISI button
%should launch in separate figure a display showing the ISI of each unita s
%well as the rasters triggered of one another
hSpikeGui = getappdata(0,'hSpikeGui');
spikeTimes = getappdata(hSpikeGui,'spikeTimes')
if isempty(spikeTimes)
    errordlg('No Spikes');
    return;
end
hF = isiSubplotNd(spikeTimes)

end

function [] = launchInteractiveWaveView(src,eventData,updateFlag)
%callback from the push button to launch the interactive waveform picker
%if updateFlag is 0- then create fig from scratch if ot then plot into
%existing axis
hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
spikeTimes = getappdata(hSpikeGui,'spikeTimes');
numUnits = length(spikeTimes);
if updateFlag
    hIntWave = getappdata(0,'hIntWave');
    figure(hIntWave.hF); %make current figure
    delete(hIntWave.hL{:}); %delete existig lines
    %delete existing marked lines
    delete(findobj(hIntWave.hA,'tag','tClickedSpike'))
    visUnits = get(hIntWave.butPan.tickBox.tb,'value');
    visUnits = find(cellfun(@(x) x==1, visUnits));
else
    visUnits = 1:2;
%create fig and axis
hIntWave.hF = figure;
hIntWave.hA = axes('parent',hIntWave.hF,'units','normalized','position',[0.05 0.05 0.9 0.7]);
%need button panel with checkboxes for each unit
hIntWave.butPan.pan = uipanel('parent',hIntWave.hF,'units','normalized',...
    'position',[0.15 0.8 0.6 0.15],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Visible Units');
for unitNum = 1:3
hIntWave.butPan.tickBox.tb(unitNum) = uicontrol('style','checkbox','parent',hIntWave.butPan.pan,...
    'units','normalized',...
    'position',[0.15*(unitNum) 0.2 0.1 0.4],...
    'string',num2str(unitNum),'value',1,'Callback',{@changeVisUnits,unitNum});
end
set(hIntWave.butPan.tickBox.tb(3),'value',0,'enable','off'); %TODO disable until needed
%now want that updates the plot
hIntWave.butPan.pB = uicontrol('style','push','parent',hIntWave.butPan.pan,...
    'units','normalized',...
    'position',[0.7 0.2 0.25 0.4],...
    'string','Update','callback',{@updateIntWave});
end
%get data and loop over units to extract the spike shapes from the data
%TODO somehow should decide which units and what time window
processedData = getappdata(hSpikeGui,'processedData');

spikeShapes = cell(numUnits,1);
for unitNum = 1:numUnits
    
    
    if ~isempty(spikeTimes{unitNum})
        [spikeShapes{unitNum} timeVector] = extractSpikeShapes(processedData.data,spikeTimes{unitNum},[1 4],processedData.fs)
    end
end

timeVector = timeVector*1000;

%now plot each unit in the correct colour
hL = cell(numUnits,1);
for unitNum = 1:numUnits
    if ~isempty(spikeShapes{unitNum})
       hL{unitNum} = line(timeVector(:),spikeShapes{unitNum},'parent',hIntWave.hA,'color',cc(unitNum,:))
    end
end

%even if invisible at the moment the lines should be plotted
invisUnits = setdiff(1:2,visUnits);
if ~isempty(invisUnits)
for unitNum = invisUnits
    set(hL{unitNum},'visible','off')
end
end
hIntWave.hL =hL;
hIntWave.spikeShapes = spikeShapes;
hIntWave.timeVector = timeVector;
set(hIntWave.hA,'xlim',[timeVector(1) timeVector(end)],'ButtonDownFcn',{@clickWaveMark})
setappdata(0,'hIntWave',hIntWave)
% setappdata(hIntWave,'spikeShapes',spikeShapes) %moved to struct approach
% setappdata(hIntWave,'timeVector',timeVector)
end

function [] = updateIntWave(src,eventData)
%callback from update button on interactive wave marking. should delte all
%lines from axis and ten redraw using current spike times
% a= 1;
% hIntWave = getappdata(0,'hIntWave');
% delete(hIntWave.hF)
launchInteractiveWaveView([],[],1)
end

function [] = changeVisUnits(hTb,eventData,unitNum)
%function called when the checkboxes for visible units in the interactive
%waveform gui are pressed.
%it needs to turn them visible or invisible and also eliminate them from
%teh click search
% hTb = handle to calling tickbox
hIntWave = getappdata(0,'hIntWave');
visStat = get(hTb,'value');
if visStat
    set(hIntWave.hL{unitNum},'visible','on')
else
    set(hIntWave.hL{unitNum},'visible','off')
end
end

function [] =  clickWaveMark(hAwave,eventData)
%hAwave is axis of wavefrom plot
%callback from clicking on the waveform axis
%should find the closest waveform and highlight it
%also should only search through waeforms visible on screen
hIntWave = getappdata(0,'hIntWave');
spikeShapes = hIntWave.spikeShapes;
timeVector = hIntWave.timeVector;
hSpikeGui = getappdata(0,'hSpikeGui');
spikeTimes = getappdata(hSpikeGui,'spikeTimes');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
%now get user to click axis and get position
[xSel ySel]  = ginput(1);
%nowfind the closest sample
[~,closeSample] = min(abs(xSel-timeVector));
%spikeShapes-
visUnits = find(cellfun(@(x) x==1,get(hIntWave.butPan.tickBox.tb,'value')));
spikeMat = [spikeShapes{visUnits}];


[i spikeInd] =  min(abs(spikeMat(closeSample,:)-ySel))
%now we find which unit it was
cumSpikeCount = cumsum(cellfun(@(x) length(x),spikeTimes(visUnits)));
dSc = cumSpikeCount-spikeInd;
dSc(dSc<0) = nan;
[i unitNum] = nanmin(dSc);
if unitNum>1
    spikeInd = spikeInd - cumSpikeCount(unitNum-1);
end


%TODO should really delete old one
line(timeVector,spikeShapes{visUnits(unitNum)}(:,spikeInd),'parent',hAwave,'color',cc(visUnits(unitNum),:),'linewidth',2,'tag','tClickedSpike')

%now in the main gui we need to jump the time and centre on this spike
newPlotLims = [spikeTimes{visUnits(unitNum)}(spikeInd)-0.1 spikeTimes{visUnits(unitNum)}(spikeInd)+0.1]
updateMainPlots(newPlotLims);
drawnow;
end



function [] = launchInteractiveRaster(src,eventData)
%callback from itneractive raster button
%this shoudl create a figure in which there is one axis and a listbox
%the listbox will say which unit x v unit y to show and the axis will display a
%raster of unit x triggered off unit y which should be clickable. When
%clicking on the axis the nearest spike will be selected both in this gui
%and in the main gui
hSpikeGui = getappdata(0,'hSpikeGui');
spikeTimes = getappdata(hSpikeGui,'spikeTimes');

numUnits = length(spikeTimes);
%if only one unit then no point
if numUnits ==1
    errordlg('Only 1 unit!');
    return;
end
hIntRast = figure;
hA = axes('parent',hIntRast,'ButtonDownFcn',{@clickRastMark});

triggerUnit = 2;
respUnit = 1;
rastWinWidth = [0.2 0.2];
%for now by default we will just use unit 1 triggered off unit 2
unitSpikes = spikeTimes{respUnit};
triggerSpikes = spikeTimes{triggerUnit};
[hRast alignedRaster] = createRaster(unitSpikes,triggerSpikes,rastWinWidth,'axHandle',hA);

%
setappdata(0,'hIntRast',hIntRast)
setappdata(hIntRast,'triggerSpikes',triggerSpikes)
setappdata(hIntRast,'unitSpikes',unitSpikes)
end

function [] = clickRastMark(src,eventData)

hIntRast  = getappdata(0,'hIntRast');
triggerSpikes = getappdata(hIntRast,'triggerSpikes');
unitSpikes = getappdata(hIntRast,'unitSpikes');
%function called when clicking on the interactive raster axis
[xSel ySel]  = ginput(1);

%now need to find the actual time of the spike in question
unitNum = 1; %TODO for now we know that it will be asimple spike

%first find which line it is
triggerSpikeNum = round(ySel);
triggerSpikeTime = triggerSpikes(triggerSpikeNum);
%now find the closest spike to triggerSpikeTime+xSel

[i targetSpikeNum] = min(abs(unitSpikes-(triggerSpikeTime+xSel)));

%now need to mark it in the raster
relativeTargetSpikeTime = unitSpikes(targetSpikeNum)-triggerSpikeTime;

line([nan relativeTargetSpikeTime],[nan triggerSpikeNum],'linestyle','none','markerfacecolor','y',...
    'marker','s','parent',src)

%now in the main gui we need to jump the time and centre on this spike
newPlotLims = [unitSpikes(targetSpikeNum)-0.1 unitSpikes(targetSpikeNum)+0.1]
updateMainPlots(newPlotLims);
drawnow;
end

%% SPIKE CONTROL FUNCTIONS
function [] = markNewSpike(src,eventData)
%callback from mark new spike button
%should launch ginput and then mark the time selected in the signal plot as
%a spike, the unit number for the spike shoudl be taken from the same popup
%menu as when moving spikes


hSpikeGui = getappdata(0,'hSpikeGui');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
hMarkerPlot = findobj('tag','tMarkerPlot');
allCurrentSpikes = getappdata(hSpikeGui,'spikeTimes');

targetUnitNum = get(findobj('tag','tUnitTarget'),'value');
numUnits = size(allCurrentSpikes,2)
if targetUnitNum>numUnits
    %if user wants to introduce another unitNumber
    %create space for new spik times in new unit
    allCurrentSpikes = [allCurrentSpikes {[]}];
    
    %set listbox strings to include new unit
    a = 1:numUnits+1;
    listStrings = cellstr(num2str(a(:)));
    listStrings = vertcat(listStrings,'+');
    set(findobj('tag','tUnitTarget'),'string',listStrings);
    
    %increase y scale on marker plot
    hMarkerPlot = findobj('tag','tMarkerPlot');
    set(hMarkerPlot,'ylim',[0 numUnits+2])
    
end

[xSel ySel]  = ginput(1);

spikeTime = xSel

%get selected spikes and all spikes for this unit
%selectedSpikeInds = allSelectedSpikeInds{targetUnitNum};
currentSpikes = allCurrentSpikes{targetUnitNum};
currentSpikes = currentSpikes(:);

%add to target unit and sort by time
currentSpikes = sort([currentSpikes' spikeTime]);   %add to unit

allCurrentSpikes{targetUnitNum} = currentSpikes;

%delete marker line and redraw
delete(findobj('tag',['tUnit ' num2str(targetUnitNum)]))


lineHeight = targetUnitNum;
line(currentSpikes,lineHeight*ones(length(currentSpikes),1),'parent',hMarkerPlot,'color',cc(targetUnitNum,:),...
    'linestyle','none','marker','s','MarkerFaceColor',cc(targetUnitNum,:),...
    'tag',['tUnit ' num2str(targetUnitNum)]);
setappdata(hSpikeGui,'spikeTimes',allCurrentSpikes)
end


function [] = clSelSpikes(src,eventData)
%clears al selected spikes without changing anything
disp('clear selected spikes')
hSpikeGui = getappdata(0,'hSpikeGui');
allSelectedSpikeInds = getappdata(hSpikeGui,'selectedSpikeInds');
numUnits = size(allSelectedSpikeInds,2); %find total number of current assigned unit
selectedSpikeInds = cell(1,numUnits);
setappdata(hSpikeGui,'selectedSpikeInds',selectedSpikeInds);
%delte patches
hMarkerPlot = findobj('tag','tMarkerPlot');
%delete(findobj(hMarkerPlot,'type','patch'))
%replaced with line object for each unit
delete(findobj(hMarkerPlot,'-regexp','tag','tSelUnit '))
end

function [] = moveSpikes(varargin)
%move selected spikes to unit in listbox or keypress (will be third arg in
hSpikeGui = getappdata(0,'hSpikeGui');
allSelectedSpikeInds = getappdata(hSpikeGui,'selectedSpikeInds');
hMarkerPlot = findobj('tag','tMarkerPlot');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
numUnits = size(allSelectedSpikeInds,2); %find total number of current assigned units
if nargin==2
    targetUnitNum = get(findobj('tag','tUnitTarget'),'value'); %get drop box value for target of moving spikes
else
    targetUnitNum = varargin{3};
end
allCurrentSpikes = getappdata(hSpikeGui,'spikeTimes');
if targetUnitNum>numUnits
    
    %if user wants to introduce another unitNumber
    %create space for new spik times in new unit
    allCurrentSpikes = [allCurrentSpikes {[]}];
    
    %set listbox strings to include new unit
    a = 1:numUnits+1;
    listStrings = cellstr(num2str(a(:)));
    listStrings = vertcat(listStrings,'+');
    set(findobj('tag','tUnitTarget'),'string',listStrings);
    
    %increase y scale on marker plot
    hMarkerPlot = findobj('tag','tMarkerPlot');
    set(hMarkerPlot,'ylim',[0 numUnits+2])
    
    
    
end
%find selected spikes across all units
effectedUnits = find(~cellfun('isempty',allSelectedSpikeInds));
%don't do this to spikes already asigned to the same unit
effectedUnits = setdiff(effectedUnits,targetUnitNum);
%if ~isempty(effectedUnits)
for unitNum = effectedUnits;
    %get selected spikes and all spikes for this unit
    selectedSpikeInds = allSelectedSpikeInds{unitNum};
    currentSpikes = allCurrentSpikes{unitNum};
    currentSpikes = currentSpikes(:);
    if unitNum ~= targetUnitNum
        %add to target unit and sort by time
        allCurrentSpikes{targetUnitNum}  = sort([allCurrentSpikes{targetUnitNum}(:)' currentSpikes(selectedSpikeInds)']);   %add to unit 2
        %remove from unit
        currentSpikes(selectedSpikeInds) = [];
    end
    allCurrentSpikes{unitNum} = currentSpikes;
    
    %delete marker line and redraw
    delete(findobj('tag',['tUnit ' num2str(unitNum)]))
    
    
    lineHeight = unitNum;
    line(currentSpikes,lineHeight*ones(length(currentSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
        'linestyle','none','marker','s','MarkerFaceColor',cc(unitNum,:),...
        'tag',['tUnit ' num2str(unitNum)]);
    
    %blank selectedSpikeInds
    delete(findobj('tag',['tSelUnit ' num2str(unitNum)]));
    allSelectedSpikeInds{unitNum} = [];
end
% end
currentSpikes = allCurrentSpikes{targetUnitNum};
%redraw target unit line
lineHeight = targetUnitNum;
line(currentSpikes,lineHeight*ones(length(currentSpikes),1),'parent',hMarkerPlot,'color',cc(targetUnitNum,:),...
    'linestyle','none','marker','s','MarkerFaceColor',cc(targetUnitNum,:),...
    'tag',['tUnit ' num2str(targetUnitNum)]);
allSelectedSpikeInds{targetUnitNum} = [];
%delete selected spike lines
delete(findobj('tag',['tSelUnit ' num2str(targetUnitNum)]));
setappdata(hSpikeGui,'spikeTimes',allCurrentSpikes)
setappdata(hSpikeGui,'selectedSpikeInds',allSelectedSpikeInds)



%end
end

function [] = deleteSpikes(src,eventData)
%button to get all selected spike locations and delete
hSpikeGui = getappdata(0,'hSpikeGui');
allSelectedSpikeInds = getappdata(hSpikeGui,'selectedSpikeInds');
hMarkerPlot = findobj('tag','tMarkerPlot');

%find selected spikes across all units
effectedUnits = find(~cellfun('isempty',allSelectedSpikeInds));
allCurrentSpikes = getappdata(hSpikeGui,'spikeTimes');

cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
for unitNum = effectedUnits
    
    selectedSpikeInds = allSelectedSpikeInds{unitNum};
    
    currentSpikes = allCurrentSpikes{unitNum};
    
    currentSpikes(selectedSpikeInds) = [];
    unitSpikes = currentSpikes;
    %now overwits this cell and add to appdata
    allCurrentSpikes{unitNum} = currentSpikes;
    
    %delte existing line and redraw
    delete(findobj('tag',['tUnit ' num2str(unitNum)]))
    lineHeight = unitNum;
    line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
        'linestyle','none','marker','s','MarkerFaceColor',cc(unitNum,:),...
        'tag',['tUnit ' num2str(unitNum)]);
    
    %delete selected marker lines and blank selectedSpikeInds
    delete(findobj('tag',['tSelUnit ' num2str(unitNum)]));
    allSelectedSpikeInds{unitNum} = [];
end


setappdata(hSpikeGui,'spikeTimes',allCurrentSpikes)
setappdata(hSpikeGui,'selectedSpikeInds',allSelectedSpikeInds)
end



function [] = clickMarker(hAx,eventData)
%button to start ginput on marker axis, callback of tht axis btdfcn
hSpikeGui = getappdata(0,'hSpikeGui');
selectedSpikeInds = getappdata(hSpikeGui,'selectedSpikeInds');

[xSel ySel]  = ginput(1);
%if want to do it correctly:
%r=sqrt(((X-x(j))./XScale).^2+((Y-y(j))./YScale).^2);

%get closest y value representing a unit
%get tags and then split thte numbers from teh end of each tag
possibleUnitStrings = get(findobj('-regexp','tag','tUnit'),'tag');
tempStrings = cellfun(@(x) textscan(x, '%s %d'),possibleUnitStrings,'uniformoutput',false);
possibleUnitNums = double(cell2mat(cellfun(@(x) x(1,2),tempStrings)));

[dummy yCloseSel] = min(abs(ySel-possibleUnitNums));

unitNum = possibleUnitNums(yCloseSel)
% hL = findobj('tag',['tUnit ' num2str(unitNum)]);
% spikeTimes = get(hL,'xdata'); %this could isntead be taken from the spikeTimes variable
allCurrentSpikes = getappdata(hSpikeGui,'spikeTimes');
spikeTimes = allCurrentSpikes{unitNum}
[dummy xCloseSel] = min(abs(xSel-spikeTimes)); %find minimum distance

selectedSpikeInds{unitNum} = [selectedSpikeInds{unitNum} xCloseSel];

%selectedSpikeTime = spikeTimes(xCloseSel);
%draw black patch around selecteed point to indicate
% patch([selectedSpikeTime+0.001 selectedSpikeTime-0.001 selectedSpikeTime-0.001 selectedSpikeTime+0.001],...
%     [unitNum-0.5 unitNum-0.5 unitNum+0.5 unitNum+0.5],...
%     [0 0 0],'parent',hAx,'facealpha',0.5)
%repalced with a second line with bigger markers
selectedTimes = spikeTimes(selectedSpikeInds{unitNum});
hL = findobj('tag',['tSelUnit ' num2str(unitNum)]);
if isempty(hL)
    hL = line(selectedTimes,unitNum*ones(1,length(selectedTimes)),'color','k',...
        'marker','s','markerfacecolor','k','markersize',10,'parent',hAx,...
        'linestyle','none','tag',['tSelUnit ' num2str(unitNum)]);
else
    
    set(hL,'xdata',selectedTimes,'ydata',unitNum*ones(1,length(selectedTimes)))
end


setappdata(hSpikeGui,'selectedSpikeInds',selectedSpikeInds)

end
%% MISC FUNCTIONS
function [] = importSpikes(src,eventData)
%TODO should enable importing spike times either from workspace
%variables/from file/auto from AO converted template match

choice = questdlg('Please Choose Data Import Method', ...
    'Import Data','Load from disk','Load from other workspace','Cancel','Cancel');
if strcmp(choice,'Load from disk')
    %TODO must get specific varialbe in file
    [fileName,pathName,filterIndex] = uigetfile('.mat','Load Data',...
        ['-editedSpikeTimes.mat']);
    if filterIndex
        load([pathName fileName],'spikeTimes')
        hSpikeGui = getappdata(0,'hSpikeGui');
        setappdata(hSpikeGui,'spikeTimes',spikeTimes); %assing it to the gui
        
        %clear selcted spikes
        hL = findobj(hSpikeGui,'-regexp','tag',['tSelUnit \d']);
        if ~isempty(hL)
            delete(hL)
            
        end
        
        
        % replace lines for spikes
        cPalette = getappdata(hSpikeGui,'cPalette');
        numUnits = numel(spikeTimes);
        setappdata(hSpikeGui,'selectedSpikeInds',cell(numUnits,1))
        cc = cPalette.spikeColours;
        delete(findobj('-regexp','tag','tUnit '))
        guiHands = getappdata(hSpikeGui,'guiHands');
        hMarkerPlot = guiHands.hMarkerPlot;
        
        %plot line
        for unitNum = 1:numUnits
            lineHeight = unitNum;
            unitSpikes = spikeTimes{unitNum};
            line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
                'linestyle','none','marker','s','MarkerFaceColor',cPalette.spikeColours(unitNum,:),...
                'tag',['tUnit ' num2str(unitNum)]);
        end
        set(hMarkerPlot,'ylim',[0 numUnits+1]);
        set(hMarkerPlot,'yTick',[1:numUnits])
        
        
        %update unit list box in spike zoomer
        a = 1:numUnits;
        listStrings = cellstr(num2str(a(:)));
        hUnitSel = guiHands.hUnitSel;
        set(hUnitSel,'string',listStrings)
        %update max spike number
        updateMaxSpikeInd(hUnitSel,[],guiHands.hSpikeMaxSel)
        %update listbox for moving spikes to
        listStrings = vertcat(listStrings,'+');
        
        set(guiHands.hMoveSpikesTo,'string',listStrings)
    else
        errordlg('incorrect file name');
        return
    end
    
elseif strcmp(choice,'Load from other workspace')
    listVariablesInBaseBox
else
    %TODO error
    return;
end
end


function listVariablesInBaseBox
%produce figure with listbox of variables in base workspace and load or
%cancle buttons
figure('units','normalized','position',[0.4 0.4 0.2 0.2]);
lb = uicontrol('Style','listbox','units','normalized',...
    'Position',[0.05 0.25 0.9 0.7]);
%need select and cancel buttons
pb = uicontrol('Style','push','units','normalized',...
    'Position',[0.05 0.05 0.45 0.2],'string','Load Var',...
    'callback',{@selectVar,lb});
pb2 = uicontrol('Style','push','units','normalized',...
    'Position',[0.5 0.05 0.45 0.2],'string','Cancel',...
    'callback','close(gcf); return;');
%'Callback',@update_listBox

vars = evalin('base','who');
set(lb,'String',vars)


end


function [] = selectVar(src,eventData,hList)
%load selected varialbe from base workspace

%TODO there is no error checking at all on loading, should be spikeTimes
%style cell
hSpikeGui = getappdata(0,'hSpikeGui');
spikeTimes = getappdata(hSpikeGui,'spikeTimes');
cPalette = getappdata(hSpikeGui,'cPalette');
cc = cPalette.spikeColours;
variableList = get(hList,'string');
varSelection = get(hList,'value');

newSpikeTimes = evalin('base',variableList{varSelection});



%find existing units adn add as higher number to prevent overwirting
%ask what to assign it to
curActiveUnits = find(~cellfun('isempty',spikeTimes));

endOfOldUnits = max(curActiveUnits);
newUnitNums = find(~cellfun('isempty',newSpikeTimes));

newUnitNums = newUnitNums+endOfOldUnits;
newUnitsString = ['Load to units: ' num2str(newUnitNums)];
choice = questdlg(['Units currently active: ' num2str(curActiveUnits)], ...
    'Import Data',newUnitsString,'Overwrite Existing','Cancel','Cancel');

if strcmp(choice,newUnitsString)
    %disp('should add new spikes to correct units')
    
    spikeTimes = [spikeTimes newSpikeTimes];
    
    
elseif strcmp(choice,'Overwrite Existing')
    %disp('should overwirte existing units')
    spikeTimes = newSpikeTimes;
else
    %cancel
    return;
end


%set listbox strings to include new units
curActiveUnits = find(~cellfun('isempty',spikeTimes));
numUnits = max(curActiveUnits);
a = 1:numUnits+1;
listStrings = cellstr(num2str(a(:)));
listStrings = vertcat(listStrings,'+');
set(findobj('tag','tUnitTarget'),'string',listStrings);


hMarkerPlot = findobj('tag','tMarkerPlot');
curActiveUnits = curActiveUnits(:)';
%now also need to add to marker plot
for unitNum = curActiveUnits
    unitSpikes = spikeTimes{unitNum};
    
    
    
    %delte existing line and redraw (if exists)
    hL = findobj('tag',['tUnit ' num2str(unitNum)])
    if ~isempty(hL)
        delete(hL)
    end
    lineHeight = unitNum;
    line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',hMarkerPlot,'color',cc(unitNum,:),...
        'linestyle','none','marker','s','MarkerFaceColor',cc(unitNum,:),...
        'tag',['tUnit ' num2str(unitNum)]);
    
    %delete selected marker lines and blank selectedSpikeInds
    %delete(findobj('tag',['tSelUnit ' num2str(unitNum)]));
    
    
end
%increase y scale on marker plot

set(hMarkerPlot,'ylim',[0 numUnits+1],'ytick',[1:numUnits],'yticklabel',listStrings(1:end-1))

selectedSpikeInds = cell(1,numUnits);
setappdata(hSpikeGui,'spikeTimes',spikeTimes)
setappdata(hSpikeGui,'selectedSpikeInds',selectedSpikeInds)
end

function [] = saveAndExit(src,eventData)
%function called to save eitehr to disk or to ws
hSpikeGui = getappdata(0,'hSpikeGui');
spikeTimes = getappdata(hSpikeGui,'spikeTimes');
processedData = getappdata(hSpikeGui,'processedData');
channelData = getappdata(hSpikeGui,'channelData');
choice = questdlg('Please Choose Data Export Method', ...
    'Export Data','Save to disk','Assign to Caller','Cancel','Cancel');
if strcmp(choice,'Save to disk')
    disp('save data to disk')
    if ~isempty(channelData.multipleFileInfo)
        %need tod etect if mutliple files and then split spikes into relevant
        %ones
        
        numFiles  = size(channelData.multipleFileInfo,1);
        allProcessedData = processedData;
        allSpikeTimes = spikeTimes;
        clear spikeTimes processedData
        %need to save under these names so celr first
        for fileNum =  1:numFiles
            
            %processedDat and spikeTimes need splitting
            thisFileTimes = channelData.multipleFileInfo{fileNum,2};
            spikeTimes = cellfun(@(x) x(x>thisFileTimes(1) & x<thisFileTimes(2)),allSpikeTimes,'uniformoutput',false);
            %spiketimes for each file should be relative to start of that
            %file
            spikeTimes = cellfun(@(x) x-thisFileTimes(1),spikeTimes,'uniformoutput',false);
            
            processedData.fs =  allProcessedData.fs;
            startSample = floor(thisFileTimes(1)*processedData.fs)+1;
            endSample =  ceil(thisFileTimes(2)*processedData.fs)-1;
            processedData.data = allProcessedData.data(startSample:endSample);
            %fill in filename automatically
            [fileName,pathName,filterIndex] = uiputfile('.mat','Save Data',...
                [channelData.multipleFileInfo{fileNum,1}(1:end-4) '-editedSpikeTimes.mat']);
            if filterIndex
                save([pathName fileName],'spikeTimes','processedData') %should also save a copy of teh data sorted
            else
                errordlg('incorrect file name');
                return
            end
            clear spikeTimes processedData
        end
    else
        [fileName,pathName,filterIndex] = uiputfile('.mat','Save Data',...
            ['-editedSpikeTimes.mat']);
        if filterIndex
            save([pathName fileName],'spikeTimes','processedData') %should also save a copy of teh data sorted
        else
            errordlg('incorrect file name');
            return
        end
        
    end
elseif strcmp(choice,'Assign to Caller')
    disp('assign to caller')
    %assignin('caller','channelData',channelData)
    assignin('caller','editedSpikeTimes',spikeTimes)
    assignin('caller','processedData',processedData)
else
    return;
    
end

closeMainWindow %TODO auto exit?
end

function [] = closeMainWindow(src,eventData)
% userResponse = questdlg('Close?',...
%             'Close?','Yes','No','Yes');
% switch userResponse
%      case 'Yes'
%            delete(varargin{1})
%      case 'No'
%             return
%      otherwise
%             return
% end
uiresume;
delete(gcf)
end