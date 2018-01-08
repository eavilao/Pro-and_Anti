function [hArtGui] = launchPreProcessGui(data,fs)
%gui for preproccessing data before spike sorting. Should do a butterworth
%filter between given frequencies and also manually and automaitcally find
%artefacts then record the artefact periods. Shell is borrowed from
%spikeRemovalGui including slider.
% p.holland@erasmusmc.nl
%INPUTS:
%   data -  vector of data (data which is not to be considered should be
%   nan)
%   fs - sample rate in Hz
%OUTPUTS:
%    hArtGui - handle to main figure
%       data structure below is eitehr assigned out or saved
%    channelData.rawSignal - vecotr of raw signal
%                                       .fs - sampling frequency (Hz)
%                                       .filtSignal - filteredSignal
%                                       .lpFiltFreq - low pass freq in Hz
%                                       .hpFiltFreq - high pass freq in Hz
%                                       .timesToInclude - Currently unused so set to [-inf inf]                                       
%                                       .timesToExclude - table with times to exclude form further analysis (periodNum,[start end])  


%TODO
%filtered signal and artefacts shoudl be represented in time slider as well
%% TAGS %TODO remove findobj cals to tehse and repalce with an appdata handels structure
%'tArtSliderPlot' - axis for slider
%'tArtFullPlot' - axis for full plot
%'tArtTimeSliderBar' - patch object for time slider
%'tArtZoomEdit' - edit box for step size
%'tFiltSignal' - filtered signal 
%'tSignal'  -  raw signal
%'tArtefactTable' - talbe containing list of artefact times
%'tAutoYScaleTb' - tick box to control auto y axis scaling
%colours:
cPalette.subPanelbg = [0.8 0.8 0.8];
cPalette.spikeColours = lines(12); %TODO SNR colour scheme?
%settings
startDispTime = 5; %number of seconds to start displaying of signal

%% INPUTS
%TODO input parser and run silent



%% INITIALISE FIGURE
hArtGui = figure('units','normalized',...
    'position',[0.1 0.1 0.8 0.8],'CloseRequestFcn',{@closeMainWindow},...
    'toolbar','figure','tag','tArtGui','KeyPressFcn', {@keyPress});

setappdata(0,'hArtGui',hArtGui)

set(gcf,'renderer','zBuffer') %TODO why zbuffer?, transparency problems ?

%main axis for signal
hFullPlot = axes('Parent',hArtGui,'units','normalized',...
    'position',[0.05 0.2 0.795 0.65],'tag','tArtFullPlot');
%full signal to enable zoom control of main plot
hSliderPlot = axes('Parent',hArtGui,'units','normalized',...
    'position',[0.05 0.05 0.8 0.1],'tag','tArtSliderPlot');  


set(get(hFullPlot,'xLabel'),'string','Time (s)')
set(get(hFullPlot,'yLabel'),'string','Amplitude (\muV)')

%plot signal
fullTimeVector = linspace(0,length(data)/fs,length(data));
%now get logical to show which times to display
timeToDisplay = fullTimeVector>0 & fullTimeVector<startDispTime;
hSignal = line(fullTimeVector(timeToDisplay),...
    data(timeToDisplay),'color',[0 0 0],'tag','tSignal','Parent',hFullPlot,'hittest','off');
initialYLims = 1.5*[nanmin(data) nanmax(data)]; %Y axis limits for main plot
set(hFullPlot,'ylim',initialYLims);

%plot downsampled version of signal in the slider plot
dsFactor = 20;
dsData = data(1:dsFactor:end);
dsTimeVetor = linspace(0,length(data)/fs,length(dsData));

hDsSig = line(dsTimeVetor,...
    dsData,'color',[0 0 0],'Parent',hSliderPlot);
set(hSliderPlot,'xlim',[0 length(data)/fs])
set(hSliderPlot,'ylim',initialYLims)

%draw patch for slider main button between vertical cursors
hTimeSlider = patch([startDispTime startDispTime 0 0],2*[nanmin(data) nanmax(data) nanmax(data) nanmin(data)],...
             [0.8 0.8 0.8],'FaceAlpha',0.5,'tag','tArtTimeSliderBar','hittest','on',...
             'yliminclude','off');
         
set(hTimeSlider,'ButtonDownFcn',{@clickTimeSlider,hArtGui})

%horizontal cursor 
hHcursor1 = line([fullTimeVector(1) fullTimeVector(end)],[3*std(data) 3*std(data)],...
    'Color',[1 0 0],'LineWidth',2,'tag','tHorArtLine','Parent',hFullPlot);
set(hHcursor1,'ButtonDownFcn',{@clickLine,hArtGui,1})

%below is still needed for threshold lines etc
set(hFullPlot,'xlim',[0 startDispTime]) %TODO shoudl be replaced with setting xdata
%% BUTTONS

%2 edit boxes for hp and lp to button call butterworth
%TODO cahnge to button group
hFiltButPan = uipanel('parent',hArtGui,'units','normalized',...
    'position',[0.875 0.725 0.11 0.1],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Filter')
hHpFreqBox = uicontrol('style','edit',...
     'Parent',hFiltButPan,'units','normalized',...
    'position',[0.05 0.5 0.45 0.45],...
    'tag','tHpFreqBox','tooltipstring','High Pass Freq (Hz)');
hLpFreqBox = uicontrol('style','edit',...
     'Parent',hFiltButPan,'units','normalized',...
    'position',[0.5 0.5 0.45 0.45],...
    'tag','tLpFreqBox','tooltipstring','Low Pass Freq (Hz)');
% change to 4 buttons for HP, LP, BP, BS
hDoButter = uicontrol('style','push',...
    'Parent',hFiltButPan,'units','normalized',...
    'position',[0.05 0.05 0.45 0.45],...
    'string','Do Butter','Callback',{@doButter,hLpFreqBox,hHpFreqBox,fs});

%TODO add comb filter and PSD buttons

hArtListButPan = uipanel('parent',hArtGui,'units','normalized',...
    'position',[0.875 0.375 0.11 0.225],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Artefact Periods')
%table to store times of periods
hArtefactTable = uitable('Parent',hArtListButPan,...
    'units','normalized', 'Position',[0.05 0.15 0.9 0.8],...
    'ColumnName',{'Start','End'},...
    'ColumnFormat',{'numeric','numeric'},...
    'ColumnWidth','auto',...
    'RowName',[],...
    'CellSelectionCallback',{@selectArtefactPeriod},...
    'tag','tArtefactTable'); 
%force columns to take up half
%tabSize = getpixelposition(hArtefactTable);
%set(hArtefactTable,'ColumnWidth',{(floor(tabSize(3))/2)-1}); 
set(hArtefactTable,'ColumnWidth',{60})
%delete all marked periods
hDeleteAll = uicontrol('style','push',...
    'Parent',hArtListButPan,'units','normalized',...
    'position',[0.5 0.05 0.45 0.1],...
    'string','Delete All','Callback',{@deleteLast});
%delete marker periods
hDeleteOne = uicontrol('style','push',...
    'Parent',hArtListButPan,'units','normalized',...
    'position',[0.05 0.05 0.45 0.1],...
    'string','Delete Sel','Callback',{@deleteSelectedArtefact,hArtefactTable});


%threshold based detection
hThreshButPan = uipanel('parent',hArtGui,'units','normalized',...
    'position',[0.875 0.625 0.11 0.1],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Threshold Crossings')
%TODO tick box for positive or negative crossing artefact thresh or both or
%either
hPosNegTheshBg = uibuttongroup('parent',hThreshButPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.45],'tag','tPosNegTheshBg');
hPosThreshRb = uicontrol('parent',hPosNegTheshBg,'style','radio',...
    'units','normalized',...
    'position',[0.05 0.1 0.45 0.6],'string','pos')
hNegThreshRb = uicontrol('parent',hPosNegTheshBg,'style','radio',...
    'units','normalized',...
    'position',[0.55 0.1 0.45 0.6],'string','neg')
set(hPosNegTheshBg,'selectedobject',hPosThreshRb)
%push button with find artefacts above level line
hFindLineCrosses = uicontrol('style','push',...
    'Parent',hThreshButPan,'units','normalized',...
    'position',[0.05 0.5 0.9 0.45],...
    'string','Find Line Crosses','Callback',{@findLineCrosses,fs});


%manual artefact period marking
hManMarkButPan = uipanel('parent',hArtGui,'units','normalized',...
    'position',[0.875 0.175 0.11 0.05],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Manual Marking')
%mark begin
hMarkPeriod = uicontrol('style','push',...
    'Parent',hManMarkButPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.9],...
    'string','Mark Artefact Period','Callback',{@markBegin});


%automatic artefact perido finding based on rms of signal
%stablility finding button - will need expanding to include settings
hAutoArtButPan = uipanel('parent',hArtGui,'units','normalized',...
    'position',[0.875 0.275 0.11 0.1],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Auto Marking')
hFindArtefacts = uicontrol('style','push',...
    'Parent',hAutoArtButPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.25],...
    'string','Find Artefacts','Callback',{@artefactFinder});

%TODO edit boxes with defaults for stability settings
stabilitySettings.RmsWindowStepSize = 50; %size of step between consecutive windows (ms)
stabilitySettings.RmsWindowSize = 100; %window size in ms, is actually 2sd's value of gaussian window
stabilitySettings.RmsSmoothing = 4;  %minimum number of bad adjacent bad windows to count as bad sector
stabilitySettings.rmsSDfactor = 1.2; %number of sd's away from median RMS is threshold
%stabilitySettings.windowType = 1; %1 guassian, 2 box;
hStabWinWidthBox = uicontrol('style','edit',...
     'Parent',hAutoArtButPan,'units','normalized',...
    'position',[0.05 0.35 0.45 0.25],...
    'tag','tStabWinWidthBox','String',num2str(stabilitySettings.RmsWindowSize),...
    'TooltipString','window size (ms)');
hStabWinStepBox = uicontrol('style','edit',...
     'Parent',hAutoArtButPan,'units','normalized',...
    'position',[0.5 0.35 0.45 0.25],...
    'tag','tStabWinStepBox','String',num2str(stabilitySettings.RmsWindowStepSize),...
    'TooltipString','window step (ms)');
hStabConsecBox = uicontrol('style','edit',...
     'Parent',hAutoArtButPan,'units','normalized',...
    'position',[0.05 0.65 0.45 0.25],...
    'tag','tStabConsecBox','String',num2str(stabilitySettings.RmsSmoothing),...
    'TooltipString','number of consecutive bad windows');
hStabThreshBox = uicontrol('style','edit',...
     'Parent',hAutoArtButPan,'units','normalized',...
    'position',[0.5 0.65 0.45 0.25],...
    'tag','tStabThreshBox','String',num2str(stabilitySettings.rmsSDfactor),...
    'TooltipString','threshold multiplier');



% amplitude axis scale control
hYaxButPan = uipanel('parent',hArtGui,'units','normalized',...
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


%zoom buttons and time step buttons
hMainPlotPanel.zoomBg = uibuttongroup(hArtGui,'units','normalized',...
    'BackGroundColor',cPalette.subPanelbg,...
    'position',[0.875 0.075 0.11 0.1],'title','Time Axis Control');
hMainPlotPanel.zoomIn.pb = uicontrol(hMainPlotPanel.zoomBg,'style','push',...
    'units','normalized','position',[0.2 0.675 0.6 0.3],...
    'string','+','tag','artZoomUpBig','callback',{@timeAxisControl,1});
hMainPlotPanel.zoomOut.pb = uicontrol(hMainPlotPanel.zoomBg,'style','push',...
    'units','normalized','position',[0.2 0.025 0.6 0.3],...
    'string','-','tag','artZoomDownBig','callback',{@timeAxisControl,2});
hMainPlotPanel.stepForward.pb = uicontrol(hMainPlotPanel.zoomBg,'style','push',...
    'units','normalized','position',[0.675 0.35 0.3 0.3],...
    'string','->','tag','artZoomUpSmall','callback',{@timeAxisControl,3});
hMainPlotPanel.stepBackward.pb = uicontrol(hMainPlotPanel.zoomBg,'style','push',...
    'units','normalized','position',[0.025 0.35 0.3 0.3],...
    'string','<-','tag','artZoomDownSmall','callback',{@timeAxisControl,4});
hMainPlotPanel.zoomDownAmount.ed = uicontrol(hMainPlotPanel.zoomBg,'style','edit',...
    'units','normalized','position',[0.35 0.35 0.3 0.3],...
    'string','1','tag','tArtZoomEdit');


%save and exit
hImpExpButPan = uipanel('parent',hArtGui,'units','normalized',...
    'position',[0.875 0.025 0.11 0.05],...
    'BackGroundColor',cPalette.subPanelbg,...
    'title','Import/Export');
hSaveAndExit = uicontrol('style','push',...
    'Parent',hImpExpButPan,'units','normalized',...
    'position',[0.05 0.05 0.9 0.9],...
    'string','Save and Exit','Callback',{@saveAndExit});
%% FINALISE
%initiliase data structure   channelData.rawSignal - vecotr of raw signal
%                                       .fs - sampling frequency (Hz)
%                                       .filtSignal - filteredSignal
%                                       .lpFiltFreq - low pass freq in Hz
%                                       .hpFiltFreq - high pass freq in Hz
%                                       .timesToInclude - Currently unused so set to [-inf inf]                                       
%                                       .timesToExclude - table with times to exclude form further analysis (periodNum,[start end])  
channelData.fullTimeVector = fullTimeVector;
channelData.rawSignal = data;
channelData.fs = fs;
channelData.lpFiltFreq = [];
channelData.hpFiltFreq = [];
channelData.filtSignal = [];
channelData.timesToInclude = [-inf inf];
channelData.timesToExclude = [];
setappdata(hArtGui,'channelData',channelData)

hArtPatches = []; %this will contain the handles to the patches of each of the artefacts
setappdata(hArtGui,'hArtPatches',hArtPatches)
end



%% DISPLAY FUNCTIONS
%TODO display each signal tick box

function [] = ampAxisMan(src,eventData,hFullPlot)
%callback from edit boxes to set ylims n main plot
hArtGui = getappdata(0,'hArtGui');
set(hFullPlot,'ylimmode','manual')
loLimStr = get(findobj(hArtGui,'tag','tYaxLimLow'),'string');
hiLimStr = get(findobj(hArtGui,'tag','tYaxLimHigh'),'string');

%TODO need to check if both proper numbers and that loLim<hiLim
loLim = str2num(loLimStr);
hiLim = str2num(hiLimStr);

if loLim >= hiLim
    errordlg('Low limit must be less than High limit')
    return;
end
set(hFullPlot,'ylim',[loLim hiLim])
%set check box to 0 (manual)
set(findobj(hArtGui,'tag','tAutoYScaleTb'),'value',0)
end


function [] = ampAxisControl(hTb,eventData,hFullPlot)
%callback from the tickbox that switches between manual and automatic yaxis
%scaling on main plot
%hTb - handle to tick box or edit box if called from there
hArtGui = getappdata(0,'hArtGui');

switch get(hTb,'value')
    case 1
        %switch to auto limit mode
        set(hFullPlot,'ylimmode','auto')
    case 0
        set(hFullPlot,'ylimmode','manual')
        loLimStr = get(findobj(hArtGui,'tag','tYaxLimLow'),'string');
        hiLimStr = get(findobj(hArtGui,'tag','tYaxLimHigh'),'string');
        
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


function [] = clickLine(hLineMoved,eventData,hArtGui,unusedCaseArg)
%function called when the threshold line is clicked before being dragged

set(hArtGui,'WindowButtonMotionFcn',{@dragLine,1,hLineMoved},...
    'WindowButtonDownFcn','drawnow;',...
    'WindowButtonUpFcn',{@dragLine,2,hLineMoved});
uiwait(hArtGui)
end

function [] = dragLine(hArtGui,eventData,moveType,hTargetLine)
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
%             else
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
    set(hArtGui,'WindowButtonMotionFcn','',...
    'WindowButtonDownFcn','',...
    'WindowButtonUpFcn','');
    uiresume(hArtGui)
else
    disp('?')
end    
end

function [] = clickTimeSlider(hSlider,eventData,hArtGui)
%function for selecting timeSlider- should be able to drag to set new time
%points
%args: 1 - handle to slider
%3- handle to mainfig

set(hArtGui,'WindowButtonMotionFcn',{@dragTimeSlider,1,hSlider},...
            'WindowButtonDownFcn','drawnow;',...
            'WindowButtonUpFcn',{@dragTimeSlider,2,hSlider});
uiwait(hArtGui)

end

function [] = dragTimeSlider(hArtGui,eventData,moveType,hSlider)
%args: 1 - main fig handle
%2 - []
%3 - case
%4 - line handle
switch moveType
    case 1
        disp('enable drag')
        
        
        %targetLine = varargin{4};
        targetAxis = get(hSlider,'Parent');
        curPoint = get(targetAxis,'CurrentPoint');
        
        %cursorLimits = get(targetLine,'UserData') ; %TODO set limits at
        %min/max of signal
        
        curVert = get(hSlider,'Vertices');
        curWidth = abs(curVert(1,1)-curVert(4,1)); %width of slider
        
        newXVert = [curPoint(1,1)-curWidth/2;curPoint(1,1)-curWidth/2;curPoint(1,1)+curWidth/2;curPoint(1,1)+curWidth/2];
        newVert = [newXVert curVert(:,2)];
        set(hSlider,'Vertices',newVert);
        
        newPlotLims = [newXVert(1) newXVert(3)]
       updateMainPlots(newPlotLims);
       drawnow;
    case 2
        disp('stop drag')
        set(hArtGui,'WindowButtonMotionFcn','',...
                'WindowButtonDownFcn','',...
                'WindowButtonUpFcn','');
                uiresume(hArtGui)
end


end


function [] = updateMainPlots(newLimits)
%function that is called to update the limits of the main signal plot
hArtGui = getappdata(0,'hArtGui');
hFullPlot = findobj(hArtGui,'tag','tArtFullPlot');
channelData = getappdata(hArtGui,'channelData');
timeToDisplay = channelData.fullTimeVector>newLimits(1) & channelData.fullTimeVector<newLimits(2);

hSignal = findobj(hArtGui,'tag','tSignal');
set(hSignal,'xdata',channelData.fullTimeVector(timeToDisplay),'ydata',channelData.rawSignal(timeToDisplay));
%do the same for the filtered data if it exists
hFiltSignal = findobj(hArtGui,'tag','tFiltSignal');
if ~isempty(hFiltSignal)
    
   set(hFiltSignal,'xdata',channelData.fullTimeVector(timeToDisplay),'ydata',channelData.filtSignal(timeToDisplay));
 
end
%TODO this should be replaced with setting the xdata and ydata to a
%correctly downslampled(save memory) version of the signal
set(hFullPlot,'xlim',newLimits)

%also need to update time slider position
newSliderX = [newLimits(2); newLimits(2); newLimits(1); newLimits(1)];
hSlider = findobj('tag','tArtTimeSliderBar');
set(hSlider,'Xdata',newSliderX)
end

function [] = timeAxisControl(src,eventData,zoomOption)
%function called to step or zoom in time, can be called either by button
%presses or by arrow key keyboard shortcuts (via keyPress)
%hArtGui = getappdata(0,'hArtGui');
hFullPlot = findobj('tag','tArtFullPlot');

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
        stepSize = str2num(get(findobj('tag','tArtZoomEdit'),'string'))
        
        if ~isnumeric(stepSize)
            disp('please enter an appropriate number into the edit box')
            return;
        end
       
        newXlim = curXlim+stepSize;
        updateMainPlots(newXlim);
    case 4  %zoom down small
        disp('shift time neg')
        stepSize = str2num(get(findobj('tag','tArtZoomEdit'),'string'))
        
        if ~isnumeric(stepSize)
            disp('please enter an appropriate number into the edit box')
            return;
        end
        newXlim = curXlim-stepSize;
        updateMainPlots(newXlim);
end



end

function [] = keyPress(src,eventData)
%callback from keybaord press in main figure window
%used for keyboard shortcuts

disp('keypress detected')
switch eventData.Key
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
    
    otherwise
        disp('unrecognised keyboard command')
end
end
%% PROCESSING FUNCTIONS
function [] = doButter(src,eventData,hLpFreqBox,hHpFreqBox,fs)
%perform a butterworth filter using teh values displayed in the edit boxes
%need to check if lp<hP
%width [0.9 1.1]
%TODO get fs from appdata
hArtGui = getappdata(0,'hArtGui');
channelData = getappdata(hArtGui,'channelData');
data = channelData.rawSignal;

lpFreq = get(hLpFreqBox,'string');
hpFreq = get(hHpFreqBox,'string');
switch num2str(cellfun('isempty',isstrprop({lpFreq; hpFreq},'digit'))')
    case '1  0' %do high pass
        disp('do high pass')
        %hpFreq
        filtFreq = str2num(hpFreq);
        filtType = 'high';
    case '0  1' %do low pass
        disp('do low pass')
        %[b a] = butter(4,lpFreq(fs/2),'low');
        filtType = 'low';
        filtFreq = str2num(lpFreq);
    case'0  0' %do band pass
        disp('do band pass')
        filtFreq = [str2num(hpFreq) str2num(lpFreq)];
       
        
        filtType = 'bandpass';
        %[b a] = butter(4,[lpFreq*0.9/(fs/2) hpFreq*1.1/(fs/2)],'stop');
    otherwise
        
        disp('choose correct filter')
        return;
        %isempty(str2num(lpFreq))
        %isempty(str2double(hpFreq))
        
end

%use buttord to determine poles
[b a] = butter(4,filtFreq/(fs/2),filtType);
filtData = filtfilt(b,a,data);

%plot
hFiltFullPlot = findobj('tag','tArtFullPlot');
%delete existing plot
delete(findobj(get(hFiltFullPlot,'Children'),'tag','tFiltSignal'));
%get current limits and only plot the relevant data
newLimits = get(hFiltFullPlot,'xlim');

timeToDisplay = channelData.fullTimeVector>newLimits(1) & channelData.fullTimeVector<newLimits(2);


hFiltSignal = line(channelData.fullTimeVector(timeToDisplay),...
    filtData(timeToDisplay),'color',[0 0 1],'tag','tFiltSignal','Parent',hFiltFullPlot,...
    'hittest','off');

channelData.lpFiltFreq = str2num(lpFreq);
channelData.hpFiltFreq = str2num(hpFreq);
channelData.filtSignal = filtData;
setappdata(hArtGui,'channelData',channelData);
end
%% ARTEFACT FUNCTIONS

function [] = deleteSelectedArtefact(src,eventData,hTable)
%should delete the selected artefact in the table (if any selected)

hArtGui = getappdata(0,'hArtGui');
channelData = getappdata(hArtGui,'channelData')
hArtPatches = getappdata(hArtGui,'hArtPatches')



eventData = get(hTable,'UserData');
artefactNum = eventData.Indices(1);

%
channelData.timesToExclude(artefactNum,:) = [];

%find patch to delete %TODO won't work if not in the order they were added.
%need separate handles structure with each artefact patch handle in.
%delete(findobj('tag',['tArtefactPeriod ' num2str(artefactNum)]))
delete(hArtPatches(artefactNum));
hArtPatches(artefactNum) = [];
%replace table data with new list and remove the highlighted selection
set(hTable,'data',channelData.timesToExclude)

setappdata(hArtGui,'channelData',channelData)
setappdata(hArtGui,'hArtPatches',hArtPatches)
end

function [] = selectArtefactPeriod(hTable,eventData)
%callback fromartefact table when a cell is selected
%TODO jump main plot to this artefact
hArtGui = getappdata(0,'hArtGui');
% channelData = getappdata(hArtGui,'channelData');

if isempty(eventData.Indices)
    disp('No artefact periods')
    return;
end


set(hTable,'UserData',eventData) %this is trick to make sure that the selection is available to other functions

curData = get(hTable,'data');
artefactNum = eventData.Indices(1);
pointInArt = eventData.Indices(2); %beginning or end of artefact
thisArtTime = curData(artefactNum,:);

hFullPlot = findobj(hArtGui,'tag','tArtFullPlot');
currentVisibleWidth = diff(get(hFullPlot,'xlim'));

%set to begining or end of artefact (depends on which column of table clicked) period with ucrrent width (time resolution)
newLimits = [thisArtTime(pointInArt)-currentVisibleWidth/2 thisArtTime(pointInArt)+currentVisibleWidth/2];


updateMainPlots(newLimits)
%TODO enable deletion of single periods
end

function [] = artefactFinder(src,eventData)
%launches automatic artefact detecction using settings in boxes and then
%adds to table and patches

patchLim = [-1000 1000]; %TODO make sure always enough to stay away from ylims
hArtGui = getappdata(0,'hArtGui');
channelData = getappdata(hArtGui,'channelData');
fs = channelData.fs; 

if isempty(channelData.filtSignal)
    data = channelData.rawSignal;
else
    data = channelData.filtSignal;
end


%TODO remove previous found artefacts

stabilitySettings.RmsWindowStepSize = str2num(get(findobj(hArtGui,'tag','tStabWinStepBox'),'String')); %size of step between consecutive windows (ms)
stabilitySettings.RmsWindowSize = str2num(get(findobj(hArtGui,'tag','tStabWinWidthBox'),'String')); %window size in ms, is actually 2sd's value of gaussian window
stabilitySettings.RmsSmoothing = str2num(get(findobj(hArtGui,'tag','tStabConsecBox'),'String'));  %minimum number of bad adjacent bad windows to count as bad sector
stabilitySettings.rmsSDfactor = str2num(get(findobj(hArtGui,'tag','tStabThreshBox'),'String')); %number of sd's away from median RMS is threshold
stabilitySettings.windowType = 1; %1 guassian, 2 box;
%stabilitySettings.samplingRate = fs;

% windowSizeSamples = floor(stabilitySettings.samplingRate/(1000/stabilitySettings.RmsWindowSize));  %calculate window size in samples
% windowStepSizeSamples = floor(stabilitySettings.samplingRate/(1000/stabilitySettings.RmsWindowStepSize));  %calculate window step size in samples
%               
% gausWinSize = length([-ceil(windowSizeSamples/1.5):ceil(windowSizeSamples/1.5)]);

 
artefactPeriods = findArtefacts(data,fs,stabilitySettings);


%get current artefact times and patches  
hArtPatches = getappdata(hArtGui,'hArtPatches');
currentArtefactPeriods = channelData.timesToExclude;


%now mark all with patch
numArtefacts = size(artefactPeriods,1);
hFullPlot = findobj(hArtGui,'tag','tArtFullPlot');
hNewArtPatches = nan(1,numArtefacts);
for artefactNum = 1:numArtefacts
       hNewArtPatches(artefactNum) = patch([artefactPeriods(artefactNum,1) artefactPeriods(artefactNum,1) artefactPeriods(artefactNum,2) artefactPeriods(artefactNum,2)],...
        [patchLim(1) patchLim(2) patchLim(2) patchLim(1)],...
        [0 1 0],'FaceAlpha',0.5,'Parent',hFullPlot,'EdgeColor',[0 1 0],'EdgeAlpha',0.5,...
        'tag',['tArtefactPeriod ' num2str(size(currentArtefactPeriods,1)+artefactNum)],...
        'yliminclude','off');
end

%sort by start time order
[currentArtefactPeriods sInds] = sortrows([currentArtefactPeriods ;artefactPeriods]);
%combine new and old patches and arrange in same order as in table
hArtPatches = [hArtPatches hNewArtPatches];
hArtPatches = hArtPatches(sInds);

%add to table
set(findobj('tag','tArtefactTable'),'Data',currentArtefactPeriods)
channelData.timesToExclude = currentArtefactPeriods;
setappdata(hArtGui,'channelData',channelData)
setappdata(hArtGui,'hArtPatches',hArtPatches); 
end

function [] = markBegin(src,eventData)
%mark artefact period using two impoints, add to table and draw patch
hFullPlot = findobj('tag','tArtFullPlot');
patchLim = [-1000 1000]; %TODO make sure always enough to stay away from ylims
%create two impoints to select artefact periods and get positions
% hFirstSel = impoint(hFullPlot);
% firstSelPos = getPosition(hFirstSel);
% hSecondSel = impoint(hFullPlot);
% secondSelPos = getPosition(hSecondSel);

[clickedXcoords t2] = ginput(2);


%get existing artefact periods and patch handles
hArtGui = getappdata(0,'hArtGui');
hArtPatches = getappdata(hArtGui,'hArtPatches');
channelData = getappdata(hArtGui,'channelData');
%add to current data and sort by start time
currentPeriods = channelData.timesToExclude;
%[currentPeriods sInds] = sortrows([currentPeriods; [firstSelPos(1) secondSelPos(1)]]);
[currentPeriods sInds] = sortrows([currentPeriods; clickedXcoords']);
%add to table
set(findobj('tag','tArtefactTable'),'Data',currentPeriods);

%draw patch around points
% hP = patch([firstSelPos(1) firstSelPos(1) secondSelPos(1) secondSelPos(1)],[patchLim(1) patchLim(2) patchLim(2) patchLim(1)],...
%     [1 0 0],'FaceAlpha',0.5,'Parent',hFullPlot,...
%     'tag',['tArtefactPeriod ' num2str(size(currentPeriods,1))],'yliminclude','off');
hP = patch([clickedXcoords(1) clickedXcoords(1) clickedXcoords(2) clickedXcoords(2)],[patchLim(1) patchLim(2) patchLim(2) patchLim(1)],...
    [1 0 0],'FaceAlpha',0.5,'Parent',hFullPlot,...
    'tag',['tArtefactPeriod ' num2str(size(currentPeriods,1))],'yliminclude','off');
%add to handles list and put into order
hArtPatches = [hArtPatches hP];
hArtPatches = hArtPatches(sInds);
setappdata(hArtGui,'hArtPatches',hArtPatches);

%clear points     
%delete([hFirstSel; hSecondSel])    
%add back to the data
channelData.timesToExclude = currentPeriods;
setappdata(hArtGui,'channelData',channelData);
end

function [] = deleteLast(src,eventData)
% %TODO should be able to delete those
%selected in the table for no just clears all
%get full plot axis
hArtGui = getappdata(0,'hArtGui');
channelData = getappdata(hArtGui,'channelData');
currentPeriods = channelData.timesToExclude;
hArtPatches = getappdata(hArtGui,'hArtPatches'); 

switch size(currentPeriods,1)
    case 0 
        disp('none to delete')
        errordlg('no marked period to delete','NO MARKED PERIODS')
        return
   otherwise
         delete(hArtPatches); %delete patch   
         %delete(findobj('tag',[currentPeriods 'Filt'])); %fitlered patch
    %more than 1 period
        %currentPeriods = sortrows(currentPeriods,-1);%descending nnumber order
        %delete(findobj('tag',currentPeriods{1})); %delete patch 
        %delete(findobj('tag',[currentPeriods{1} 'Filt'])); %fitlered patch
end
%delete last row of table
%currentTablePeriods = get(findobj('tag','tArtefactTable'),'Data');
%currentTablePeriods(end,:) = [];

channelData.timesToExclude = [];
set(findobj('tag','tArtefactTable'),'Data',[]);
setappdata(hArtGui,'channelData',channelData);
hArtPatches = [];
setappdata(hArtGui,'hArtPatches',hArtPatches); 
end

function [] = findLineCrosses(src,eventData,fs)
%function that detects crossings of teh artefact line and then finds
%consecutive sections and assigns them as artefacts
hArtGui = getappdata(0,'hArtGui');
channelData = getappdata(hArtGui,'channelData');


patchLim = [-1000 1000]; %TODO make sure always enough to stay away from ylims

artLevel = get(findobj('tag','tHorArtLine'),'YData');
artLevel = artLevel(1);
%use fitlered data if available, must work with nans
%hFiltSignal = findobj('tag','tFiltSignal');


if isempty(channelData.filtSignal)
    data = channelData.rawSignal;
else
    data = channelData.filtSignal;
end

% is threshold positive or negative. %TODO should also detect using both or
% either as options
switch get(get(findobj('tag','tPosNegTheshBg'),'selectedobject'),'string')
    case 'pos'    
        badTimes = data>artLevel;
    case 'neg'
        badTimes = data<artLevel;  
end

%find consecutive sections
[c,valnew] = accumconncomps(badTimes); %TODO replace
if c(1)==1
    c = [1 c];
    valnew = [1 valnew];
end    

if c(end)==0
    c = c(1:end-1);
    valnew = valnew(1:end-1);
end 
%get table of times (s) with start in first column end in second column
artefactPeriods = reshape(cumsum(valnew),2,[])';
artefactPeriods = artefactPeriods/fs;
currentArtefactPeriods = channelData.timesToExclude;
%now mark all with patch
numArtefacts = length(artefactPeriods);
hFullPlot = findobj('tag','tArtFullPlot');
hNewArtPatches = nan(1,numArtefacts);
for artefactNum = 1:numArtefacts
       hNewArtPatches(artefactNum) = patch([artefactPeriods(artefactNum,1) artefactPeriods(artefactNum,1) artefactPeriods(artefactNum,2) artefactPeriods(artefactNum,2)],...
        [patchLim(1) patchLim(2) patchLim(2) patchLim(1)],...
        [0 1 0],'FaceAlpha',0.5,'Parent',hFullPlot,'EdgeColor',[0 1 0],'EdgeAlpha',0.5,...
        'tag',['tArtefactPeriod ' num2str(size(currentArtefactPeriods,1)+artefactNum)],...
        'yliminclude','off');
end

%add to table

[currentArtefactPeriods sInds] = sortrows([currentArtefactPeriods ;artefactPeriods]);
%get list of current patches, add new nes and sort by start time order to
%match table
hArtPatches = getappdata(hArtGui,'hArtPatches');
hArtPatches = [hArtPatches hNewArtPatches];
hArtPatches = hArtPatches(sInds);
setappdata(hArtGui,'hArtPatches',hArtPatches);

set(findobj('tag','tArtefactTable'),'Data',currentArtefactPeriods)
channelData.timesToExclude = currentArtefactPeriods;
setappdata(hArtGui,'channelData',channelData);
end
%% MISC FUNCTIONS

function [] = saveAndExit(src,eventData)
%function to save (to disk or workspace) and exit
hArtGui = getappdata(0,'hArtGui');
channelData = getappdata(hArtGui,'channelData');
%ask user to assign or save

choice = questdlg('Please Choose Data Export Method', ...
    'Export Data','Save to disk','Assign to Caller','Cancel','Cancel');
if strcmp(choice,'Save to disk')
    disp('save data to disk')
    [fileName,pathName,filterIndex] = uiputfile('.mat','Save Data',...
        ['-cleanedData.mat']);
    if filterIndex
         save([pathName fileName],'channelData')
    else
        errordlg('incorrect file name');
        return
    end
elseif strcmp(choice,'Assign to Caller')
    disp('assign to caller')
    assignin('caller','channelData',channelData)
else
    return;
    
end
set(hArtGui,'userdata','exported')
%closeMainWindow %TODO atuo exit?
end

function [] = closeMainWindow(varargin)
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
%hArtGui = getappdata(0,'hArtGui');
uiresume;
delete(gcf)
end