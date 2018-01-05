function [kuniRes kuniStats axHand] = kunimatsuAnalysisNew(trialStructure,trialsToPlot,timeWindow,pdfSaveName,unitNum,varargin)
%function to plot one trial of Erics data with spikes, eye data and time
%points of each event
%updated to deal with the 2016 kunimatsu paper with a 400ms window starting 300ms before the target appears 

p = inputParser;
p.addParamValue('sigLevel',0.05,@(x) isnumeric(x));
p.parse(varargin{:});
%displayFlag = p.Results.displayFlag;
sigLevel = p.Results.sigLevel;

if isempty(pdfSaveName)
 doPlotFlag = 0;
else
    doPlotFlag = 1; %flag to control trial by trial plotting
end

numTrials = length(trialsToPlot);
stInit = cell(1,numTrials);
kuniRes = struct('baselineRate',stInit,'typeInstRate',stInit,...
    'dirInstRate',stInit,'saccadeRate',stInit,'trialInFile',stInit,...
    'baselineCV2',stInit,'typeInstCV2',stInit,'dirInstCV2',stInit,'saccadeCV2',stInit);

%loop over trials and create plot of the trial progression
for trialNum = 1:numTrials
    
    
    kuniRes(trialNum).trialInFile = trialsToPlot(trialNum);
    
    thisTrial = trialStructure(trialsToPlot(trialNum));
    
    %get relative timings
    targetAppears = thisTrial.bit2time-thisTrial.trialStart;
    goCue = targetAppears+0.1;
    targetHit = thisTrial.bit8time-thisTrial.trialStart;

    
     %kunimatsu analysis %TODO this will give errors if no saccade present
    %timings are in s
    %baseline is -0.4:0  relative to when he enters fixation (time 0)
    %typeInstruction is -0.1:0 relative to  targetAppears
    %directionInstruction is ):0.1 relative to  targetAppears
    %saccadePeriod is -0.05:0.05 relative to saccade onset
    
   
%  Old   Times to Analyse:
% Inst1 - bit6:bit6+300ms - he knows what type of trial but not which direction
% Inst2 - bit2:bit2+100ms - he knows target direction but can't move
% Sacc - onset-25ms:onset+25ms

    trialTestTimes.baseLine = [-0.4 0]; %bit6-400ms:bit6
    trialTestTimes.typeInst = [-0.3 0.1]; %relative to bit 2
    trialTestTimes.dirInst = [0 0.1]; %this is relative to bit2 
    trialTestTimes.saccPeriod = [-0.025 0.25]; %this is relative to saccade onset
   
    
     %calcualte spike rate in each period  and mark on bottom axis
    baselineSpikes = cellfun(@(x) x(x>trialTestTimes.baseLine(1) & x<trialTestTimes.baseLine(2)),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    typeInstSpikes = cellfun(@(x) x(x>(targetAppears+trialTestTimes.typeInst(1)) & x<(targetAppears+trialTestTimes.typeInst(2))),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    dirInstSpikes = cellfun(@(x) x(x>(targetAppears+trialTestTimes.dirInst(1)) & x<(targetAppears+trialTestTimes.dirInst(2))),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    saccadeSpikes = cellfun(@(x) x(x>(thisTrial.saccadeTime+trialTestTimes.saccPeriod(1)) & x<(thisTrial.saccadeTime+trialTestTimes.saccPeriod(2))),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    
    
    %instead of calculating CV2 like this in isolation it is best to
    %calculate the CV2 of the whole trace then get the values for those
    %spikes within range
%     kuniRes(trialNum).baselineCV2 = nanmedian(calculateCV2([baselineSpikes{:}]));
%     kuniRes(trialNum).typeInstCV2 = nanmedian(calculateCV2([typeInstSpikes{:}]));
%     kuniRes(trialNum).dirInstCV2 = nanmedian(calculateCV2([dirInstSpikes{:}]));
%     kuniRes(trialNum).saccadeCV2 = nanmedian(calculateCV2([saccadeSpikes{:}]));
    
    kuniRes(trialNum).baselineCV2 = nan;
    kuniRes(trialNum).typeInstCV2 = nan;
    kuniRes(trialNum).dirInstCV2 = nan;
    kuniRes(trialNum).saccadeCV2 = nan;


    thisCV2trace = [calculateCV2(thisTrial.alignedSpikes{unitNum}); nan]; %get the CV2 of teh whole trace and then add a nan to the end so it is the same length as the spikes
    baselineSpikeInds = cellfun(@(x) find(x>trialTestTimes.baseLine(1) & x<trialTestTimes.baseLine(2)),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    typeInstSpikesInds = cellfun(@(x) find(x>(targetAppears+trialTestTimes.typeInst(1)) & x<(targetAppears+trialTestTimes.typeInst(2))),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    dirInstSpikesInds = cellfun(@(x) find(x>(targetAppears+trialTestTimes.dirInst(1)) & x<(targetAppears+trialTestTimes.dirInst(2))),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    saccadeSpikesInds = cellfun(@(x) find(x>(thisTrial.saccadeTime+trialTestTimes.saccPeriod(1)) & x<(thisTrial.saccadeTime+trialTestTimes.saccPeriod(2))),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    
    
    if length(thisCV2trace)>2
    if ~cellfun('isempty',baselineSpikeInds) 
    kuniRes(trialNum).baselineCV2 = nanmedian(thisCV2trace([baselineSpikeInds{:}]));
    end
    if ~cellfun('isempty',typeInstSpikesInds)
     kuniRes(trialNum).typeInstCV2 = nanmedian(thisCV2trace([typeInstSpikesInds{:}]));
    end
    if ~cellfun('isempty',dirInstSpikesInds)
    kuniRes(trialNum).dirInstCV2 = nanmedian(thisCV2trace([dirInstSpikesInds{:}]));
    end
    if ~cellfun('isempty',saccadeSpikesInds)
    kuniRes(trialNum).saccadeCV2 = nanmedian(thisCV2trace([saccadeSpikesInds{:}]));
    end
    end
    
    
    %now convert all to firing rate in Hz
    kuniRes(trialNum).baselineRate = length(baselineSpikes{:})/(trialTestTimes.baseLine(2)-trialTestTimes.baseLine(1));
    kuniRes(trialNum).typeInstRate = length(typeInstSpikes{:})/(trialTestTimes.typeInst(2)-trialTestTimes.typeInst(1));
    kuniRes(trialNum).dirInstRate = length(dirInstSpikes{:})/(trialTestTimes.dirInst(2)-trialTestTimes.dirInst(1));
    kuniRes(trialNum).saccadeRate = length(saccadeSpikes{:})/(trialTestTimes.saccPeriod(2)-trialTestTimes.saccPeriod(1));
    
    if doPlotFlag
        
        hF = figure;
        hAx1 = subplot(3,1,1); %top plot for eye movement
        hAx2 = subplot(3,1,2); %for spikes
        hAx3 = subplot(3,1,3); %trial timings
        
        set([hAx1 hAx2 hAx3],'xlim',[-timeWindow(1) timeWindow(2)])
        
        %plot both eye traces %TODO also plot target position
        eyeT = linspace(-timeWindow(1),timeWindow(2),length(thisTrial.eyePositionX));
        
        line(eyeT,thisTrial.eyePositionX,'color','r','linewidth',2,'parent',hAx1)
        line(eyeT,thisTrial.eyePositionY,'color','b','linewidth',2,'parent',hAx1)
        legend('X','Y')
        
        %plot spike timing for both units %TODO include option for more units
        line(thisTrial.alignedSpikes{1},1*ones(1,length(thisTrial.alignedSpikes{1})),'color','k','parent',hAx2,...
            'linestyle','none','marker','d','MarkerFaceColor','k')
        line(thisTrial.alignedSpikes{2},2*ones(1,length(thisTrial.alignedSpikes{2})),'color','r','parent',hAx2,...
            'linestyle','none','marker','d','MarkerFaceColor','r')
        legend('unit 1','unit 2')
        ylim([0.5 2.5])
        
        
        
        %now mark the various trial timing aspects
        line([0 0],[0 200],'color','k','parent',hAx3,...
            'linestyle','-') %trial start
        line([targetAppears targetAppears],[0 200],'color','r','parent',hAx3,...
            'linestyle','-') %target apears
        line([goCue goCue],[0 200],'color','g','parent',hAx3,...
            'linestyle','-') %go Cue
        
        if ~isempty(targetHit)
        line([targetHit targetHit],[0 200],'color','g','parent',hAx3,...
            'linestyle','-.','yliminclude','off') %target hit according to stim computer
        end
        ylim([0 150])
        
        
        %saccade time
        if ~isnan(thisTrial.saccadeTime)
            
            line(thisTrial.saccadeTime,1,'color','g','parent',hAx1,...
                'linestyle','none','marker','d','MarkerFaceColor','g','markersize',5)
            
        end
        
        %now plot lines representing the rate at each time point
        line([trialTestTimes.baseLine(1) trialTestTimes.baseLine(2)],repmat(kuniRes(trialNum).baselineRate,1,2),'parent',hAx3)
        line([trialTestTimes.typeInst(1) trialTestTimes.typeInst(2)],repmat(kuniRes(trialNum).typeInstRate,1,2),'parent',hAx3)
        line([targetAppears+trialTestTimes.dirInst(1) targetAppears+trialTestTimes.dirInst(2)],repmat(kuniRes(trialNum).dirInstRate,1,2),'parent',hAx3)
        line([thisTrial.saccadeTime+trialTestTimes.saccPeriod(1) thisTrial.saccadeTime+trialTestTimes.saccPeriod(2)],repmat(kuniRes(trialNum).saccadeRate,1,2),'parent',hAx3)
        %save figure and close
        export_fig(pdfSaveName, '-pdf','-append', hF);
        delete(hF)
        
    end
   
    
   
end


kuniCell = squeeze(struct2cell(kuniRes));

%perform paired ttest for each trial, for both instruction
%periods and the saccade period
%[hTypeInst pTypeInst] = ttest([kuniCell{1,:}],[kuniCell{2,:}]);
%[hDirInst pDirInst] = ttest([kuniCell{1,:}],[kuniCell{3,:}]);
%[hSacc pSacc] = ttest([kuniCell{1,:}],[kuniCell{4,:}]);
%now switched to nonparametric version
[pTypeInst hTypeInst] = signrank([kuniCell{1,:}],[kuniCell{2,:}],'alpha',sigLevel);
[pDirInst hDirInst] = signrank([kuniCell{1,:}],[kuniCell{3,:}],'alpha',sigLevel);
[pSacc hSacc] = signrank([kuniCell{1,:}],[kuniCell{4,:}],'alpha',sigLevel);



%calcualte mean and sem for errorbar plot - no longer needed as
%ditribuaiton is plotted
muRates = nanmean(cell2mat(kuniCell),2);
semRates = nanstd(cell2mat(kuniCell),0,2)./sqrt(size(kuniCell,2));


kuniStats.muRates = muRates;
kuniStats.semRates = semRates;
kuniStats.pValues = [pTypeInst pDirInst pSacc];


pTypeInstCV2 = nan;
hTypeInstCV2 = nan;
pDirInstCV2 = nan;
hDirInstCV2 = nan;
pSaccCV2 = nan;
hSaccCV2 = nan;
%now alos get the signifcantly different CV2 periods and return
kuniTemp = cell2mat(kuniCell);


if sum(~isnan(sum(kuniTemp([6 7],:))))>1
    [pTypeInstCV2 hTypeInstCV2] = signrank([kuniCell{6,:}],[kuniCell{7,:}],'alpha',sigLevel);
end
    
if sum(~isnan(sum(kuniTemp([6 8],:))))>1
    [pDirInstCV2 hDirInstCV2] = signrank([kuniCell{6,:}],[kuniCell{8,:}],'alpha',sigLevel);
end
if sum(~isnan(sum(kuniTemp([6 9],:))))>1
    [pSaccCV2 hSaccCV2] = signrank([kuniCell{6,:}],[kuniCell{9,:}],'alpha',sigLevel);
end

kuniStats.pValuesCV2 = [pTypeInstCV2 pDirInstCV2 pSaccCV2];

%for each period mark with 0 for nonsig, 1 for sig increase and -1 for sig
%decrease %TODO must be a nicer way of doing this but I'm tired


goneDown = (muRates(2:4)<muRates(1));
goneUp = (muRates(2:4)>muRates(1));
sigChange = (kuniStats.pValues<sigLevel)';

modulationDirection = nan(1,3);
modulationDirection(goneUp & sigChange) = 1;
modulationDirection(goneDown & sigChange) = -1;
modulationDirection(isnan(modulationDirection)) = 0;
kuniStats.modulationDirection = modulationDirection;
%output actual rates for all trials to enable testing against other
%types %TODO can probably be replaced by using squeeze outside fcn
kuniStats.saccRates = [kuniCell{4,:}];
kuniStats.typeInstRates = [kuniCell{2,:}];
kuniStats.dirInstRates = [kuniCell{3,:}];
kuniStats.baselineRates = [kuniCell{1,:}];

kuniStats.saccCV2s = [kuniCell{9,:}];
kuniStats.typeInstCV2s = [kuniCell{7,:}];
kuniStats.dirInstCV2s = [kuniCell{8,:}];
kuniStats.baselineCV2s = [kuniCell{6,:}];

%errorbar plot of the 4 time points
forDistPlot =  cell2mat(kuniCell);



hFb = figure;
hAb = subplot(2,1,1,'parent',hFb); %for firing rate
handles = distributionPlot(hAb,forDistPlot(1:4,:)','xNames',{'base','typeInst','dirInst','sacc'});
axHand{1} = handles{3};
set(get(handles{3},'title'),'string','Rate')
sigChange = (kuniStats.pValues<sigLevel)';
curYLim = get(handles{3},'ylim');
yLevel = 1*curYLim(2);
yPoints = nan(1,4);
xPoints = 1:4;
yPoints(find(sigChange)+1) = yLevel;
line(xPoints,yPoints,'linestyle','none','color','r','marker','*',...
    'markersize',15,'parent',handles{3})
hAb2 = subplot(2,1,2,'parent',hFb); %for CV2
handles = distributionPlot(hAb2,forDistPlot(6:9,:)','xNames',{'base','typeInst','dirInst','sacc'});
axHand{2} = handles{3};
set(get(handles{3},'title'),'string','CV2')
set(handles{3},'ylim',[0 1.5])
sigChange = (kuniStats.pValuesCV2<sigLevel)';
curYLim = get(handles{3},'ylim');
yLevel = 1*curYLim(2);
yPoints = nan(1,4);
xPoints = 1:4;
yPoints(find(sigChange)+1) = yLevel;
line(xPoints,yPoints,'linestyle','none','color','r','marker','*',...
    'markersize',15,'parent',handles{3})

%now also makr those that are significantly different from baseline with stars

 
           
