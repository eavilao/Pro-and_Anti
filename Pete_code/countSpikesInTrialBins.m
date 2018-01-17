function [countResults axHand] = countSpikesInTrialBins(trialStructure,trialsToPlot,timeWindow,pdfSaveName,unitNum,varargin)
%function to plot one trial of Erics data with spikes, eye data and time
%points of each event
 %this section is way out of date but the wholoe of this function was taken
 %from kunimatsuAnalysis.m. The idea is that it will take all trials and
 %then count the number of spikes (convert to rate to take into account differeing sacc lengthts) within given bins.
 %the bins we want for this time round are:
%1 - instruction - 0:150ms post fixation
%2 - presaccade - -200:0ms saccade onset
%3 - saccade - during saccade, need to find end points
%4 - post-saccade - 0:200ms post end of saccade, need to find endpoints
%
%the last 3 are all relative to sacade sowill do that first.

    trialTestTimes.instruction = [0 0.15]; %relative to fixation
    trialTestTimes.preSacc = [-0.2 0]; %relative saccade onset
    trialTestTimes.periSacc = [0 0.1]; %this is start and end of saccade (not really 100ms)
    trialTestTimes.postSacc = [0 0.2]; %this is relative to saccade offset
    
p = inputParser;
p.addParamValue('sigLevel',0.05,@(x) isnumeric(x));
p.parse(varargin{:});
%displayFlag = p.Results.displayFlag;
sigLevel = p.Results.sigLevel;
axHand = nan;
if isempty(pdfSaveName)
    doPlotFlag = 0;
else
    doPlotFlag = 1; %flag to control trial by trial plotting
end

numTrials = length(trialsToPlot);
stInit = cell(1,numTrials);
% kuniRes = struct('baselineRate',stInit,'typeInstRate',stInit,...
%     'dirInstRate',stInit,'saccadeRate',stInit,'trialInFile',stInit,...
%     'baselineCV2',stInit,'typeInstCV2',stInit,'dirInstCV2',stInit,'saccadeCV2',stInit);
countResults = struct('instructionBin',stInit,'preSaccBin',stInit,'periSaccBin',stInit,'postSaccBin',stInit,'conditionCode',stInit);
%loop over trials and create plot of the trial progression
for trialNum = 1:numTrials
    
    
   % kuniRes(trialNum).trialInFile = trialsToPlot(trialNum);
    
    thisTrial = trialStructure(trialsToPlot(trialNum));
    countResults(trialNum).conditionCode = thisTrial.conditionCode;
    %get relative timings
    targetAppears = thisTrial.bit2time-thisTrial.trialStart;
    goCue = targetAppears+0.1;
    targetHit = thisTrial.bit8time-thisTrial.trialStart;

    
    saccEndTime = thisTrial.saccadeTime+thisTrial.saccadeDuration/1000;
      preSaccSpikes = cellfun(@(x) x(x>(trialTestTimes.preSacc(1)+thisTrial.saccadeTime) & x<(trialTestTimes.preSacc(2)+thisTrial.saccadeTime)),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
    
    periSaccSpikes = cellfun(@(x) x(x>thisTrial.saccadeTime & x<saccEndTime),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
 
 postSaccSpikes = cellfun(@(x) x(x>(saccEndTime+trialTestTimes.postSacc(1)) & x<(saccEndTime+trialTestTimes.postSacc(2))),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);

instructionSpikes = cellfun(@(x) x(x>(trialTestTimes.instruction(1)) & x<(trialTestTimes.instruction(2))),thisTrial.alignedSpikes(unitNum),'uniformOutput',false);
   

%now convert all to firing rate in Hz
countResults(trialNum).instructionBin = numel([instructionSpikes{:}])/sum(trialTestTimes.instruction);
 countResults(trialNum).postSaccBin = numel([postSaccSpikes{:}])/sum(trialTestTimes.postSacc);
    countResults(trialNum).periSaccBin = numel([periSaccSpikes{:}])/(thisTrial.saccadeDuration/1000);
     countResults(trialNum).preSaccBin = numel([preSaccSpikes{:}])/sum(abs(trialTestTimes.preSacc));
   

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
        line([0.15 0.15],[0 200],'color','r','parent',hAx3,...
            'linestyle','-.') %instruction
        line([thisTrial.saccadeTime thisTrial.saccadeTime],[0 200],'color','g','parent',hAx3,...
            'linestyle','-') %saccade start
         line([saccEndTime saccEndTime],[0 200],'color','g','parent',hAx3,...
            'linestyle','-') %saccade end
        line([thisTrial.saccadeTime thisTrial.saccadeTime],[0 200],'color','g','parent',hAx3,...
            'linestyle','-') %saccade start
        
        if ~isempty(targetHit)
        line([targetHit targetHit],[0 200],'color','m','parent',hAx3,...
            'linestyle','-.','yliminclude','off') %target hit according to stim computer
        end
        ylim([0 150])
        
        
       
        
        %now plot lines representing the rate at each time point
        line([trialTestTimes.instruction(1) trialTestTimes.instruction(2)],repmat(countResults(trialNum).instructionBin,1,2),'parent',hAx3)
        line([trialTestTimes.preSacc(1)+thisTrial.saccadeTime trialTestTimes.preSacc(2)+thisTrial.saccadeTime],repmat(countResults(trialNum).preSaccBin,1,2),'parent',hAx3)
        line([thisTrial.saccadeTime saccEndTime],repmat(countResults(trialNum).periSaccBin,1,2),'parent',hAx3)
        line([saccEndTime+trialTestTimes.postSacc(1) saccEndTime+trialTestTimes.postSacc(2)],repmat(countResults(trialNum).postSaccBin,1,2),'parent',hAx3)
        axHand = [hAx3 hAx2 hAx1];
        linkaxes(axHand,'x')
        xlim([-0.5 1])
        %save figure and close
        export_fig(pdfSaveName, '-pdf','-append', hF);
        delete(hF)
        
    end
   
    
   
end

% 
% kuniCell = squeeze(struct2cell(kuniRes));
% 
% %perform paired ttest for each trial, for both instruction
% %periods and the saccade period
% %[hTypeInst pTypeInst] = ttest([kuniCell{1,:}],[kuniCell{2,:}]);
% %[hDirInst pDirInst] = ttest([kuniCell{1,:}],[kuniCell{3,:}]);
% %[hSacc pSacc] = ttest([kuniCell{1,:}],[kuniCell{4,:}]);
% %now switched to nonparametric version
% [pTypeInst hTypeInst] = signrank([kuniCell{1,:}],[kuniCell{2,:}],'alpha',sigLevel);
% [pDirInst hDirInst] = signrank([kuniCell{1,:}],[kuniCell{3,:}],'alpha',sigLevel);
% [pSacc hSacc] = signrank([kuniCell{1,:}],[kuniCell{4,:}],'alpha',sigLevel);
% 
% 
% 
% %calcualte mean and sem for errorbar plot - no longer needed as
% %ditribuaiton is plotted
% muRates = nanmean(cell2mat(kuniCell),2);
% semRates = nanstd(cell2mat(kuniCell),0,2)./sqrt(size(kuniCell,2));
% 
% 
% kuniStats.muRates = muRates;
% kuniStats.semRates = semRates;
% kuniStats.pValues = [pTypeInst pDirInst pSacc];
% 
% 
% pTypeInstCV2 = nan;
% hTypeInstCV2 = nan;
% pDirInstCV2 = nan;
% hDirInstCV2 = nan;
% pSaccCV2 = nan;
% hSaccCV2 = nan;
% %now alos get the signifcantly different CV2 periods and return
% kuniTemp = cell2mat(kuniCell);
% 
% 
% if sum(~isnan(sum(kuniTemp([6 7],:))))>1
%     [pTypeInstCV2 hTypeInstCV2] = signrank([kuniCell{6,:}],[kuniCell{7,:}],'alpha',sigLevel);
% end
%     
% if sum(~isnan(sum(kuniTemp([6 8],:))))>1
%     [pDirInstCV2 hDirInstCV2] = signrank([kuniCell{6,:}],[kuniCell{8,:}],'alpha',sigLevel);
% end
% if sum(~isnan(sum(kuniTemp([6 9],:))))>1
%     [pSaccCV2 hSaccCV2] = signrank([kuniCell{6,:}],[kuniCell{9,:}],'alpha',sigLevel);
% end
% 
% kuniStats.pValuesCV2 = [pTypeInstCV2 pDirInstCV2 pSaccCV2];
% 
% %for each period mark with 0 for nonsig, 1 for sig increase and -1 for sig
% %decrease %TODO must be a nicer way of doing this but I'm tired
% 
% 
% goneDown = (muRates(2:4)<muRates(1));
% goneUp = (muRates(2:4)>muRates(1));
% sigChange = (kuniStats.pValues<sigLevel)';
% 
% modulationDirection = nan(1,3);
% modulationDirection(goneUp & sigChange) = 1;
% modulationDirection(goneDown & sigChange) = -1;
% modulationDirection(isnan(modulationDirection)) = 0;
% kuniStats.modulationDirection = modulationDirection;
% %output actual rates for all trials to enable testing against other
% %types %TODO can probably be replaced by using squeeze outside fcn
% kuniStats.saccRates = [kuniCell{4,:}];
% kuniStats.typeInstRates = [kuniCell{2,:}];
% kuniStats.dirInstRates = [kuniCell{3,:}];
% kuniStats.baselineRates = [kuniCell{1,:}];
% 
% kuniStats.saccCV2s = [kuniCell{9,:}];
% kuniStats.typeInstCV2s = [kuniCell{7,:}];
% kuniStats.dirInstCV2s = [kuniCell{8,:}];
% kuniStats.baselineCV2s = [kuniCell{6,:}];
% 
% %errorbar plot of the 4 time points
% forDistPlot =  cell2mat(kuniCell);
% 
% 
% 
% hFb = figure;
% hAb = subplot(2,1,1,'parent',hFb); %for firing rate
% handles = distributionPlot(hAb,forDistPlot(1:4,:)','xNames',{'base','typeInst','dirInst','sacc'});
% axHand{1} = handles{3};
% set(get(handles{3},'title'),'string','Rate')
% sigChange = (kuniStats.pValues<sigLevel)';
% curYLim = get(handles{3},'ylim');
% yLevel = 1*curYLim(2);
% yPoints = nan(1,4);
% xPoints = 1:4;
% yPoints(find(sigChange)+1) = yLevel;
% line(xPoints,yPoints,'linestyle','none','color','r','marker','*',...
%     'markersize',15,'parent',handles{3})
% hAb2 = subplot(2,1,2,'parent',hFb); %for CV2
% handles = distributionPlot(hAb2,forDistPlot(6:9,:)','xNames',{'base','typeInst','dirInst','sacc'});
% axHand{2} = handles{3};
% set(get(handles{3},'title'),'string','CV2')
% set(handles{3},'ylim',[0 1.5])
% sigChange = (kuniStats.pValuesCV2<sigLevel)';
% curYLim = get(handles{3},'ylim');
% yLevel = 1*curYLim(2);
% yPoints = nan(1,4);
% xPoints = 1:4;
% yPoints(find(sigChange)+1) = yLevel;
% line(xPoints,yPoints,'linestyle','none','color','r','marker','*',...
%     'markersize',15,'parent',handles{3})

%now also makr those that are significantly different from baseline with stars

 
           
