function newTrialStructure = appendRastersToTrials(trialData,spikeTimes,saccadeDetectOut,eyeData,rasterTimesWindow,pdfFileName)
%function to append the eye data and spike times to the trialData structure
%for Eric's saccade experiemnt


%INPUTS:
%rasterTimesWindow = [0.5 4]; %times in seconds either side of trail start to get spikes for

reactionTimeLimit = [-.15 0.5]; %how many seconds after the go cue [atLeast atMost] to count as valid reaction times
eyePosition = hypot(saccadeDetectOut.traces.position.raw(:,1),saccadeDetectOut.traces.position.raw(:,2));

doTrialPlot = 0; %flag to control if trials are plotted individually

%create copy of trial data
newTrialStructure = trialData;

%get start tiems of all saccades
allSaccStarts = saccadeDetectOut.saccades.landmarks.v(:,1)/1000;

numTrials = size(trialData,2);
 
%  if numTrials ==81
%      numTrials = [112 137 149 150 157]-83
%  else
%      numTrials = 1:numTrials
%  end
 
%loop over every trial
for trialNum =     1:numTrials
    
    disp(['processing trial number ' num2str(trialNum) ' of ' num2str(numTrials)])
    %get start and end time of this trial
    thisTrialStart = trialData(trialNum).trialStart;
    thisTrialEnd = trialData(trialNum).trialEnd;
    goCueTime = trialData(trialNum).bit3time;
    %if it is a 'valid' trial (ie both a start, end and go cue) then check if tehre
    %is a saccade(s) within it
    if ~isempty(thisTrialStart) && ~isempty(thisTrialEnd) && ~isempty(goCueTime)
        
       %also record reward (bit4) time if there

        
        rewardBitIndex = find(trialData(trialNum).trialBits(5,:)==1);
        
        if ~isempty(rewardBitIndex) %TODO currently takes first reward only if more than 1
            rewardBitTime = trialData(trialNum).trialBits(1,rewardBitIndex(1));
            newTrialStructure(trialNum).relativeRewardTime = rewardBitTime-thisTrialStart;
        else 
            newTrialStructure(trialNum).relativeRewardTime = nan;
        end
        %first extract all spikes within window of this trial start(for both units)
       theseSpikes = cellfun(@(x) x(x>(thisTrialStart-rasterTimesWindow(1)) & x<(thisTrialStart+rasterTimesWindow(2))),spikeTimes,'uniformOutput',false); 
       %now align them to the start of the trial
        newTrialStructure(trialNum).alignedSpikes = cellfun(@(x) (x-thisTrialStart),theseSpikes,'uniformoutput',false);
       
        
        %now get eye movement aligned to start of trial
        %TODO currently we drop trials that are too clsoe to the beginning
        %or end of file but should actually keep them
        
        eyeStartSample = floor((thisTrialStart-rasterTimesWindow(1))*eyeData.eyeFsCal);
        eyeSegmentLength = ceil(sum(rasterTimesWindow)*eyeData.eyeFsCal);
        
        if eyeStartSample>0 && (eyeStartSample+eyeSegmentLength)<length(eyeData.eyeXcal)
            
            newTrialStructure(trialNum).eyePositionX = eyeData.eyeXcal(eyeStartSample:eyeStartSample+eyeSegmentLength-1);
            newTrialStructure(trialNum).eyePositionY = eyeData.eyeYcal(eyeStartSample:eyeStartSample+eyeSegmentLength-1);
            %new traces to add to structure
            newTrialStructure(trialNum).eyeVelocityX =  eyeData.eyeVelX(eyeStartSample:eyeStartSample+eyeSegmentLength-1);
            newTrialStructure(trialNum).eyeVelocityY =  eyeData.eyeVelY(eyeStartSample:eyeStartSample+eyeSegmentLength-1);
            newTrialStructure(trialNum).eyeVelocityV =  eyeData.eyeVelV(eyeStartSample:eyeStartSample+eyeSegmentLength-1);
            newTrialStructure(trialNum).eyeAngle = eyeData.eyeAng(eyeStartSample:eyeStartSample+eyeSegmentLength-1);
        else
            newTrialStructure(trialNum).eyePositionX = NaN(1,eyeSegmentLength);
            newTrialStructure(trialNum).eyePositionY = NaN(1,eyeSegmentLength);
            %new traces to add to structure
            newTrialStructure(trialNum).eyeVelocityX =  NaN(1,eyeSegmentLength);
            newTrialStructure(trialNum).eyeVelocityY =  NaN(1,eyeSegmentLength);
            newTrialStructure(trialNum).eyeVelocityV =  NaN(1,eyeSegmentLength);
            newTrialStructure(trialNum).eyeAngle = NaN(1,eyeSegmentLength);
            
        end
        
        
        
        %now find all teh detected saccades within this trial
         thisTrialSaccadeNumbers = find(allSaccStarts>thisTrialStart & allSaccStarts<thisTrialEnd);
        
         %now get all the times of these saccades and align them to trial
         %start
         thisTrialSaccadeStarts = allSaccStarts(thisTrialSaccadeNumbers)-thisTrialStart;
         
       %now trim to saccades that are after the go cue but within the
       %reactionTimeLimit
       thisTrialValidSaccadeNumbers = thisTrialSaccadeNumbers(allSaccStarts(thisTrialSaccadeNumbers)>goCueTime+reactionTimeLimit(1)...
           & allSaccStarts(thisTrialSaccadeNumbers)<goCueTime+reactionTimeLimit(2));
    
        thisTrialValidSaccadeStarts = allSaccStarts(thisTrialValidSaccadeNumbers)-thisTrialStart;
        relativeGoCueTime = goCueTime-thisTrialStart;
        %TODO now do plot of trial to cehck saccade detection
        if doTrialPlot
        
        
            hF = figure;
            hAx1 = subplot(3,1,1); %top plot for eye movement
            hAx2 = subplot(3,1,2); %for spikes
            hAx3 = subplot(3,1,3); %trial timings
            
            set([hAx1 hAx2 hAx3],'xlim',[-rasterTimesWindow(1) rasterTimesWindow(2)])
            
            %plot both eye traces %TODO also plot target position
            eyeT = linspace(-rasterTimesWindow(1),rasterTimesWindow(2),eyeSegmentLength);
            line(eyeT,newTrialStructure(trialNum).eyePositionX,'color','r','linewidth',2,'parent',hAx1)
            line(eyeT,newTrialStructure(trialNum).eyePositionY,'color','b','linewidth',2,'parent',hAx1)
            legend('X','Y')
        
            %plot the detected saccade times
            line(thisTrialSaccadeStarts,1*ones(1,length(thisTrialSaccadeStarts)),'color','k','parent',hAx1,...
                'linestyle','none','marker','s','MarkerFaceColor','k')
            %on top mark the valid saccades
            line(thisTrialValidSaccadeStarts,1*ones(1,length(thisTrialValidSaccadeStarts)),'color','g','parent',hAx1,...
                'linestyle','none','marker','s','MarkerFaceColor','g','markersize',2)
            
            %plot spike timing for both units %TODO include option for more units
            line(newTrialStructure(trialNum).alignedSpikes{1},1*ones(1,length(newTrialStructure(trialNum).alignedSpikes{1})),'color','k','parent',hAx2,...
                'linestyle','none','marker','d','MarkerFaceColor','k')
            line(newTrialStructure(trialNum).alignedSpikes{2},2*ones(1,length(newTrialStructure(trialNum).alignedSpikes{2})),'color','r','parent',hAx2,...
                'linestyle','none','marker','d','MarkerFaceColor','r')
            legend('unit 1','unit 2')
            ylim([0.5 2.5])
            
            
            
            
            %get the time reported by stim computer as a hit (if exists)
            if isempty(newTrialStructure(trialNum).bit8time)
                targetHit = NaN;
            else
                targetHit = newTrialStructure(trialNum).bit8time-thisTrialStart; 
            end
            %now mark the various trial timing aspects
            line([0 0],[0 1],'color','k','parent',hAx3,...
                'linestyle','-') %trial start
            line([relativeGoCueTime relativeGoCueTime]-0.1,[0 1],'color','m','parent',hAx3,...
                'linestyle','-') %target apears
            line([relativeGoCueTime relativeGoCueTime],[0 1],'color','r','parent',hAx3,...
                'linestyle','-') %go cue
            line([targetHit targetHit],[0 1],'color','g','parent',hAx3,...
                'linestyle','-') %target hit according to stim computer
            ylim([0 1])
            
            suptitle(['Trial ' num2str(trialNum) ', condition ' num2str(newTrialStructure(trialNum).conditionCode)])
            linkaxes([hAx1 hAx2 hAx3],'x')
            try 
                export_fig(pdfFileName, '-pdf','-append', hF);
            catch
                 export_fig(pdfFileName, '-pdf', hF);           
            end
            delete(hF)
        end
        
                  
        
        
        if isempty(thisTrialValidSaccadeNumbers)
            
            newTrialStructure(trialNum).goCueTime = relativeGoCueTime;
            %if no valid saccades then everything should be Nan
            newTrialStructure(trialNum).saccadeTime = NaN;
            newTrialStructure(trialNum).reactionTime = NaN;
            newTrialStructure(trialNum).saccadeAmplitude = NaN;
            newTrialStructure(trialNum).saccadeDuration = NaN;
            newTrialStructure(trialNum).saccadePeakVel = NaN;
            newTrialStructure(trialNum).saccadeOriginX = NaN;
            newTrialStructure(trialNum).saccadeOriginY = NaN;
            newTrialStructure(trialNum).saccadeLandmarks = nan(1,6);
            newTrialStructure(trialNum).saccadeTailDuration = NaN;
            newTrialStructure(trialNum).saccadeTailAmpltiude = NaN;
             newTrialStructure(trialNum).saccadeAngle = NaN;
        else
            %TODO for now we only take the first valid saccade but we should do
            %something with the rest
            if length(thisTrialValidSaccadeNumbers)>1
                %if more than 1 saccade choose the one with the one closest
                %to the targetHit time
                
                % throw out sac to fix
                 
                for i = thisTrialValidSaccadeNumbers'
                     saccadeOnAndOffset  = eyePosition(int32(saccadeDetectOut.saccades.landmarks.v(i,[1,4]))); % eye position at on and offset saccsade
                     if saccadeOnAndOffset(1) > saccadeOnAndOffset(2) % if offset saccade closer to fixation than onset: remove
                         thisTrialValidSaccadeNumbers(thisTrialValidSaccadeNumbers==i)=[];
                     end
                end
               
                     
                % check how many sacs left
                if numel(thisTrialValidSaccadeNumbers) == 1
                   saccInd = 1;
                    % do nothing 
                % Find sac closest to gocue/bit8
                
                % if this proves unreliable find biggest sac
                
                elseif isempty(newTrialStructure(trialNum).bit8time)
                    %if incorrect trial (therefore no bit8 for correct) then
                    %use biggest saccade
                    
                    relativeSaccadeTimes = allSaccStarts(thisTrialValidSaccadeNumbers) - goCueTime ;
                    
                    % trow out saccades with reaction time bigger than 
                    
                    if numel(relativeSaccadeTimes) > 1
                    relativeSaccadeTimes(relativeSaccadeTimes>reactionTimeLimit(2))=[];
                    end 
                    thisTrialSaccadeOnAndOffset  = eyePosition(saccadeDetectOut.saccades.landmarks.v(thisTrialValidSaccadeNumbers,[1,4]));
                    
                   % find first sac after go-cue! 
                   
                   [val saccInd] = min(relativeSaccadeTimes);
                   
%                    [val2 saccInd2]    =  max(abs(thisTrialSaccadeOnAndOffset(:,2) - thisTrialSaccadeOnAndOffset(:,1))); % biggest saccade!

                   % biggest saccade is not the first one after goCue 
                                 
                    %[val saccInd] = max(saccadeDetectOut.saccades.amplitude(thisTrialValidSaccadeNumbers,1));
                else
                    [val saccInd] = min(abs(allSaccStarts(thisTrialValidSaccadeNumbers)-newTrialStructure(trialNum).bit8time));
                    
                end
                % [val saccInd] = max(saccadeDetectOut.saccades.amplitude(thisTrialValidSaccadeNumbers,1));
            else
                saccInd = 1;
            end
            thisSaccadeStartTime = allSaccStarts(thisTrialValidSaccadeNumbers(saccInd));
            
            %we need a final check that there is no eye movement in the
            %baselines
            %find index of saccade in eye trace. then check from 20ms before
            %to 100ms before that bothX and Y positon of the ey stay within
            %1 degree
            saccadeTimeIndex = round((rasterTimesWindow(1)+thisSaccadeStartTime-thisTrialStart)*eyeData.eyeFsCal);
            if saccadeTimeIndex<=eyeSegmentLength
                if  max(abs(newTrialStructure(trialNum).eyePositionX(saccadeTimeIndex-100:saccadeTimeIndex-20)))>1.5 ...
                    || max(abs(newTrialStructure(trialNum).eyePositionY(saccadeTimeIndex-100:saccadeTimeIndex-20)))>1.5
                
                    newTrialStructure(trialNum).goCueTime = relativeGoCueTime;
                    %if no valid saccades then everything should be Nan
                    newTrialStructure(trialNum).saccadeTime = NaN;
                    newTrialStructure(trialNum).reactionTime = NaN;
                    newTrialStructure(trialNum).saccadeAmplitude = NaN;
                    newTrialStructure(trialNum).saccadeDuration = NaN;
                    newTrialStructure(trialNum).saccadePeakVel = NaN;
                    newTrialStructure(trialNum).saccadeOriginX = NaN;
                    newTrialStructure(trialNum).saccadeOriginY = NaN;
                    newTrialStructure(trialNum).saccadeLandmarks = nan(1,6);
                    newTrialStructure(trialNum).saccadeTailDuration = NaN;
                    newTrialStructure(trialNum).saccadeTailAmpltiude = NaN;
                    newTrialStructure(trialNum).saccadeAngle = NaN;
                else
                    %save all of teh relevant saccade parameters here
                    newTrialStructure(trialNum).saccadeAmplitude = saccadeDetectOut.saccades.amplitude(thisTrialValidSaccadeNumbers(saccInd),1);
                    newTrialStructure(trialNum).saccadeAmplitudeX = saccadeDetectOut.saccades.amplitude(thisTrialValidSaccadeNumbers(saccInd),2);
                    newTrialStructure(trialNum).saccadeAmplitudeY = saccadeDetectOut.saccades.amplitude(thisTrialValidSaccadeNumbers(saccInd),3);
                    newTrialStructure(trialNum).saccadeOriginX = saccadeDetectOut.saccades.origin(thisTrialValidSaccadeNumbers(saccInd),1);
                    newTrialStructure(trialNum).saccadeOriginY = saccadeDetectOut.saccades.origin(thisTrialValidSaccadeNumbers(saccInd),2);
                    newTrialStructure(trialNum).saccadeDuration = saccadeDetectOut.saccades.duration(thisTrialValidSaccadeNumbers(saccInd));
                    newTrialStructure(trialNum).saccadePeakVel = saccadeDetectOut.saccades.peakV.v(thisTrialValidSaccadeNumbers(saccInd));
                    %the following stats will all be relative to the trialStart time
                    newTrialStructure(trialNum).goCueTime = relativeGoCueTime;
                    newTrialStructure(trialNum).saccadeTime = thisSaccadeStartTime-thisTrialStart;
                    newTrialStructure(trialNum).reactionTime = newTrialStructure(trialNum).saccadeTime-(newTrialStructure(trialNum).goCueTime-0.1);
                    newTrialStructure(trialNum).saccadeLandmarks = saccadeDetectOut.saccades.landmarks.v(thisTrialValidSaccadeNumbers(saccInd),:);
                    %add saccade landmarks and saccade tail paramters
                    newTrialStructure(trialNum).saccadeTailDuration = saccadeDetectOut.saccades.tail.duration(thisTrialValidSaccadeNumbers(saccInd));
                    newTrialStructure(trialNum).saccadeTailAmpltiude = saccadeDetectOut.saccades.tail.amplitude(thisTrialValidSaccadeNumbers(saccInd));
                    %TODO also add end point angle
                     xEnd =   newTrialStructure(trialNum).saccadeAmplitudeX+ newTrialStructure(trialNum).saccadeOriginX;
                    yEnd =  newTrialStructure(trialNum).saccadeAmplitudeY+ newTrialStructure(trialNum).saccadeOriginY;
                    newTrialStructure(trialNum).saccadeAngle = atan2(yEnd,xEnd);
                end
            else
                newTrialStructure(trialNum).goCueTime = relativeGoCueTime;
                %if no valid saccades then everything should be Nan
                newTrialStructure(trialNum).saccadeTime = NaN;
                newTrialStructure(trialNum).reactionTime = NaN;
                newTrialStructure(trialNum).saccadeAmplitude = NaN;
                newTrialStructure(trialNum).saccadeAmplitudeY = NaN;
                newTrialStructure(trialNum).saccadeAmplitudeX = NaN;
                newTrialStructure(trialNum).saccadeDuration = NaN;
                newTrialStructure(trialNum).saccadePeakVel = NaN;
                newTrialStructure(trialNum).saccadeOriginX = NaN;
                newTrialStructure(trialNum).saccadeOriginY = NaN;
                newTrialStructure(trialNum).saccadeLandmarks = nan(1,6);
                newTrialStructure(trialNum).saccadeTailDuration = NaN;
                newTrialStructure(trialNum).saccadeTailAmpltiude = NaN;
                newTrialStructure(trialNum).saccadeAngle = NaN;
            end
            
        end
        
         %need to store the information about saccade the monkey makes to get into the fixation window 
       [val thisTrialFixationSaccadeNumber] = min(abs(allSaccStarts-thisTrialStart));
       relFixSacTime = allSaccStarts(thisTrialFixationSaccadeNumber)-thisTrialStart;
       if relFixSacTime<0 && relFixSacTime>-0.75 %only if this saccade was within 500ms before fixation time
          newTrialStructure(trialNum).fixSacAmplitude = saccadeDetectOut.saccades.amplitude(thisTrialFixationSaccadeNumber,1);
          newTrialStructure(trialNum).fixSacAmplitudeX = saccadeDetectOut.saccades.amplitude(thisTrialFixationSaccadeNumber,2);
          newTrialStructure(trialNum).fixSacAmplitudeY = saccadeDetectOut.saccades.amplitude(thisTrialFixationSaccadeNumber,3);
          newTrialStructure(trialNum).fixSacDuration = saccadeDetectOut.saccades.duration(thisTrialFixationSaccadeNumber);
          newTrialStructure(trialNum).fixSacPeakVel = saccadeDetectOut.saccades.peakV.v(thisTrialFixationSaccadeNumber);
          newTrialStructure(trialNum).fixSacOriginX = saccadeDetectOut.saccades.origin(thisTrialFixationSaccadeNumber,1);
          newTrialStructure(trialNum).fixSacOriginY = saccadeDetectOut.saccades.origin(thisTrialFixationSaccadeNumber,2);
          newTrialStructure(trialNum).fixSaccadeLandmarks = saccadeDetectOut.saccades.landmarks.v(thisTrialFixationSaccadeNumber,:);
        newTrialStructure(trialNum).fixSacStartTimeRelToTrialStart = relFixSacTime;
         %TODO also add end point angle
           xEnd =   newTrialStructure(trialNum).fixSacAmplitudeX+ newTrialStructure(trialNum).fixSacOriginX;
                    yEnd =  newTrialStructure(trialNum).fixSacAmplitudeY+ newTrialStructure(trialNum).fixSacOriginY;
               
         newTrialStructure(trialNum).fixSacAngle  = atan2(yEnd,xEnd);
       else
           
           newTrialStructure(trialNum).fixSacAmplitude = NaN;
          newTrialStructure(trialNum).fixSacAmplitudeX = NaN;
          newTrialStructure(trialNum).fixSacAmplitudeY = NaN;
          newTrialStructure(trialNum).fixSacDuration = NaN;
          newTrialStructure(trialNum).fixSacPeakVel = NaN;
          newTrialStructure(trialNum).fixSacOriginX = NaN;
          newTrialStructure(trialNum).fixSacOriginY = NaN;
            newTrialStructure(trialNum).fixSaccadeLandmarks = nan(1,6);
           newTrialStructure(trialNum).fixSacStartTimeRelToTrialStart = NaN;
           
           newTrialStructure(trialNum).fixSacAngle = NaN;
       end
    end
    
    
   
end