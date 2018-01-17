function [calData] = calibrateEyeData2016(thisFileInfo,eyeX,eyeY,eyeFs,posThresh,trialData,conditionLegend,calibFileNameEnd)
%function that calibrates monkey eye data and cleans it for saccades and
%saturation points
%either loads previosuly done calibration or does it and saves it.
%INPUTS
%conditionLegend.test contrains a cell array giving the condition codes ad positions

%TODO for more tahn 16 codntions calibration should extend

fileNameBase = thisFileInfo.fileName(1:end-4);
numConditions = size(conditionLegend.test,1);

%check if calibration already done and either load or do
%%calibration of eye trace using strucutres in cell above
%TODO run cleaning first???
if exist([fileNameBase calibFileNameEnd],'file')==2
    load([fileNameBase calibFileNameEnd],'xOffset','yOffset','xMagRatio','yMagRatio')
else
    %TODO finish copying from other script ericElectro.. line 588
    numTrials = size(trialData,2);
     %conditionLegend = load('condtionLegend.mat'); %contrains cell with target locations
    cc = lines(12);
    
   
    
    figure; %for displaying stat and end of saccades
    hAx1 = axes;
    hAx2 = axes('XAxisLocation','top','YAxisLocation','right','color','none',...
        'Xcolor',[0 0.5 0],'Ycolor',[0 0.5 0],'XLim',[-15 15],'YLim',[-15 15]);
    calMatrix = zeros(6,numTrials); %TODO make nan so ignores bad trials?
    
    for trialNum = 1:numTrials
        %get the time of bit 8 and bit 3 and then get the eye trace
        %+50:+100 samples and -200:-100 samples respectively
         bit8Ts = round(trialData(trialNum).bit8time*eyeFs)
        if ~isempty(bit8Ts)
           
            calMatrix(1,trialNum) = mean(eyeX(bit8Ts+100:bit8Ts+200));
            calMatrix(2,trialNum) = mean(eyeY(bit8Ts+100:bit8Ts+200));
            calMatrix(3,trialNum) = trialData(trialNum).conditionCode;
            bit3Ts = round(trialData(trialNum).bit3time*eyeFs)
            calMatrix(4,trialNum) = mean(eyeX(bit3Ts-200:bit3Ts-100));
            calMatrix(5,trialNum) = mean(eyeY(bit3Ts-200:bit3Ts-100));
            
            %draw a line represent start and end and label with condition code
            text(calMatrix(1,trialNum),calMatrix(2,trialNum),num2str(calMatrix(3,trialNum)),...
                'Color',[0.5 0 0],'Parent',hAx1);
            line([calMatrix(4,trialNum) calMatrix(1,trialNum)],...
                [calMatrix(5,trialNum) calMatrix(2,trialNum)],'Parent',hAx1);
            
            
            %now on the other axes draw the full saccade path
            samplesAroundSaccade = [200 500]; %[prebit3 postbit8]
            saccadePath = [eyeX(bit3Ts-samplesAroundSaccade(1):bit8Ts+samplesAroundSaccade(2)); eyeY(bit3Ts-samplesAroundSaccade(1):bit8Ts+samplesAroundSaccade(2))];
            
            %   line(saccadePath(1,:),saccadePath(2,:),'parent',hA(1))
            %    x = 0:.05:2*pi;
            %y = sin(x);
%             hF2 = figure; %for displaying full saccade paths
%             hA(1) = axes('parent',hF2);
%             z = zeros(1,size(saccadePath,2));
%             col = linspace(0,2*pi,size(saccadePath,2));  % This is the color, vary with x in this case.
%             surface([saccadePath(1,:);saccadePath(1,:)],[saccadePath(2,:);saccadePath(2,:)],[z;z],[col;col],...
%                 'facecol','no',...
%                 'edgecol','interp',...
%                 'linew',2,'parent',hA(1));
%             %now add markers for bit 3 and bit 8
%             line([nan saccadePath(1,samplesAroundSaccade(1))],[nan saccadePath(2,samplesAroundSaccade(1))],'marker','s','parent',hA(1),'markerfacecolor','k','markersize',5,'color','k')
%             line([nan saccadePath(1,size(saccadePath,2)-samplesAroundSaccade(2))],[nan saccadePath(2,size(saccadePath,2)-samplesAroundSaccade(2))],'marker','^','parent',hA(1),'markerfacecolor','k','markersize',5,'color','k')
%             
%             %now mark current guesss of finishing and starting point from calmatrix
%             line([nan calMatrix(1,trialNum)],[nan calMatrix(2,trialNum)],'marker','^','parent',hA(1),'markerfacecolor','r','markersize',5,'color','r')
%             line([nan calMatrix(4,trialNum)],[nan calMatrix(5,trialNum)],'marker','s','parent',hA(1),'markerfacecolor','r','markersize',5,'color','r')
        else
            calMatrix(:,trialNum) = nan;
            calMatrix(3,trialNum) = trialData(trialNum).conditionCode;
        end
    end
    
    
     %loop through conditions and plot target on second axis
     
    for conditionNum =1:numConditions
        rectangle('position',[conditionLegend.test{conditionNum,2} conditionLegend.test{conditionNum,3},...
            1 1],'Parent',hAx2,'FaceColor',cc(conditionLegend.test{conditionNum,4},:));
        text(conditionLegend.test{conditionNum,2}+0.5,conditionLegend.test{conditionNum,3}+0.5,...
            num2str(conditionLegend.test{conditionNum,1}),'Parent',hAx2,'Color',[1 1 1])
        
    end
    
    corTrials  = find([trialData.correctResponse]==2);
    %loop over correct trials
   numCorTrials = length(corTrials)
    H = zeros(1,numCorTrials);
    V = H;
    x = H;
    y = H;
    for trialNum = 1:numCorTrials
        x(1,trialNum) = calMatrix(4,corTrials(trialNum));
        y(1,trialNum) = calMatrix(5,corTrials(trialNum));
        H(1,trialNum) = conditionLegend.test{calMatrix(3,corTrials(trialNum)),2};
        V(1,trialNum) = conditionLegend.test{calMatrix(3,corTrials(trialNum)),3};
     
    end
    
    %make 0,0 by minus mean of fixation
    calMus = nanmean(calMatrix,2);
%     fixX  = nanmean(calMatrix(4,:)); %this gave an error of averaging
%     over the worng dimension
%     fixY  = nanmean(calMatrix(5,:));
    fixX  = calMus(4);
    fixY  = calMus(5);
    calMatrix([1 4],:) = calMatrix([1 4],:)-fixX;
    calMatrix([2 5],:) = calMatrix([2 5],:)-fixY;
    %replot with new 0,0
    
    figure;
    hAxR1 = axes;
    
    hAxR2 = axes('XAxisLocation','top','YAxisLocation','right','color','none',...
        'Xcolor',[0 0.5 0],'Ycolor',[0 0.5 0],'XLim',[-15 15],'YLim',[-15 15]);
    
    for trialNum = 1:numTrials
        text(calMatrix(1,trialNum),calMatrix(2,trialNum),num2str(calMatrix(3,trialNum)),...
            'Color',[0.5 0 0],'Parent',hAxR1);
        line([calMatrix(4,trialNum) calMatrix(1,trialNum)],...
            [calMatrix(5,trialNum) calMatrix(2,trialNum)],'Parent',hAxR1);
    end
    for conditionNum =1: numConditions
        rectangle('position',[conditionLegend.test{conditionNum,2} conditionLegend.test{conditionNum,3},...
            1 1],'Parent',hAxR2,'FaceColor',cc(conditionLegend.test{conditionNum,4},:));
        text(conditionLegend.test{conditionNum,2}+0.5,conditionLegend.test{conditionNum,3}+0.5,...
            num2str(conditionLegend.test{conditionNum,1}),'Parent',hAxR2,'Color',[1 1 1])
        
    end
    
    
    
    
    %scale by dividing endpoints
    x = [];
    y = [];
    trialCount = 1;
    for trialNum = corTrials
        x(trialCount) = calMatrix(1,trialNum);
        y(trialCount) = calMatrix(2,trialNum);
        trialCount = trialCount +1;
    end
    xRatios = bsxfun(@rdivide,x,H);
    xRatios(abs(xRatios)==Inf) = NaN;
    yRatios = bsxfun(@rdivide,y,V);
    yRatios(abs(yRatios)==Inf) = NaN;
    calMatrix([1 4],:) = calMatrix([1 4],:)./nanmean(abs(xRatios),2); %PH
    calMatrix([2 5],:) = calMatrix([2 5],:)./nanmean(abs(yRatios),2); %PH
    
    
    
     figure;
    hAxR1 = axes;
    %xlim([-2000 600])
    %ylim([-2000 1200])
    hAxR2 = axes('XAxisLocation','top','YAxisLocation','right','color','none',...
        'Xcolor',[0 0.5 0],'Ycolor',[0 0.5 0],'XLim',[-15 15],'YLim',[-15 15]);

    for trialNum = 1:numTrials
        text(calMatrix(1,trialNum),calMatrix(2,trialNum),num2str(calMatrix(3,trialNum)),...
            'Color',[0.5 0 0],'Parent',hAxR1);
        line([calMatrix(4,trialNum) calMatrix(1,trialNum)],...
            [calMatrix(5,trialNum) calMatrix(2,trialNum)],'Parent',hAxR1);
    end
    for conditionNum = 1:numConditions
        rectangle('position',[conditionLegend.test{conditionNum,2} conditionLegend.test{conditionNum,3},...
            1 1],'Parent',hAxR2,'FaceColor',cc(conditionLegend.test{conditionNum,4},:));
        text(conditionLegend.test{conditionNum,2}+0.5,conditionLegend.test{conditionNum,3}+0.5,...
            num2str(conditionLegend.test{conditionNum,1}),'Parent',hAxR2,'Color',[1 1 1])
        
    end
    linkaxes([hAxR2 hAxR1])
    
     xOffset = fixX;
    yOffset = fixY;
    xMagRatio = nanmean(abs(xRatios),2);
    yMagRatio = nanmean(abs(yRatios),2);
    
    
    
      save([fileNameBase calibFileNameEnd],'xOffset','yOffset','xMagRatio','yMagRatio')
    %return
end



%now use calibration values to calibrate eye movement signals and cahnge on
%amster plot
%-mean and divede by ratio for each direction
eyeX = (eyeX- xOffset)./xMagRatio;
eyeY = (eyeY- yOffset)./yMagRatio;


%now clean saccades by finding time points in which position
%exceeds a thrshold


eyeX2 = eyeX;
noiseSegments = find(abs(eyeX2)>posThresh);
%remove sections of noise near beginning and end of file
noiseSegments(noiseSegments<4) = [];
noiseSegments(noiseSegments>(length(eyeX2)-3)) =[];

%replace saccades with NaNs then interpolate over NaN sections
eyeX2(noiseSegments) = NaN;
eyeX2(noiseSegments-3) = NaN;
eyeX2(noiseSegments+3) = NaN;
eyeX2(1) = 0; %zeros first and last to prevent NaNs in interpolation
eyeX2(end) = 0;
bd = isnan(eyeX2);
gd = find(~bd);
bd([1:(min(gd)-1) (max(gd)+1):end])=0;
eyeX2(bd) = interp1(gd,eyeX2(gd),find(bd));

%do the smae for eyeY %TODO these shoudl work interactively
noiseSegments = find(abs(eyeY)>posThresh);
noiseSegments(noiseSegments<4) = [];
noiseSegments(noiseSegments>(length(eyeY)-3)) =[];
eyeY2 = eyeY;

eyeY2(noiseSegments) = NaN;
eyeY2(noiseSegments-3) = NaN;
eyeY2(noiseSegments+3) = NaN;
eyeY2(1) = 0; %zeros first and last to prevent NaNs in interpolation
eyeY2(end) = 0;
bd = isnan(eyeY2);
gd = find(~bd);
bd([1:(min(gd)-1) (max(gd)+1):end])=0;
eyeY2(bd) = interp1(gd,eyeY2(gd),find(bd));

% now downsample to 1K
eyeXds = resample(eyeX2,1000,round(eyeFs));
eyeYds = resample(eyeY2,1000,round(eyeFs));

%simple filter
windowWidth = 5; %TODO move to input arg
eyeXf = filtfilt(ones(1,windowWidth)/windowWidth,1,eyeXds);
eyeYf = filtfilt(ones(1,windowWidth)/windowWidth,1,eyeYds);

calData.eyeXcal = eyeXf;
calData.eyeYcal = eyeYf;
calData.eyeFsCal = 1000;
calData.eyeTcal = linspace(0,length(eyeXf)/1000,length(eyeXf));
