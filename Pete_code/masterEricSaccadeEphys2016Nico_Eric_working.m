%masterEricSaccadeElectrophys - just for testing nicos new files

%master script to do analysis of all Eric' saccade monkey data
%this is a trimmed down version of masterEricSaccadeEphys

%% Main VARIALBE LIST
% analysisResults(fileNum).newTrialStructure - output from
%   appendRasterToTrials includes alignedSpikes, event info and saccade
%   detection info with everything split into trials and aligned to their
%   onset

%% PATHS

% Different pcs -  different paths
set(0,'defaultfigurecolor','w')
set(0,'defaultaxescolor','w')


[~, cn] = system('hostname');
switch strtrim(cn)
    case 'NIN673'  % analysis pc at nin 

cellList = 'full';
addpath(genpath('D:\CODE\AntiSaccadeCode\_github'));

externalHdDriveLetter = 'D';
rawDataDir = [externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\_to_run'] ;%'I:\Raw Data\Pro v Anti Saccade Task';
processDir = [externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\_process']; %for holding intermediate files, -dataStructure, -trialInfo, -saccDetect, -calib

resortDir = [externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\_to_run']; %loaction of resorted spike files, just neesd adding to path
neuronListDir =  [externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\neuronLists'];



addpath(genpath([externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\_to_run\Mickey Lateral\_spikeSorted']))
addpath(genpath([externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\_to_run\Mickey Vermis\_spikeSorted']))
addpath(genpath([externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\_to_run\Moshe Lateral\_spikeSorted']))
addpath(genpath([externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\_to_run\Moshe Vermis\_spikeSorted']))


addpath(genpath(resortDir))
 addpath(genpath(rawDataDir))
  
 newList =1;
 if newList
 load([neuronListDir '\neuronList_complete5.mat'])
 end
%directory where all results will be saved
if strcmp(cellList,'full')
resultsDir = [externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\_results']; %for full list
elseif strcmp(cellList,'clean')
resultsDir = [externalHdDriveLetter ':\DATA\monkey\Pro_and_Anti\_PeteAnalysis\data\results_clean']; %for clean units only
end
%list

%create results directory if not already there
if exist(resultsDir,'dir')~=7
    mkdir(resultsDir)
end
if exist(processDir,'dir')~=7
    mkdir(processDir)
end

%this section makes sure mario's asccade detection fucntions are on the
% %path
% diskPath = 'D:\Nico\_PeteAnalysis\pro_anti-git\Pete's_code\Code\Mario_saccade_detection';
% codeDir = ['D:\Nico\_PeteAnalysis\AnalysisCodePete'];
% addpath(genpath(codeDir))
% rmpath([codeDir '\new mario'])
% rmpath([codeDir '\old mario'])
% rmpath(genpath([diskPath '\FEX\chronux']))
% addpath([diskPath '\FEX\Info Toolbox']) %for MI
%  addpath([diskPath '\FEX\outlier_library_external'])
%  addpath([diskPath '\FEX\mi'])
%  addpath([diskPath '\FEX\distributionPlot\distributionPlot'])
%  addpath([diskPath '\FEX\xticklabel_rotate'])
%  addpath([diskPath '\FEX\countmember'])
%   addpath([diskPath '\FEX\SubAxis'])
  
  
  
    case '...' % Here goes you pc-name ( to find pc name excute: system('hostname')

        
externalHdDriveLetter = 'C';
rawDataDir = [externalHdDriveLetter ':\Users\erico\Box Sync\NIN\data\N1'] ;%'I:\Raw Data\Pro v Anti Saccade Task';
processDir = [externalHdDriveLetter ':\Users\erico\Box Sync\NIN\data']; %for holding intermediate files, -dataStructure, -trialInfo, -saccDetect, -calib

resortDir = [externalHdDriveLetter ':\Users\erico\Box Sync\NIN\data']; %loaction of resorted spike files, just neesd adding to path
addpath(genpath(resortDir))
%addpath(genpath(rawDataDir))
%directory where all results will be saved
if strcmp(cellList,'full')
resultsDir = [externalHdDriveLetter ':\Users\erico\Box Sync\NIN\Results']; %for full list
elseif strcmp(cellList,'clean')
resultsDir = [externalHdDriveLetter ':\Users\erico\Box Sync\NIN\Results\results_clean']; %for clean units only
end
%list

%create results directory if not already there
if exist(resultsDir,'dir')~=7
    mkdir(resultsDir)
end
if exist(processDir,'dir')~=7
    mkdir(processDir)
end

%this section makes sure mario's asccade detection fucntions are on the
%path

diskPath = 'C:\Users\erico\Box Sync\NIN\Code';
codeDir = [diskPath '\Mario_saccade_detection'];
addpath(genpath(codeDir))

    otherwise
        disp('Fill in your PC-name at line at :  case(...)')
        return  
end

% userName = getenv('USERNAME'); %Makes it work on PC and Mac 
% %codeDir = ['/Users/' userName '/Backup Erasmus/NIN/Matfiles & Analysis/Analysis/Active Plot Generator'];
% codeDir = ['C:\Users\' userName '\Google Drive\Current Matlab Code\EricSaccadeElectrophys'];
% addpath(genpath(codeDir))
% addpath(['C:\Users\' userName '\Google Drive\Current Matlab Code\FEX\mi']) %toolbox for MIcalculation
% rmpath(['C:\Users\' userName '\Google Drive\Current Matlab Code\EricSaccadeElectrophys\new mario'])
% rmpath(['C:\Users\' userName '\Google Drive\Current Matlab Code\EricSaccadeElectrophys\old mario'])
% addpath(genpath(['C:\Users\' userName '\CloudStation\Current Analysis\FEX\CircStat2012a']))
%% RECORDING LIST
%list of all neurons is contained within the following .m file, running it
%as a script loads the list as a struct

%recordingListEricPvAElectrophys %final list
%recordingListEricPvAElectrophysClean %final list of only clean cs pause neurons
 
if strcmp(cellList,'full') && ~newList
recordingListEricPvAElectrophysNicoRec%recordingListEricPvAElectrophysLocResort %recordingListEricPvAElectrophys %for full list
elseif strcmp(cellList,'clean') && ~newList
recordingListEricPvAElectrophysClean %for clean units only %TODO won't work any more as no resortFLag field
end
%find recordings to analyse by looking for non-empty recording names
  
   
recList = squeeze(struct2cell(neuronList));
 
if newList
recList2= recList;

recList(1,:) = recList2(2,:);% name
recList(2,:) = recList2(7,:);% fileList
recList(3,:) = [];
recList(4,:) = recList2(3,:);% area
recList(5,:) = []; %depth
recList(6,:) = [];% gridloc
recList(7,:) = recList2(6,:);
recList(8,:) = recList2(1,:);% path
recList(9,:) = recList2(4,:);% monkey
recList(10,:) = recList2(5,:); %datetime
recList(11,:) = recList2(8,:); % stimListPath
% recList(12,:) = 
 
recordingsToAnalyse = find(cellfun(@(x) ~isempty(x),recList(1,:)));
end
exampleNeurons = [14 31]; %neurons for which eps of the rasters will be exported
%find where which are neurons are recorded from using field in recording list
%recList(4,find(cellfun(@(x) isempty(x),recList(4,:))))={nan};
% vermNeurons  = find(strcmpi(recList(4,recordingsToAnalyse),'vermis'));
% latNeurons = find(strcmpi(recList(4,recordingsToAnalyse),'lateral')); 
vermNeurons  = find(strcmpi(recList(4,:),'vermis'));
latNeurons = find(strcmpi(recList(4,:),'lateral')); 

%recordingsToAnalyse = 1; %TODO remove after testing
%% SAVED FILE NAMES
if strcmp(cellList,'full')
dataFileNameEnd = '-dataStructureResort.mat';
elseif strcmp(cellList,'clean')
dataFileNameEnd = '-dataStructureClean.mat'; %data for all files in one neuron saved into one file
end

%
mainPlotFlag = 1; %flag to control plotting of main plot with data, eye and trial info
analyseISI = 1; %flag to control plotting of ISIs for each unit
checkCSpikePause = 1; %flag to control the plotting and saving of a CS triggered SS raster
figPdfFileNameEnd = '-outputPdfFig2016.pdf'; %pdf file containing results from indiviudla files
analysisFileNameEnd = '-anFile2016fixsacc2016v2-wholeNeuronCalib.mat';  %'-anFile2015.mat';%data for saving analysis to
wholeNeuronNameEnd = 'wholeNeuronResultsResort2016realignfixsacc2016-fixedbit8.mat'; %store the final trial structures after removing unwanted unstabel areas. Make sure to overwrite or create new
% allNeuronSumFigName = [resultsDir filesep 'ericAllNeuronsFinalNewTimesB0508.pdf'];
% indNeuronSumFigEnd = '-newStyleSummaryBurstFullRunNewTimesB0508.pdf';
% burstSumFigName = [resultsDir filesep 'burstSummaryNewTimes1008test3.pdf'];
burstSumPlotFlag = 0; %TODO may have to go after the correlations

%population response 
doPop = 1; %flag to turn on and off
popResPdfName = [resultsDir filesep 'ericPopulationResults2016moreSacs50binAngle350to450vel-csonlabelled-selNeurons.pdf'];
popResFileName = 'populationResults2016moreSacs50bin-selNeurons.mat'; %saved results, if adding more neurons make sure to overwrite
doPopCV2 = 0;
cvPopResPdfName = [resultsDir filesep 'cvPopulationResults2016.pdf'];

bitNum = 8; %bit to trigger from in eventData(6,:) 
secondBitNum = 3; %second trigger to mark on rasters eventData(7,:)


%burst saccade correlation
corrPdfName = 'ericCorTestPArealign2015Resort.pdf';
doCor = 0;

%MI
doMI = 0;
miSumFigName = 'ericMiTest2016.pdf';
doPopMi = 0;
doSlidWinMI = 0;
doSlidWinMIinst = 0;
burstFileNameEnd = '-burstDetectionFile2015Resort.mat';
burstFigPdfFileNameEnd = '-burstDetectionTes2015Resort.pdf';
burstPlotFlag = 0;

%kunimatsu filenames
indNeuronSumFigEnd = '-kuniResults2017-selectedNeurons.pdf';
allNeuronSumFigName = [resultsDir filesep 'kuniResults201realign2017-selectedNeurons.pdf'];
allKuniResFile = [resultsDir filesep 'kunimatsu2017-selectedNeurons.mat'];

%max/min
doMaxMin = 1; %only applies to pro v anti not directional
maxMinFigNameEnd = 'maxminResort2016-selNeurons.pdf';
maxMinByDirPdfName = 'maxminspider2016-selNeurons.pdf';
doMaxMinByDir = 1;

doRasters = 1; %saccade rasters
doInstRasters = 1; %isntruction period rasters
doCsByDir = 1; %counts of Cs in each direction
allCsPausePdfName = [resultsDir filesep 'csPause2016-selNeurons.pdf']; %cs pause for whole neuron

doSpikesInBins = 1; %spike counting indifferent bins
allNeuronPdfName = [resultsDir filesep 'alignToMaxProCSsummaries2016-selNeurons.pdf']; %filename for above

%code sand times
trialConditionCodes = [14 6;15 7;8 16;4 12;NaN NaN;10 2;11 3;5 13;1 9]; %this is for telling which trials are pro and anti and which go together, in subplot order
timeWindow = [4 4]; %time window (s) for extraction of each trial (relative to trial start)
baseInd = [3000 3250]; %baseline here will be taken as ?[-0.4 0] to trialStart in kunimatsu so
%[-1.25 -0.75] relative to saccade. 

csStaFigNameEnd = 'csSta-test2016.pdf';
wholeNeurCsSta = 0;
wholeNeurCsPause = 0;

doEndPoints = 1;
erAngPdfName = [resultsDir filesep 'errorAngle2016-FixEp-selNeurons.pdf'];
endPointErrorsFileName = 'csOnDirsFixEp.mat';
doEndPointsFix = 0;
endPointErrorsFixationFileName = 'csOnDirsFixation.mat';
erAngFixationPdfName = [resultsDir filesep 'errorAngleFixation2016-Fix.pdf'];


doCorMatrix = 0;
doPopCorMatrix = 0;
doCorKin = 0;
corKinFigNameStart = 'kinematicsCorrelationstest';

csTimingPdfName = 'timingOfCs-20msGaus.pdf';


rasterByDirPdfName = 'by direction rasters.pdf';


doPopResByDir = 0;
popByDirFigName = [resultsDir filesep 'popByDirFixEp-comp.pdf'];
popByDirFigNameRot = [resultsDir filesep  'popByDirFixEpRotated-comp-fix2.pdf'];
%% CHRIS SUMMARY FIGURE PNG NAMES
muSdfPngNameEnds = {'-muSdfSS2016.png','-muSdfCS2016.png'}; %for each unit
csPauseFigNameEnd = '-csPauseTest2016.png';

ssDirTuningPngNameEnd = '-SStuning2016.png'; %also has a prefix of proNadir,proZenith,proLowTime,proHighTime,proHigh or proLow eg name-prefix-ssDirTuningPngName

csDirTuningPngNameEnd = '-CStuning2016.png';

instRastPngNameEnd = '-instructionRast2016.png'; %also has prefix

unsortedRasterPngName = '-trialOrderRaster2016.png'; %also has prefix
%TODO CS algined to instruction, SS tuned to instruction
saccadeSdfPngNameEnds = {'Saccade-muSdfSS2016.png','Saccade-muSdfCS2016.png'}; %for each unit
instructSdfPngNameEnds = {'Instruct-muSdfSS2016.png','Instruct-muSdfCS2016.png'}; %for each unit
%% COLOURs AND LABELS

cc = lines(12); %default coulourmap for figures

rastToPlot = {'corSacTrials','corProTrials','corAntiTrials'};
byDirLabel = {1,'Pro';2,'Anti'};
ccProAnti = [0 1 0;1 0 0];

%% DATA EXTRACTION
%TODO eyeDelay = 0;%0.013; %13ms hardware delay
%loop over strucutre and test for non converted files
%recordingsToAnalyse = recordingsToAnalyse(recordingsToAnalyse>=29);
 
% % recordingsToAnalyse = [27 28 29 38 40 42];
 
notEmptyNeurons = find(~cellfun(@isempty, recList(4,:)));

% recordingsToAnalyse = [72] ;

% recordingsToAnalyse =  notEmptyNeurons(1:20);
% recordingsToAnalyse =  notEmptyNeurons(21:40);
% recordingsToAnalyse =  notEmptyNeurons(41:50);
% recordingsToAnalyse =  notEmptyNeurons(51:60);
% recordingsToAnalyse =  notEmptyNeurons(61:70);
% recordingsToAnalyse =  notEmptyNeurons(81:90);
% recordingsToAnalyse =  notEmptyNeurons(91:100);
% recordingsToAnalyse =  notEmptyNeurons(101:110);
recordingsToAnalyse =  notEmptyNeurons(111:120);

% flag 667

if ~exist('error_log.mat')
    error_log = cell(size(neuronList,2),3);
    for i= 1:size(neuronList,2)
        error_log{i,1} =  neuronList(i).neuronName;
        error_log{i,2} = neuronList(i).path;
    end
    
    save([neuronListDir '/error_log.mat'])
    
else 
    load([neuronListDir '/error_log.mat'])
end

    
    
for neuronNum = recordingsToAnalyse;
    try
    if exist([processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd],'file')~=2
          
        disp(['stability marking for ' neuronList(neuronNum).neuronName])
        %convert file list to structure
        stimDescript{neuronNum} = cell2struct(neuronList(neuronNum).fileList, ...
                    {'fileName', 'spkChannelNumber', 'unitNumbers', 'startTime', 'stopTime','stimStartFile','resortFlag'}, ...
                    2);
        %launch stable period marking
        %PH 171016 added switch to enable use of different eye channels
        if isempty(neuronList(neuronNum).eyeChannels)
        hFig = markStablePeriods(stimDescript{neuronNum},[17 18]);
        drawnow;
        waitfor(hFig)
        else
            
            hFig = markStablePeriods(stimDescript{neuronNum},neuronList(neuronNum).eyeChannels);
        drawnow;
        waitfor(hFig)
            
        end

        neuronList(neuronNum).fileList = editedFileList;
        save([processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd],'dataStructure','editedFileList','-v7.3')
        clear dataStructure editedFileList
    else
        disp(['stability already marked for ' neuronList(neuronNum).neuronName])
    
    end
    
    catch em
        
        error_log{neuronNum,3} = em ;
            save([neuronListDir '/error_log.mat'],'error_log')

    end 
    
end

 keyboard
%% FILE BY FILE ANALYSIS
%loop over all neurons and files within those neurons
for neuronNum = recordingsToAnalyse;
    close all %TODO be careful!!
    disp(['    processing file ' num2str(neuronNum) ' of ' num2str(max(recordingsToAnalyse))])
    
    %create file name for the summary figure pdf for this neuron
    figureFileName = [resultsDir filesep neuronList(neuronNum).neuronName figPdfFileNameEnd];
    
     
    %load matfile into memory
     thisMatFileName = [processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd];
     mLoadedData = matfile(thisMatFileName); %create mat file
    
   
     
     
     %check for each neuron if the analysis file has already been saved, if
     %not then analyse
     if exist([resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd],'file')~=2
         
         
         numFiles = size(mLoadedData.editedFileList,1);
         %create blank structure for holdind analysis results for this neuron
        structInitializer = cell(numFiles,1); %has 2 comlums as second is for the sceond unit (complex spikes)
         analysisResults = struct('saccadeDetection',structInitializer,'newTrialStructure',structInitializer,'alignedSpikeTimesCell',structInitializer);
         %loop over all files and load the required data
         
         
        
         
         for fileNum = 1:numFiles
 
             disp(['    processing file ' num2str(fileNum) ' of ' num2str(numFiles)])
             
             %load data from mat file
             thisData = mLoadedData.dataStructure(fileNum,1);
             thisFileInfo = mLoadedData.editedFileList(fileNum,1);
             thisFileInfo.neuronNum = neuronNum;
             thisFileContents = whos(mLoadedData);
             
              
     fileNameBase = thisFileInfo.fileName(1:end-8);

%      foldercontent = getAllFiles([resortDir filesep neuronList(neuronNum).neuronName]);
     
%      foldercontent = getAllFiles([neuronList(neuronNum).neuronName]);
%           foldercontent = getAllFiles([neuronList(neuronNum).path]);
% %           load(cell2mat(foldercontent( cellfun(@isempty, regexp(foldercontent,'\._'))))) % load folder contents

      
%      load(cell2mat(foldercontent( ~cellfun(@isempty, regexp(foldercontent,'Moshe'))))) % load folder contents
% load(cell2mat(foldercontent))

% xx=vertcat(foldercontent{cellfun(@isempty, regexp(foldercontent,'\._')),:})
    load([neuronList(neuronNum).stimListFile])

%     for i = 1: 312
%         if ~isempty( neuronList(i).path)
%             prefix = 'D:\DATA\monkey\Pro_and_Anti';
%             neuronList(i).path = [prefix neuronList(i).path(8:end)];
%         end
%         if ~isempty( neuronList(i).stimListFile)
%             prefix = 'D:\DATA\monkey\Pro_and_Anti';
%             neuronList(i).stimListFile = [prefix neuronList(i).stimListFile(8:end)]
%         end
%     end
%     
            
            
     stimfilename = whos('-file',neuronList(neuronNum).stimListFile);
     load(neuronList(neuronNum).stimListFile)
     
     if ~strcmp(stimfilename,'stimListFile')
         eval(['stimListFile ='  stimfilename.name]);
     end
     
       
%      save([processDir filesep neuronList(neuronNum).neuronName '-stimList.mat'], 'stimListFile')
%     stimListFile=denoise_eye(stimListFile);

              save([processDir filesep  thisFileInfo.fileName(1:end-9) '-stimList.mat'], 'stimListFile')
         
         
             %get trial information
             thisTrialInfoName = [processDir filesep thisFileInfo.fileName(1:end-4) '-trialInfo.mat'];
             
             if exist(thisTrialInfoName,'file')~=2
                  
                 trialData = extractTrialData(thisFileInfo,thisData.spkDataStartTime,bitNum,secondBitNum,processDir,neuronList);
                 save(thisTrialInfoName,'trialData');
             else
                 load(thisTrialInfoName);
             end
             
             
             
             wholeNeuronCalibFileNameEnd = '-wholeNeurCal.mat';
             %below is the old calibration method that get's a different
             %cal for each file, instead we should check that all of the
             %files have a calib file and then take the mean of all as the
             %final calib measure
             conditionLegend = load('conditionLegend.mat');
             %calibrate eye movements or load from structure
            posThresh =12; % 12; %positon threshold for saccades
                         
            
             
            calibFileNameEnd = '-calib2016.mat';
            if fileNum==1
               
                allFilesCalib = cell(numFiles,1);
                for  tempFileNum = 1:numFiles
                    thisTempFileInfo = mLoadedData.editedFileList(tempFileNum,1);
                    fileNameBase = thisTempFileInfo.fileName(1:end-4);
                    if exist([fileNameBase calibFileNameEnd],'file')==2  %&& do_sac_detect==0
                        tempCal =  load([fileNameBase calibFileNameEnd],'xOffset','yOffset','xMagRatio','yMagRatio')
                        allFilesCalib{tempFileNum} = tempCal;
                    else
                          [calData] = calibrateEyeData2016(thisTempFileInfo,thisData.eyeX,thisData.eyeY,thisData.eyeFs,posThresh,trialData,conditionLegend,calibFileNameEnd);
                        tempCal = load([fileNameBase calibFileNameEnd],'xOffset','yOffset','xMagRatio','yMagRatio')
                        allFilesCalib{tempFileNum} = tempCal;
                        %don't we need to save?
                    end
                end
                
                
             wholeNeuronCalibTemp = [allFilesCalib{:}];
             wholeNeuronCalib.xOffset = mean([wholeNeuronCalibTemp.xOffset]);
              wholeNeuronCalib.yOffset = mean([wholeNeuronCalibTemp.yOffset]);
               wholeNeuronCalib.xMagRatio = mean([wholeNeuronCalibTemp.xMagRatio]);
                wholeNeuronCalib.yMagRatio = mean([wholeNeuronCalibTemp.yMagRatio]);
                
                %TODO now save the wholeNeuron calibration
                
                save([processDir filesep neuronList(neuronNum).neuronName wholeNeuronCalibFileNameEnd],'wholeNeuronCalib')
            else
                %TODO now load the wholeNeuronCalib
                 load([processDir filesep neuronList(neuronNum).neuronName wholeNeuronCalibFileNameEnd],'wholeNeuronCalib')
            end
            
            
            %now we want to calibrate the fileusing the those values
            [calData] = simpleEyeCalibration(wholeNeuronCalib,thisData.eyeX,thisData.eyeY,thisData.eyeFs);
            
            %this is temporary editing to force recalibration without
            %redoing dataStructure PH 220316
          %  [calData] = calibrateEyeData2016(thisFileInfo,thisData.eyeX,thisData.eyeY,thisData.eyeFs,posThresh,trialData,conditionLegend,calibFileNameEnd);
%                 
                %temporarily enable writing to file and save
                mLoadedData.Properties.Writable = true;
                mLoadedData.calibEye(fileNum,1) = calData;
                mLoadedData.Properties.Writable = false;
%             if sum(strcmp({thisFileContents.name},'calibEye'))<1 %if doesn't exist
%                 %[calData] = calibrateEyeData(thisFileInfo,thisData.eyeX,thisData.eyeY,thisData.eyeFs,posThresh,trialData,conditionLegend);
%                 [calData] = calibrateEyeData2016(thisFileInfo,thisData.eyeX,thisData.eyeY,thisData.eyeFs,posThresh,trialData,conditionLegend,calibFileNameEnd);
%                 
%                 %temporarily enable writing to file and save
%                 mLoadedData.Properties.Writable = true;
%                 mLoadedData.calibEye(fileNum,1) = calData;
%                 mLoadedData.Properties.Writable = false;
%             else
%                 if size(mLoadedData.calibEye,1)~=numFiles %if wrong size
%                    % [calData] = calibrateEyeData(thisFileInfo,thisData.eyeX,thisData.eyeY,thisData.eyeFs,posThreshhresh,trialData,conditionLegend);
%                     [calData] = calibrateEyeData2016(thisFileInfo,thisData.eyeX,thisData.eyeY,thisData.eyeFs,posThresh,trialData,conditionLegend,calibFileNameEnd);
%                
%                     %temporarily enable writing to file and save
%                     mLoadedData.Properties.Writable = true;
%                     mLoadedData.calibEye(fileNum,1) = calData;
%                     mLoadedData.Properties.Writable = false;
%                 else
%                     calData = mLoadedData.calibEye(fileNum,1);
%                 end
%             end
             
             %shift time vector for teh spikes to be 0 at start of file
            alignedSpkTimeVector = thisData.spkTimeVec-thisData.spkDataStartTime;
             %align all spikeTimes to the start of file
            spikeTimesCell = cellfun(@(x) (x-thisData.spkDataStartTime),thisData.alignedSpikes,'uniformoutput',false);
            
            %now also take into account the delay of13ms in the recording
            %of the eyeData %TODO this will also have to be taken into
            %account by any LFP analysis
            spikeTimesCell = cellfun(@(x) (x+0.013),spikeTimesCell,'uniformoutput',false);
            alignedSpkTimeVector = alignedSpkTimeVector+0.013;
            
            %detect and analyse saccades
%             if isempty(analysisResults(fileNum,1).saccadeDetection)
%                 saccadeDetectOut = saccadeDetection([calData.eyeXcal' calData.eyeYcal'],'sgolayall', 1,'sgolayfilter', 1,'sgolaywindow', 20, 'detectionwindow', [120 200]);
%                 analysisResulets(fileNum,1).saccadeDetection = saccadeDetectOut;
%                 allSaccStarts = saccadeDetectOut.saccades.landmarks.v(:,1)/1000; %vector of saccade times for plotting
%             else
%                 saccadeDetectOut = analysisResults(fileNum,1).saccadeDetection;
%             end
            %now with file instead of fields for keeping
             
            saccFileName = [processDir filesep thisFileInfo.fileName(1:end-4) '-saccDetectNewestRecalib-wholeNeurCal.mat'];
            if exist(saccFileName,'file')~=2   
                
%                filter out high freq peaks from poor tracker calibration
             calData.eyeXcal= medfilt1(calData.eyeXcal,40);
             calData.eyeYcal= medfilt1(calData.eyeYcal,40);
 
                 saccadeDetectOut = saccadeDetection([calData.eyeXcal' calData.eyeYcal'],'sgolayall', 1,'sgolayfilter', 1,'sgolaywindow', 20,'detectionwindow', [120 180],'wellformed', false); %
               
                analysisResults(fileNum,1).saccadeDetection = saccadeDetectOut;
                            save(saccFileName,'saccadeDetectOut')
                            close all
                            % 'detectionwindow', [120 200],'sgolaywindow', 20,'fixedthreshold', 20
            else
                load(saccFileName)
                %saccadeDetectOut = saccadeDetectOut;
            end
            %TODO be careful with Marios position data it seems to
            %overshoot th real position at end of saccades, velocity is
            %much smoother though (see plot on line below)
           %  figure; plot(calData.eyeYcal,'b'); hold on; plot(saccadeDetectOut.traces.position.raw(:,2),'r')
            %now apped filtered traces to calData
            calData.eyeVelX = saccadeDetectOut.traces.velocity.xy(:,1)';
             calData.eyeVelY = saccadeDetectOut.traces.velocity.xy(:,2)';
              calData.eyeVelV = saccadeDetectOut.traces.velocity.v';
              calData.eyeAng = atan2(calData.eyeYcal,calData.eyeXcal);
            %allSaccStarts = (saccadeDetectOut.saccades.landmarks.v(:,2)/1000); %vector of saccade times for plotting
           allSaccStarts = (saccadeDetectOut.saccades.landmarks.v(:,1)/1000); %in old version I used to use the second column which is now apprentyl maximum velocity, so switching to 1 which is onset
            % allSaccPeaks = (saccadeDetectOut.saccades.landmarks.v(:,2)/1000);
            %plot the figure showing all the data from this file
             if mainPlotFlag
                
                hF = figure;
                hAe = axes('Parent',hF,'units','normalized','position',[0.05 0.55 0.9 0.4],...
                    'XLim',[alignedSpkTimeVector(1) alignedSpkTimeVector(end)]);
                hAs = axes('Parent',hF,'units','normalized','position',[0.05 0.05 0.9 0.4]);
                %plot spike data
                line(alignedSpkTimeVector,thisData.spkData,'parent',hAs,'color','k');
                linkaxes([hAe hAs],'x')
                %plot eye traces
                line(calData.eyeTcal,calData.eyeXcal,'parent',hAe,'color','r','tag','tEyeX')
                line(calData.eyeTcal,calData.eyeYcal,'parent',hAe,'color','b','tag','tEyeY')
                
                %plot unit markers
                for unitNum = thisFileInfo.unitNumbers
                    %spikes for this unit number aligned to start time -
                    %now done earlier
                    %unitSpikes = thisData.alignedSpikes{unitNum}-thisData.spkDataStartTime;
                    unitSpikes = spikeTimesCell{unitNum};
                    lineHeight = 1.3*max(thisData.spkData)+unitNum*20;
                    line(unitSpikes,lineHeight*ones(length(unitSpikes),1),'parent',hAs,'color',cc(unitNum,:),...
                        'linestyle','none','marker','s','MarkerFaceColor',cc(unitNum,:))
                    
                end
                
                %plotting
                lineHeight = 20;
                line(allSaccStarts,lineHeight*ones(length(allSaccStarts),1),'parent',hAe,'color','k',...
                    'linestyle','none','marker','s','MarkerFaceColor','k')
%                lineHeight = 25;
%                line(allSaccPeaks,lineHeight*ones(length(allSaccPeaks),1),'parent',hAe,'color','k',...
%                     'linestyle','none','marker','s','MarkerFaceColor','r')
                try
                    export_fig(figureFileName, '-pdf', '-append', hF);
                catch
                    export_fig(figureFileName, '-pdf', hF);
                end
                
                %save actual .fig to enable zooming
                hgsave(hF,[resultsDir filesep thisFileInfo.fileName(1:end-4) '-mainPlot.fig'])
                delete(hF)
               
             end
            
            %plot ISI and cross-triggerd psth
            if analyseISI

                %find which units actually have spikes on them
                validUnits = intersect(thisFileInfo.unitNumbers,find(~cellfun('isempty',thisData.alignedSpikes)));
                
                if length(validUnits)>1
                    %if more than 1 then launch the grid of isi's a triggered
                    %spiek triggered rasters
                    hFi = isiSubplotNd({thisData.alignedSpikes{validUnits}},'rastWinWidth',[0.1 0.1]);
                else
                    %if not then just produce an isi histogram
                    hFi = figure;
                    hAxi = axes('parent',hFi);
                    isiVec = diff(thisData.alignedSpikes{validUnits});
                    [isiHist bins] = hist(isiVec,500);
%                     isiHist = 100*isiHist/leng 
                    hFi = isiSubplotNd({thisData.alignedSpikes{validUnits}},'rastWinWidth',[0.1 0.1]);

                    xlabel('ISI (s)')
                    ylabel('Percent Spikes')
                    title(['ISI for unit ' num2str(validUnits)])
                end
                %try and export
                try
                    export_fig(figureFileName, '-pdf', '-append', hFi);
                catch
                    export_fig(figureFileName, '-pdf', hFi);
                end
                
                delete(hFi)
            end
            
            if checkCSpikePause
                if length(thisFileInfo.unitNumbers)>1
                %hF = figure;
                [hRast alignedRaster] = createRaster(spikeTimesCell{thisFileInfo.unitNumbers(1)},...
                    spikeTimesCell{thisFileInfo.unitNumbers(2)},[0.1 0.1],'spikeMarkSize',3); %TODO also make this plot other complex spikes
                 title(thisFileInfo.fileName(1:end-4))
                export_fig(figureFileName, '-pdf', '-append', gcf);
                delete(gcf)
             
                %look at spike shapes
               % for unitNum = thisFileInfo.unitNumbers
               tempSpikeTimesCell = cellfun(@(x) (x-thisData.spkDataStartTime),thisData.alignedSpikes,'uniformoutput',false);
            
                [relabelledSpikeShapes timeVector] = extractSpikeShapesFromTimes(thisData.spkData,tempSpikeTimesCell,[0.8 1.8],thisData.spkFs,...
                    'cc',lines(12),'doPlotFlag',1,'exSpikeShapeNum',50);
                title(thisFileInfo.fileName(1:end-4))
                export_fig(figureFileName, '-pdf', '-append', gcf);
                delete(gcf)
                end
            end
            
            
            %create new appended trial structure
             
            newTrialStructure = appendRastersToTrials(trialData,spikeTimesCell,saccadeDetectOut,calData,timeWindow,figureFileName);
             
            %also add lfp to trials, %not sure why sometimes is struct and
            %other cell
            if isstruct(neuronList(neuronNum).fileList)
                chNum = neuronList(neuronNum).fileList(fileNum).spkChannelNumber;
            else
           chNum = neuronList(neuronNum).fileList{fileNum,2};
            end
            %chNum = neuronList(neuronNum).fileList(fileNum).spkChannelNumber;
            newTrialStructure = appendLfpToTrials(newTrialStructure,thisFileInfo.fileName,...
                chNum,timeWindow);
            
            analysisResults(fileNum).newTrialStructure = newTrialStructure;
           
            analysisResults(fileNum).alignedSpikeTimesCell = spikeTimesCell;
            
         end
         
          save([resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd],'analysisResults')
         
     else
         
         disp('neuron already analysed')
     end
end
 
%% FILE By FILE BURST/PAUSE -UC
% allNeuronBurstFileName = [resultsDir filesep 'allNeuronBurstAnalysis2015.mat'];
% doBursts = 1
% if doBursts
% if exist(allNeuronBurstFileName,'file')~=2
%     structInitializer = cell(max(recordingsToAnalyse),1);
%     allNeuronBurstAnalysisResults = struct('neuronBurstResults',structInitializer);
%     for neuronNum = recordingsToAnalyse;
%         close all %TODO be careful!!
%         disp('Burst Pause Detection')
%         disp(['    processing neuron ' num2str(neuronNum) ' of ' num2str(max(recordingsToAnalyse))])
%         
%         %load matfile into memory
%         mLoadedData = matfile([processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd]) %create mat file
%         
%         
%         
%         numFiles = size(mLoadedData.editedFileList,1);
%         %create blank structure for holdind analysis results for this neuron
%         structInitializer = cell(numFiles,1); %has 2 comlums as second is for the sceond unit (complex spikes)
%         burstAnalysisResults = struct('burstResults',structInitializer);
%         %loop over all files and load the required data
%         
%         for fileNum = 1:numFiles
%             disp(['    processing file ' num2str(fileNum) ' of ' num2str(numFiles)])
%             
%             
%             %load data from mat file
%             thisData = mLoadedData.dataStructure(fileNum,1);
%             thisFileInfo = mLoadedData.editedFileList(fileNum,1);
%             thisFileContents = whos(mLoadedData);
%             
%             spikeTimesCell = cellfun(@(x) (x-thisData.spkDataStartTime),thisData.alignedSpikes,'uniformoutput',false);
%             
%             validUnits = intersect(thisFileInfo.unitNumbers,find(~cellfun('isempty',thisData.alignedSpikes)))
%             
%             
%             %loop over units and get the burst fro each
%             unitNum = validUnits(1); %TODO currently only for simple spikes
%             theseSpikeTimes = spikeTimesCell{unitNum};
%             
%             
%             [spikePatterns] = detectBurstPause(theseSpikeTimes);
%             %hgsave(gcf,[thisFileInfo.fileName '-bp.fig'])
%             delete(gcf)
%             
%             burstAnalysisResults(fileNum).burstResults = spikePatterns;
%         end
%         
%         allNeuronBurstAnalysisResults(neuronNum).neuronBurstResults = burstAnalysisResults;
%     end
%     
%     
%     save([resultsDir filesep 'allNeuronBurstAnalysis.mat'],'allNeuronBurstAnalysisResults')
%     
% else
%     
%     load('allNeuronBurstAnalysis.mat')
%     
% end 
% 
% %% FILE BY FILE BURST PAUSE RASTERS -UC
% 
% 
% 
% for neuronNum = recordingsToAnalyse;
%     mLoadedData = matfile([processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd]) %create mat file
%         
%         
%         
%         numFiles = size(mLoadedData.editedFileList,1);
%     for fileNum = 1:numFiles
%         theseBurstPauseResults = allNeuronBurstAnalysisResults(neuronNum).neuronBurstResults(fileNum);
%         load([resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd],'analysisResults')
% %             
%              thisTrialStruct = analysisResults(fileNum).newTrialStructure;
%              burstFigureFileName = [neuronList(neuronNum).neuronName '-file-' num2str(fileNum) '-burstTrialFigs20155.pdf'];
%         newTrialData = appendBurstPausesToTrials(thisTrialStruct,theseBurstPauseResults,[2 2],burstFigureFileName);
%         %TODO append burstPauseResults to trials
%         analysisResults(fileNum).newTrialStructure = newTrialData;
%         
%     end
%     save([resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd],'analysisResults')
% end
% 
% end
%% FILE BY FILE BURST DETECTION - OBS


%stInit = cell(1,numNeurons);
%allNeuronResults = struct('allStableBurstTrials',stInit);
analyseBursts = 0; %be careful, only turn off if burst already analysed for all files
if analyseBursts
for neuronNum = recordingsToAnalyse;
    close all %TODO be careful!!
    disp('Burst Detection')
    disp(['    processing neuron ' num2str(neuronNum) ' of ' num2str(max(recordingsToAnalyse))])
    
    %create file name for the summary figure pdf for this neuron
    burstFigureFileName = [resultsDir filesep neuronList(neuronNum).neuronName burstFigPdfFileNameEnd];
    thisNeuronSumFigName = [resultsDir filesep neuronList(neuronNum).neuronName indNeuronSumFigEnd];
    
    %load matfile into memory
    mLoadedData = matfile([processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd]) %create mat file
    load([resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd],'analysisResults') %also load the analysis results for this neuron
    
    %TODO add her if to check for saved burst results.
    %
    
    numFiles = size(mLoadedData.editedFileList,1);
    %create blank structure for holdind analysis results for this neuron
    structInitializer = cell(numFiles,1); %has 2 comlums as second is for the sceond unit (complex spikes)
    burstAnalysisResults = struct('burstResults',structInitializer);
    %loop over all files and load the required data
    
    for fileNum = 1:numFiles
        
        disp(['    processing file ' num2str(fileNum) ' of ' num2str(numFiles)])
        
        
        %load data from mat file
        thisData = mLoadedData.dataStructure(fileNum,1);
        thisFileInfo = mLoadedData.editedFileList(fileNum,1);
        thisFileContents = whos(mLoadedData);
        
        
        thisBurstFileName = [resultsDir filesep thisFileInfo.fileName(1:end-4) burstFileNameEnd];
        if exist(thisBurstFileName,'file')~=2
            
            
            spikeTimesCell = cellfun(@(x) (x-thisData.spkDataStartTime),thisData.alignedSpikes,'uniformoutput',false);
            
            validUnits = intersect(thisFileInfo.unitNumbers,find(~cellfun('isempty',thisData.alignedSpikes)))
            
            
            %loop over units and get the burst fro each
            unitNum = validUnits(1); %TODO currently only for simple spikes
            theseSpikeTimes = spikeTimesCell{unitNum};
            
            
            [archive_burst_RS,archive_burst_length,archive_burst_start]=rsBurstDetection(theseSpikeTimes);
            
            %                    hF = figure;
            %                    hAx = axes('parent',hF)
            %                    %now plot teh spike times and start marking bursts
            %                    line(theseSpikeTimes,ones(1,length(theseSpikeTimes)),'linestyle','none','marker','.')
            numBursts = length(archive_burst_start);
            stInit = cell(1,numBursts);
            burstStruct = struct('startSpike',stInit,'endSpike',stInit,'SI',stInit,'startTime',stInit,'endTime',stInit);
            
            for burstNum = 1:numBursts
                
                line([theseSpikeTimes(archive_burst_start(burstNum)) theseSpikeTimes(archive_burst_start(burstNum)+archive_burst_length(burstNum)-1)],[1.2 1.2],'color','r')
                %store burst in terms of seconds
                burstStruct(burstNum).startSpike = archive_burst_start(burstNum);
                burstStruct(burstNum).endSpike = archive_burst_start(burstNum)+archive_burst_length(burstNum);
                burstStruct(burstNum).SI = archive_burst_RS(burstNum);
                burstStruct(burstNum).startTime = theseSpikeTimes(archive_burst_start(burstNum));
                burstStruct(burstNum).endTime = theseSpikeTimes(archive_burst_start(burstNum)+archive_burst_length(burstNum)-1);
                
                
                
                
                
            end
            
            
            %alos plot inst fr to see what it looks like
            instFr = diff(theseSpikeTimes);
            instFr = 1./instFr;
            line(theseSpikeTimes,[0 ;instFr]./100,'color','k')
            
            
            
            
            
            %now attach burst times to trialStructure
            
            
            
            
            
            thisTrialStruct = analysisResults(fileNum).newTrialStructure;
            newTrialStruct = thisTrialStruct;
            %loop ove reach trial and if there is a burst within it add it
            %to the trialStruct.
            
            numTrials  = size(newTrialStruct,2);
            
            for trialNum =1:numTrials
                
                thisTrialStart = newTrialStruct(trialNum).trialStart;
                thisTrialEnd = newTrialStruct(trialNum).trialEnd;
                goCueTime = newTrialStruct(trialNum).bit3time;
                overlapFr = nan(1); %for peak of smoth firing rate during this time
                %check if it a valid trial
                if ~isempty(thisTrialStart) && ~isempty(thisTrialEnd) && ~isempty(goCueTime)
                    
                    %find detected bursts within the trial
                    theseBursts = find([burstStruct.startTime]>thisTrialStart & [burstStruct.endTime]<thisTrialEnd);
                    
                    if ~isempty(theseBursts)
                        
                        
                        
                        %get event times for this trial
                        relativeGoCueTime = newTrialStruct(trialNum).goCueTime;
                        saccadeTime = newTrialStruct(trialNum).saccadeTime;
                        
                        
                        binSize = 0.02; %20ms gaussian - %TODO check!!
                        testDt = 0.001; %1ms resolution
                        testTimeVector = -timeWindow(1):testDt:timeWindow(2);
                        %also use ksdensity to get a spike density function
                        if ~isempty(newTrialStruct(trialNum).alignedSpikes{1})
                        smoothFiringRate = ksdensity(newTrialStruct(trialNum).alignedSpikes{1},testTimeVector,'width',binSize);
                        
                        smoothFiringRate = smoothFiringRate*length(newTrialStruct(trialNum).alignedSpikes{1});
                        
                        else
                            smoothFiringRate = zeros(1,length(newTrialStruct(trialNum).alignedSpikes{1}));
                        end
                        if burstPlotFlag
                            %if not empty then plot as patches on top of the
                            %spikes
                            hF = figure;
                            hAx = axes('parent',hF);
                            %plot spike timing for both units %TODO include option for more units
                            line(newTrialStruct(trialNum).alignedSpikes{1},1*ones(1,length(newTrialStruct(trialNum).alignedSpikes{1})),'color','k','parent',hAx,...
                                'linestyle','none','marker','d','MarkerFaceColor','k')
                            line(newTrialStruct(trialNum).alignedSpikes{2},2*ones(1,length(newTrialStruct(trialNum).alignedSpikes{2})),'color','r','parent',hAx,...
                                'linestyle','none','marker','d','MarkerFaceColor','r')
                            legend('unit 1','unit 2')
                            ylim([0.5 2.5])
                            
                            
                            
                            
                            %now mark the various trial timing aspects
                            line([0 0],[0 1],'color','k','parent',hAx,...
                                'linestyle','-') %trial start
                            line([relativeGoCueTime relativeGoCueTime]-0.1,[0 1],'color','r','parent',hAx,...
                                'linestyle','-') %target apears
                            line([relativeGoCueTime relativeGoCueTime],[0 1],'color','r','parent',hAx,...
                                'linestyle','-') %go cue
                            line([saccadeTime saccadeTime],[0 1],'color','g','parent',hAx,...
                                'linestyle','-') %target hit according to stim computer
                            %ylim([0 1])
                            
                            %now for each burst in this trial draw a patch representing teh
                            %within burst time with height representing surprise index
                            
                            
                            
                            %switch to separatey axis
                            hAx2 = axes('parent',hF,'units','normalized','pos',get(hAx,'position'),'yaxislocation','right')
                            line(testTimeVector,smoothFiringRate,'color','k','parent',hAx2)
                            
                            set(hAx2,'color','none')
                            %now find only burst that are within the saccade
                            %related period. currenty 100ms either side of
                            %saccade onset but if chanignt hte also cahge in
                            %kunimatsuAnalysis hich produces the ku iamtsu
                            %significance results
                            
                        end
                        
                        saccadePeriodTv = saccadeTime+thisTrialStart+[-0.1:0.001:0.1]; %for testing overlap of bursts
                        numTrialBursts = length(theseBursts);
                        %overlapLog = nan(1,numTrialBursts); %for recording if overalpping with saccade perid
                        
                        %overlapFrPeak = nan(1,numTrialBursts); %TODO doesn't need to be 1,x.
                        
                        for trialBurst = 1:numTrialBursts
                            
                            burstStartAligned = burstStruct(theseBursts(trialBurst)).startTime-thisTrialStart
                            burstEndAligned = burstStruct(theseBursts(trialBurst)).endTime-thisTrialStart
                            
                            %cehck if burst overlaps with saccade period
                            burstTv = [burstStruct(theseBursts(trialBurst)).startTime:0.001:burstStruct(theseBursts(trialBurst)).endTime];
                            
                            
                            
                            if ~isempty(intersect(roundTo(saccadePeriodTv,0.001),roundTo(burstTv,0.001)))
                                
                                if burstPlotFlag
                                    patch([burstEndAligned burstEndAligned burstStartAligned burstStartAligned],...
                                        [1.5 2.5 2.5 1.5],'m','parent',hAx,'facealpha',0.5)
                                end
                                %overlapLog(trialBurst) = 1;
                                
                                saccadeInds = [saccadePeriodTv(1)-thisTrialStart saccadePeriodTv(end)-thisTrialStart+0.5]*1000;
                                %TODO for now take maximum in the burst and
                                %tehn if more tahn one burst then take max
                                %of all
                                %below is to deal with strange trial that
                                %has saccde ate very end. TODO use a debug point
                                %to understand
                                if saccadeInds(2)>length(smoothFiringRate)
                                    saccadeInds(2) = length(smoothFiringRate);
                                end
                                
                                thisBurstFr = max(smoothFiringRate(saccadeInds(1):saccadeInds(2)));
                                
                                if  ~isnan(overlapFr)
                                    %TODO need to get timing as well
                                    overlapFr = max(thisBurstFr,overlapFr);
                                else
                                    overlapFr = thisBurstFr;
                                    %overlapBurst = trialBurst; %the number of burst used
                                end
                                
                            else
                                if burstPlotFlag
                                    patch([burstEndAligned burstEndAligned burstStartAligned burstStartAligned],...
                                        [1.5 2.5 2.5 1.5],'r','parent',hAx,'facealpha',0.5)
                                end
                                %overlapLog(trialBurst) = 0;
                            end
                        end
                        
                        if burstPlotFlag
                            text(0.5,1.5,num2str(overlapFr),'parent',hAx)
                            
                            
                            title(['Trial : ' num2str(trialNum)])
                            
                            export_fig(burstFigureFileName, '-pdf', '-append', hF);
                            delete(hF)
                        end
                        
                        
                        newTrialStruct(trialNum).incBurstStruct = burstStruct(theseBursts);
                        newTrialStruct(trialNum).overlapFr = overlapFr;
                        
                    end
                    
                    
                    
                end
                
                
            end
            
            burstAnalysisResults = newTrialStruct; %TODO fix name and join to other analysis
            save(thisBurstFileName,'burstAnalysisResults')
            
            
        else
            load(thisBurstFileName)
            
        end
        analysisResults(fileNum).newTrialStructure = burstAnalysisResults;
    end
    
    save([resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd],'analysisResults')
         
end 
end
%%  NEURON BY NEURON ANALYSIS
%loop over each neuron and stick all trials together, then remove those
%that are ouside the stable period marked in editedFileList


smallTimeWindow = [0.3 0.5]; %time window for actual rasters


 if exist([resultsDir filesep wholeNeuronNameEnd],'file')~=2 %delete the file ans start again every time
     delete([resultsDir filesep wholeNeuronNameEnd])
 end
    %TODO initialise to correct size
    wholeNeuronResults = struct('allStableTrials',[],'selectedTrials',[],'monkey',[],'area',[],'depth',[],'gridLoc',[]); %allStableTrails contaiins the trial structure for all trials in a neuron within the periods marked as stabel by Eric.
    
    
    
     % flag 666 
    %selectedTrials is a struct containing files with different groups of trial
    %numbers (eg correct trials)
    
    for neuronNum = recordingsToAnalyse;
    
        
        wholeNeuronResults(neuronNum).monkey = neuronList(neuronNum).monkey;
        wholeNeuronResults(neuronNum).area = neuronList(neuronNum).area;
        wholeNeuronResults(neuronNum).depth = neuronList(neuronNum).depth;
        wholeNeuronResults(neuronNum).gridLoc = neuronList(neuronNum).gridLocation;

        
        
               %load matfile into memory
     thisMatFileName = [processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd];
     mLoadedData = matfile(thisMatFileName); %create mat file
    
        %build filename for analysis
        thisAnalysisFileName = [resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd];
        
        %build figure filename
        thisNeuronSumFigName = [resultsDir filesep neuronList(neuronNum).neuronName indNeuronSumFigEnd];
        
        %     %build filename for burst results
        %     thisNeuronBurstFileName =
        
        numFiles = size(mLoadedData.editedFileList,1);
        load(thisAnalysisFileName);
        stablePeriodTrials = cell(1,numFiles);
        
        
        %loop over files and delete trials falling otuside marked times
        for fileNum = 1:numFiles
            
            %get file info
            thisFileInfo = mLoadedData.editedFileList(fileNum,1);
            
            %get start and end of stable section in this file
            if isinf(thisFileInfo.startTime)
                thisStartTime = -inf;
            else
                thisStartTime = thisFileInfo.startTime;
            end
            thisStopTime = thisFileInfo.stopTime;
            
            
            %list all trial starts falling within these
            trialsToKeep = cellfun(@(x) x>thisStartTime & x<thisStopTime,{analysisResults(fileNum).newTrialStructure.trialStart},'uniformoutput',false);
            stablePeriodTrials{fileNum} = analysisResults(fileNum).newTrialStructure([trialsToKeep{:}]);
            
            %stablePeriodTrials{fileNum} =  burstAnalysisResults(fileNum).burstResults;
        end
        
        %concatenate all together
        allStableTrials = [stablePeriodTrials{:}];
        wholeNeuronResults(neuronNum).allStableTrials  = allStableTrials;
        
        
        %find all trails with a valid saccade and a correct response
        
        tempLogCell = cellfun(@(x) x==2,{allStableTrials.correctResponse},'uniformoutput',false);
        %tempLogCell2 = cellfun(@(x) any(x),{allStableTrials.bit8time},'uniformoutput',false);
        
        correctTrialNums = find([tempLogCell{:}]); %TOOD this is an error check what is supposed to happen
        %correctTrialNums = find([tempLogCell2{:}]);
         tempLogCell = cellfun(@(x) x==1,{allStableTrials.correctResponse},'uniformoutput',false);
        %tempLogCell2 = cellfun(@(x) any(x),{allStableTrials.bit8time},'uniformoutput',false);
        
        wrongTrialNums = find([tempLogCell{:}]);
        
        
          %now we need to remove any trials which are labelled as corSac but
        %don't contain a bit8
       theseBit8Times = {allStableTrials(correctTrialNums).bit8time};
      badBit8trials =  find(cellfun('isempty',theseBit8Times));
      if ~isempty(badBit8trials)
          correctTrialNums = setdiff(correctTrialNums,correctTrialNums(badBit8trials));
         
      end
       
        find(cellfun('isempty',{allStableTrials(correctTrialNums).bit8time}))
        
        a=1;
        
        tempLogCell = cellfun(@(x) any(x),{allStableTrials.saccadeTime},'uniformoutput',false);
        
        validSaccadeTrialNums = find([tempLogCell{:}]); 
        corSacTrialNums = intersect(correctTrialNums,validSaccadeTrialNums);
        wroSacTrialNums = intersect(wrongTrialNums,validSaccadeTrialNums);
        
      
        %now split by condition (pro v anti to start)
        tempLogCell = cellfun(@(x) ismember(x,trialConditionCodes(:,1)),{allStableTrials.conditionCode},'uniformoutput',false);
        %all pro-saccades
        proTrialNums = find([tempLogCell{:}]);
        %correct pro sacccades
        corProTrialNums = intersect(corSacTrialNums,proTrialNums);
         wroProTrialNums = intersect(wroSacTrialNums,proTrialNums);
        %all anti-saccades
        antiTrialNums = find([tempLogCell{:}]~=1);
        %correct anti saccades
        corAntiTrialNums = intersect(corSacTrialNums,antiTrialNums);
       wroAntiTrialNums = intersect(wroSacTrialNums,antiTrialNums);
        
        wholeNeuronResults(neuronNum).selectedTrials.proTrials = proTrialNums;
        wholeNeuronResults(neuronNum).selectedTrials.corProTrials = corProTrialNums;
        wholeNeuronResults(neuronNum).selectedTrials.antiTrials = antiTrialNums;
        wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials = corAntiTrialNums;
        wholeNeuronResults(neuronNum).selectedTrials.correctTrials = correctTrialNums;
        wholeNeuronResults(neuronNum).selectedTrials.validSaccadeTrials = validSaccadeTrialNums;
        wholeNeuronResults(neuronNum).selectedTrials.corSacTrials = corSacTrialNums; %correct trials with valid saccade
        wholeNeuronResults(neuronNum).selectedTrials.wroProTrials = wroProTrialNums;
        wholeNeuronResults(neuronNum).selectedTrials.wroAntiTrials = wroAntiTrialNums;
        wholeNeuronResults(neuronNum).selectedTrials.wroSacTrials = wroSacTrialNums;
        clear analysisResults
        
        
    end
    save([resultsDir filesep wholeNeuronNameEnd],'wholeNeuronResults','-v7.3')
    
%else
    
 %   load([resultsDir filesep wholeNeuronNameEnd])
%end
%% CHECK BIT 8 PROBLEM
behaviourPdfName = [resultsDir filesep 'bitproblemSummaryTestnewJan-removed.pdf'];
 numDirections = 8;
 validDirectionRows = find(~isnan(trialConditionCodes(:,1)))';
allWrongBitTrials = nan(length(recordingsToAnalyse),2);
recCounter = 1;
for neuronNum = recordingsToAnalyse
   
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
   corProTrials =  wholeNeuronResults(neuronNum).selectedTrials.corProTrials
    corAntiTrials =  wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials
   theseBit8timesPro = {theseStableTrials(corProTrials).bit8time};
   
  allWrongBitTrials(recCounter,1) = numel(find(cellfun('isempty',theseBit8timesPro)))./numel(corProTrials)
    theseBit8timesAnti = {theseStableTrials(corAntiTrials).bit8time};
    allWrongBitTrials(recCounter,2) = numel(find(cellfun('isempty',theseBit8timesAnti)))./numel(corAntiTrials)
    
    
    
    
     hF = figure;
    for directionNum = 1:numDirections
        
         plNum = validDirectionRows(directionNum);
         hA(directionNum) = subaxis(3,3,plNum,'SpacingVert',0.2,'MR',0.05);
%           theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
%         theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        thisProCode = trialConditionCodes(validDirectionRows(directionNum),1);
            thisAntiCode = trialConditionCodes(validDirectionRows(directionNum),2);
            
            
        tempLogCell = cellfun(@(x) x==thisProCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirProTrialNums = find([tempLogCell{:}]);
          %  thisDirCorProTrialNums = intersect(allProTrialNums,thisDirProTrialNums);
            
          thisDirProCorVec =  [theseStableTrials(thisDirProTrialNums).correctResponse];
            
            tempLogCell = cellfun(@(x) x==thisAntiCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirAntiTrialNums = find([tempLogCell{:}]);
        %    thisDirCorAntiTrialNums = intersect(allAntiTrialNums,thisDirAntiTrialNums);
            
         thisDirAntiCorVec =    [theseStableTrials(thisDirAntiTrialNums).correctResponse];
         
          thisDirProCorTrialNums = thisDirProTrialNums(thisDirProCorVec==2);
         thisDirProCorTrialNums =  intersect(thisDirProCorTrialNums,corProTrials)
         theseBit8timesPro = {theseStableTrials(thisDirProCorTrialNums).bit8time};
         wrongBitproTrials = numel(find(cellfun('isempty',theseBit8timesPro)))/numel(thisDirProCorTrialNums)
            
         thisDirAntiCorTrialNums = thisDirAntiTrialNums(thisDirAntiCorVec==2)
          thisDirAntiCorTrialNums =  intersect(thisDirAntiCorTrialNums,corAntiTrials)
         theseBit8timesAnti = {theseStableTrials(thisDirAntiCorTrialNums).bit8time};
         wrongBitantiTrials = numel(find(cellfun('isempty',theseBit8timesAnti)))/numel(thisDirAntiCorTrialNums)
         
         
         set(hA(directionNum),'xlim',[0 4],'ylim',[0 4])
        % for corCode = 0:2 
%              text(1,corCode+1,num2str(corCode),'color','k')
%              text(2,corCode+1,num2str(sum(thisDirProCorVec==corCode)),'color','g')
             text(1,1,num2str(wrongBitproTrials),'color','g')
             text(3,1,num2str(numel(thisDirProCorTrialNums)),'color','g')
             text(1,2,num2str(wrongBitantiTrials),'color','r') 
             text(3,2,num2str(numel(thisDirAntiCorTrialNums)),'color','r')
         %end
         
         if isnan(wrongBitantiTrials)
             
             a = 1;
         end
        
    end
    
    
     suptitle(strrep(neuronList(neuronNum).neuronName,'_','-'))
    
    
    export_fig(hF,'-pdf','-append',behaviourPdfName)
        delete(hF)
    
     recCounter =  recCounter+1;
end
%% BEHAVIOURAL PERFORMANCE
 numDirections = 8;
 validDirectionRows = find(~isnan(trialConditionCodes(:,1)))';
 behaviourPdfName = [resultsDir filesep 'behaviourSummaryTestnewJan.pdf'];
for neuronNum = recordingsToAnalyse;
    
    
    %get all correct trial numbers
%         theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
%         theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        
        theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
        theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
    
    %get all pro trials
    allProTrialNums  = wholeNeuronResults(neuronNum).selectedTrials.proTrials;
    allAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.antiTrials;
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
    
    hF = figure;
    for directionNum = 1:numDirections
        
         plNum = validDirectionRows(directionNum);
         hA(directionNum) = subaxis(3,3,plNum,'SpacingVert',0.2,'MR',0.05);
%           theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
%         theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        thisProCode = trialConditionCodes(validDirectionRows(directionNum),1);
            thisAntiCode = trialConditionCodes(validDirectionRows(directionNum),2);
            
            
        tempLogCell = cellfun(@(x) x==thisProCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirProTrialNums = find([tempLogCell{:}]);
          %  thisDirCorProTrialNums = intersect(allProTrialNums,thisDirProTrialNums);
            
          thisDirProCorVec =  [theseStableTrials(thisDirProTrialNums).correctResponse];
            
            tempLogCell = cellfun(@(x) x==thisAntiCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirAntiTrialNums = find([tempLogCell{:}]);
        %    thisDirCorAntiTrialNums = intersect(allAntiTrialNums,thisDirAntiTrialNums);
            
         thisDirAntiCorVec =    [theseStableTrials(thisDirAntiTrialNums).correctResponse];
         
         set(hA(directionNum),'xlim',[0 4],'ylim',[0 4])
         for corCode = 0:2 
             text(1,corCode+1,num2str(corCode),'color','k')
             text(2,corCode+1,num2str(sum(thisDirProCorVec==corCode)),'color','g')
             text(3,corCode+1,num2str(sum(thisDirAntiCorVec==corCode)),'color','r')
             
         end
        
    end
    suptitle(strrep(neuronList(neuronNum).neuronName,'_','-'))
    
    export_fig(hF,'-pdf','-append',behaviourPdfName)
        delete(hF)
end
%% SACCADE 
doSaccPaths= 1;
if doSaccPaths
    
 conditionLegend = load('conditionLegend.mat');
for neuronNum = recordingsToAnalyse
    
    %this looks very strange so let's do it by target and check what's
    %going on.
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
    theseCorTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
    for stimCode = 1:16;
        hF = figure;
        hA = axes('parent',hF);
        tempLogCell2 = cellfun(@(x) (x==stimCode),{theseStableTrials.conditionCode},'uniformoutput',false);
        thisCodeTrialNums = find([tempLogCell2{:}]);
        corSacTrialNums = intersect(thisCodeTrialNums,theseCorTrialNums);
        theseTrialNums = corSacTrialNums;
        numTrials = length(theseTrialNums);
        for trialNum = 1:numTrials
            %whats the +4000 for????? to get back to middle of trial which
            %is -4 4 s around saccade initiation

            
            thisTrial = theseStableTrials(theseTrialNums(trialNum));
            thisTrialSaccadeInds = round([thisTrial.saccadeTime*1000:thisTrial.saccadeTime*1000+thisTrial.saccadeDuration]+4000);
            
            %also extract extended version for plotting with dotted line
            extensionSize = [500 500]; %ms either side to extend by
            thisTrialPreExtendedSaccadeInds = round([thisTrial.saccadeTime*1000-extensionSize(1):thisTrial.saccadeTime*1000+thisTrial.saccadeDuration]+4000);
            thisTrialPostExtendedSaccadeInds = round([thisTrial.saccadeTime*1000:thisTrial.saccadeTime*1000+thisTrial.saccadeDuration+extensionSize(2)]+4000);
           
            %trim to 8000 to stop going over the number of samples
            thisTrialPostExtendedSaccadeInds = thisTrialPostExtendedSaccadeInds(thisTrialPostExtendedSaccadeInds<8000);
            thisTrialPreExtendedSaccadeInds = thisTrialPreExtendedSaccadeInds(thisTrialPreExtendedSaccadeInds>0);
            
            thisTrialBit8relativeTime = thisTrial.bit8time-thisTrial.trialStart;
            thisTrialBit8timeind = round(thisTrialBit8relativeTime*1000)+4000;
            
            
            
                % bit8 should be between begin and end of saccade
            if thisTrialBit8timeind>=thisTrialSaccadeInds(1) && thisTrialBit8timeind<=thisTrialSaccadeInds(end)
                lCol = 'b';
                line([thisTrial.eyePositionX(thisTrialSaccadeInds(thisTrialBit8timeind-thisTrialSaccadeInds(1)+1)) NaN],...
                [thisTrial.eyePositionY(thisTrialSaccadeInds(thisTrialBit8timeind-thisTrialSaccadeInds(1)+1)) NaN],...
                'parent',hA,'marker','o','markerfacecolor',lCol)
           
            else
                 lCol = 'r';
                
            end
            %plot main sequence
            line(thisTrial.eyePositionX(thisTrialSaccadeInds),thisTrial.eyePositionY(thisTrialSaccadeInds),'parent',hA,'color',lCol);
            %plot extended version in dotted line, black for pre
            %saccade,line colour specified by correct trial for post
           
            line(thisTrial.eyePositionX(thisTrialPostExtendedSaccadeInds),thisTrial.eyePositionY(thisTrialPostExtendedSaccadeInds),'parent',hA,'color',lCol,'linestyle',':');
           line(thisTrial.eyePositionX(thisTrialPreExtendedSaccadeInds),thisTrial.eyePositionY(thisTrialPreExtendedSaccadeInds),'parent',hA,'color','k','linestyle',':');
           
            
            %now addd a square for start and triangle for end;
            line([thisTrial.eyePositionX(thisTrialSaccadeInds(1)) NaN],[thisTrial.eyePositionY(thisTrialSaccadeInds(1)) NaN],'parent',hA,'marker','s','markerfacecolor',lCol)
            line([thisTrial.eyePositionX(thisTrialSaccadeInds(end)) NaN],[thisTrial.eyePositionY(thisTrialSaccadeInds(end)) NaN],'parent',hA,'marker','^','markerfacecolor',lCol)
            thisTrialTarget = thisTrial.conditionCode;
            thisTargetPosition =  [conditionLegend.test{thisTrialTarget,2} conditionLegend.test{thisTrialTarget,3}];
            drawCircle(thisTargetPosition(1),thisTargetPosition(2),1,hA)
            
            
            
            
        end
        title(['Target' num2str(stimCode)])
        export_fig([resultsDir filesep neuronList(neuronNum).neuronName 'saccadeTestingbit8-newCalib-extendedBothv2-newcalib-jan.pdf'], '-pdf','-append', hF);
        hgsave(hF,[resultsDir filesep neuronList(neuronNum).neuronName '-targ' num2str(stimCode) 'saccadeTestingbit8-newCalib-extendedBothv2-newcalib.fig'])
%         delete(hF)
        
    end
end
end

%% NEURON BY NEURON RECALIBRATION
%as the calibration from earlier doesn't seem perfect I will now use the
%saccade end points from the marios saccade detection and then retune it
%now.
%use only pro correct trials 
if doSaccPaths
newCalFileNameEnd = '-postSacDetectRecal.mat';
newCalPdfNameStart = [resultsDir filesep 'calibrationTestingWholeNeuron-'];
proCodes = trialConditionCodes(~isnan(trialConditionCodes(:,1)),1);
proCodes = proCodes(:)';
targCc = lines(length(proCodes)); 
for neuronNum = recordingsToAnalyse
    
    
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
    theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
    hF =figure;
    hA = axes('parent',hF);
    stimCodeCount = 1;
    endPoints = nan(100,2,length(proCodes)); %for storing end points:(trails,[x y],)
    for stimCode = proCodes
        
        
        tempLogCell2 = cellfun(@(x) (x==stimCode),{theseStableTrials.conditionCode},'uniformoutput',false);
        thisCodeTrialNums = find([tempLogCell2{:}]);
        theseTrialNums = intersect(theseCorProTrialNums,thisCodeTrialNums);
        %now let's loop over the pro codes in turn and make a plot for each
        
        numTrials = length(theseTrialNums);
        for trialNum = 1:numTrials
            thisTrial = theseStableTrials(theseTrialNums(trialNum));
            %now we want to find the start of the saccade and check it is
            %at 0,0 and find the end of the saccade and mark it in relation
            %to the target
            saccadeEnd = [thisTrial.saccadeOriginX+thisTrial.saccadeAmplitudeX thisTrial.saccadeOriginY+thisTrial.saccadeAmplitudeY];
            line([thisTrial.saccadeOriginX saccadeEnd(1)],...
                [thisTrial.saccadeOriginY saccadeEnd(2)],...
                'color',targCc(stimCodeCount,:),'marker','s','markerfacecolor',targCc(stimCodeCount,:),...
                'parent',hA)
            
            endPoints(trialNum,:,stimCodeCount) = saccadeEnd;
        end
        stimCodeCount = stimCodeCount+1;
    end
    stimCodeCount = 1;
    %now let's draw circles representing the targets
    muMagErs = nan(length(proCodes),2); %to store the mean magntude error for x and y for each target
    for stimCode = proCodes
        
        thisTargetPosition =  [conditionLegend.test{stimCode,2} conditionLegend.test{stimCode,3}];
        hC =  drawCircle(thisTargetPosition(1),thisTargetPosition(2),1,hA)
        set(hC,'color',targCc(stimCodeCount,:))
        
        
        
         %now get the x and y absolute errror for each saccade
         theseSacEnds = endPoints(:,:,stimCodeCount)
         theseXErs = thisTargetPosition(1)./theseSacEnds(:,1);
          theseYErs = thisTargetPosition(2)./theseSacEnds(:,2);
           %if the target has 0 position in this dimension then we shou;d
         %ignore it from the end calculation
        noPos =  find(thisTargetPosition==0);
          
         muMagErs(stimCodeCount,1) = nanmean(theseXErs);
         muMagErs(stimCodeCount,2) = nanmean(theseYErs);
        muMagErs(stimCodeCount,noPos) = nan;
         text(thisTargetPosition(1),thisTargetPosition(2),num2str(muMagErs(stimCodeCount,:)))
        stimCodeCount = stimCodeCount+1;
    end
    title([neuronList(neuronNum).neuronName '-Pre Cal'])
     export_fig([newCalPdfNameStart '.pdf'], '-pdf','-append', hF);
       delete(hF)
       
       
       %now lets get the mean in each dimension and stretch the saccades by
       %that amount to see alignment then, this time include incorrect
       %saccades to see if they miss
  
   magFactors = nanmean(muMagErs);
       
   save([processDir filesep neuronList(neuronNum).neuronName newCalFileNameEnd],'magFactors')
   
   
       theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
     theseProTrialNums =  intersect(wholeNeuronResults(neuronNum).selectedTrials.proTrials,wholeNeuronResults(neuronNum).selectedTrials.validSaccadeTrials)
    %above includes incorrect trial but all have valid saccades
    hF =figure;
    hA = axes('parent',hF);
    stimCodeCount = 1;
     recalEndPoints = nan(100,2,length(proCodes)); %for storing end points:(trails,[x y],)
    for stimCode = proCodes
        
        
        tempLogCell2 = cellfun(@(x) (x==stimCode),{theseStableTrials.conditionCode},'uniformoutput',false);
        thisCodeTrialNums = find([tempLogCell2{:}]);
        theseTrialNums = intersect(theseProTrialNums,thisCodeTrialNums);
        %now let's loop over the pro codes in turn and make a plot for each
        
        numTrials = length(theseTrialNums);
        for trialNum = 1:numTrials
            thisTrial = theseStableTrials(theseTrialNums(trialNum));
            %now we want to find the start of the saccade and check it is
            %at 0,0 and find the end of the saccade and mark it in relation
            %to the target
            saccadeEnd = [thisTrial.saccadeOriginX+thisTrial.saccadeAmplitudeX thisTrial.saccadeOriginY+thisTrial.saccadeAmplitudeY];
            
             saccadeEnd(1) =  saccadeEnd(1)*magFactors(1);
             saccadeEnd(2) =  saccadeEnd(2)*magFactors(2);
            if thisTrial.correctResponse==2
                markType = 's';
            else
                 markType = 'o';
            end
            
            line([thisTrial.saccadeOriginX saccadeEnd(1)],...
                [thisTrial.saccadeOriginY saccadeEnd(2)],...
                'color',targCc(stimCodeCount,:),'marker',markType,'markerfacecolor',targCc(stimCodeCount,:),...
                'parent',hA)
            
            recalEndPoints(trialNum,:,stimCodeCount) = saccadeEnd;
        end
        stimCodeCount = stimCodeCount+1;
    end
    
     stimCodeCount = 1;
    %now let's draw circles representing the targets
    muMagErsRecal = nan(length(proCodes),2); %to store the mean magntude error for x and y for each target
    for stimCode = proCodes
        
        thisTargetPosition =  [conditionLegend.test{stimCode,2} conditionLegend.test{stimCode,3}];
        hC =  drawCircle(thisTargetPosition(1),thisTargetPosition(2),1,hA)
        set(hC,'color',targCc(stimCodeCount,:))
        
        
        
         %now get the x and y absolute errror for each saccade
         theseSacEnds =  recalEndPoints(:,:,stimCodeCount)
         theseXErs = thisTargetPosition(1)./theseSacEnds(:,1);
          theseYErs = thisTargetPosition(2)./theseSacEnds(:,2);
           %if the target has 0 position in this dimension then we shou;d
         %ignore it from the end calculation
        noPos =  find(thisTargetPosition==0);
          
         muMagErsRecal(stimCodeCount,1) = nanmean(theseXErs);
         muMagErsRecal(stimCodeCount,2) = nanmean(theseYErs);
        muMagErsRecal(stimCodeCount,noPos) = nan;
         text(thisTargetPosition(1),thisTargetPosition(2),num2str(muMagErsRecal(stimCodeCount,:)))
        stimCodeCount = stimCodeCount+1;
    end
    title([neuronList(neuronNum).neuronName '-Post Cal'])
     export_fig([newCalPdfNameStart '.pdf'], '-pdf','-append', hF);
       delete(hF)
end
end
%% RASTERS FOR CORRECT/WRONG trials
doCorWroRast = 0;
if doCorWroRast
for neuronNum = recordingsToAnalyse;
    for trialType = {'Sac','Pro','Anti'};
        
        thisSumFigName = [resultsDir filesep 'summaryRastersFor_' char(trialType) 'Resorttrialorder2016.pdf'];
        %              if exist(thisSumFigName,'file')==2
        %                delete(thisSumFigName)
        %            end
        frAxisLim = nan(2,2);
        hF = figure; % plot correct trials on left and worng on right, need subplot for eye and spikes so 2x2 (1 is eye for cor, 3 is spikes)
        hSa = nan(3,2);
        hSa(1,1) = subplot(2,2,1);
        hSa(1,2) = subplot(2,2,2);
        hSa(2,1) = subplot(2,2,3);
        hSa(2,2) = subplot(2,2,4);
        %two hidden axes to enable smoothed rate plotting over the raster
        hSa(3,1) = axes('XAxisLocation','bottom','YAxisLocation','right','color','none',...
            'Xcolor',[1 1 1],'Ycolor',[1 0 0],'XTick',[],'parent',hF,'position',get(hSa(2,1),'position'));
        
        hSa(3,2) = axes('XAxisLocation','bottom','YAxisLocation','right','color','none',...
            'Xcolor',[1 1 1],'Ycolor',[1 0 0],'XTick',[],'parent',hF,'position',get(hSa(2,2),'position'));
        
        corTrialFieldString = ['cor' char(trialType) 'Trials'];
        wroTrialFieldString = ['wro' char(trialType) 'Trials'];
        if numel(wholeNeuronResults(neuronNum).selectedTrials.(corTrialFieldString))>=1
            [hAx1 hAx2 hA2] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
                wholeNeuronResults(neuronNum).selectedTrials.(corTrialFieldString),smallTimeWindow,timeWindow,'alignTo','saccadeTime','markEvent','goCueTime','sortBy','');%,...
            
            frAxisLim(1,:) =  get(hA2,'ylim')
            newHand = copyobj(allchild(hAx2),hSa(2,1));
            newHand = copyobj(allchild(hA2),hSa(3,1));
            newHand = copyobj(allchild(hAx1),hSa(1,1));
            %TODO give a 3,1 list of axis handels to rasterFromTrials
            set(hSa(2,1),'ylim',[0 numel(wholeNeuronResults(neuronNum).selectedTrials.(corTrialFieldString))])
            delete(get(hAx2,'parent'))
        end
        if numel(wholeNeuronResults(neuronNum).selectedTrials.(wroTrialFieldString))>=1
            [hAx1 hAx2 hA2] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
                wholeNeuronResults(neuronNum).selectedTrials.(wroTrialFieldString),smallTimeWindow,timeWindow,'alignTo','saccadeTime','markEvent','goCueTime','sortBy','');%,...
            newHand = copyobj(allchild(hAx2),hSa(2,2));
            newHand = copyobj(allchild(hA2),hSa(3,2));
            newHand = copyobj(allchild(hAx1),hSa(1,2));
            
            frAxisLim(2,:) =  get(hA2,'ylim')
            set(hSa(2,2),'ylim',[0 numel(wholeNeuronResults(neuronNum).selectedTrials.(wroTrialFieldString))])
            delete(get(hAx2,'parent'))
        end
        
        if sum(isnan(frAxisLim(:)))~=4
            set(get(hSa(1,1),'title'),'string','Correct')
            set(get(hSa(1,2),'title'),'string','Wrong')
            set(hSa(:),'xlim',[-smallTimeWindow(1) smallTimeWindow(2)])
            set(hSa(1,:),'ylim',[-15 15])
            set(hSa(3,:),'ylim',[min(frAxisLim(:,1)) max(frAxisLim(:,2))])
            
            suptitle([neuronList(neuronNum).neuronName char(trialType) ' trials, align to saccades'])
            try
                export_fig(hF,'-pdf','-append',thisSumFigName)
            catch
                export_fig(hF,'-pdf',thisSumFigName)
            end
            
        end
        delete(gcf)
        
    end
    
end
end
%% HERZFELD BURST/PAUSE classification

%for the lateral neurons especially the baseline used by herzfeld captures
%the build up of activity, so will switch to be 400ms:100ms before saccade
bpClassFigureFileName = [resultsDir filesep 'bpClasstestNewbaseallsaccsBestEarly baseRealigned.pdf'];
bpClassFileName = [resultsDir filesep 'bpClassificationNico.mat'];
if exist(bpClassFileName,'file')~=2

    allNeuronClassification = cell(max(recordingsToAnalyse),1); % 'burst', 'pause', 'not sig'
%     allNeuronClassification = cell(max(vermNeurons),1); 
    
for neuronNum = recordingsToAnalyse;
    %for each trial take SS firing rate in period 200:-100ms before saccade
    %onset, compare to 150ms window centered on saccade midpoint (they use
    %time of peak velocity %TODO get from saccade detection).
    
   thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
     theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
     theseCorSacTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
     numTrials = length(theseCorSacTrialNums);
     baseAndSaccRates = nan(numTrials,2); %first baseline spike rate, second perisacccade spike rate.
     
     for trialNum = 1:numTrials
        thisTrialNum = theseCorSacTrialNums(trialNum); 
         thisTrial = theseStableTrials(thisTrialNum);
         
          thisTrialSs = thisTrial.alignedSpikes{thisNeuronUnits(1)};
         baselineSpikeCount = sum(thisTrialSs>thisTrial.saccadeTime-0.4 & thisTrialSs<thisTrial.saccadeTime-0.2);
         baseAndSaccRates(trialNum,1) = baselineSpikeCount/0.2;
         saccMidpoint =thisTrial.saccadeLandmarks(3)/1000-thisTrial.trialStart; %third saccade landmark is moment of deceleration
         %saccMidpoint =
         %thisTrial.saccadeTime+(thisTrial.saccadeDuration/2)/1000; %htis
         %was th old halfway point of saccade method
       %  periSacSpikesCount = sum(thisTrialSs>saccMidpoint-0.075 & thisTrialSs<saccMidpoint+0.075);%herzfeld uses 150ms window we will try just the first 75ms to get more sig
         periSacSpikesCount = sum(thisTrialSs>saccMidpoint-0.075 & thisTrialSs<saccMidpoint+0.075);
        baseAndSaccRates(trialNum,2) = periSacSpikesCount/0.15;
     end
     
     %now do a paired ttest to see if different and also if higher or
     %lower.
     [h p] = ttest(baseAndSaccRates(:,1),baseAndSaccRates(:,2));
     if p<0.05
         if mean(baseAndSaccRates(:,1)-baseAndSaccRates(:,2))<0
              allNeuronClassification{neuronNum} = 'burst';
         else
         allNeuronClassification{neuronNum} = 'pause';
         end
     else
         allNeuronClassification{neuronNum} = 'not sig';
     end
     
     %we are not finding many signifcant cells so let's try checking why
     hF = figure;
     hA = axes('parent',hF);
     scatter(baseAndSaccRates(:,1),baseAndSaccRates(:,2))
     line([0 150],[0 150],'color','k')
     xlabel('base')
    ylabel('sacc')

     title([neuronList(neuronNum).neuronName '-' char(allNeuronClassification{neuronNum})])
      export_fig(bpClassFigureFileName, '-pdf', '-append', gcf);
                delete(gcf)
end

 
allLabels = allNeuronClassification(recordingsToAnalyse);
countmember({'burst','pause','not sig'},allLabels)
vermLabels = allNeuronClassification(vermNeurons);
% countmember({'burst','pause','not sig'},vermLabels)
countmember({'burst','pause','not sig'},vermLabels(~cellfun(@isempty,vermLabels)))
latLabels = allNeuronClassification(latNeurons);
countmember({'burst','pause','not sig'},latLabels)
save(bpClassFileName,'allNeuronClassification')

else
    load(bpClassFileName)
    
end
%TODO now also do plot of eachmuSdf coloured by class and for each area.
%% SACCADE END POINT ERRORS
 alignmentTarget = 5; %colun to put peak in
postSaccCsTimeWindow = [0.05 0.3];


angVec = linspace(-pi,pi,9);
if doEndPoints
    
    conditionLegend = load('conditionLegend.mat');
    allNeuronsCsProbByAng = nan(max(recordingsToAnalyse),9);
%     allNeuronsCsProbByAng = nan(max(vermNeurons),9);
    
     
    for neuronNum = recordingsToAnalyse;
        %take all pro trials and get the saccade endpoint, compare to target location and get error direction.
%     keyboard
        thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
        theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        
        numTrials = length(theseCorProTrialNums);
        angsAndCs = nan(numTrials,3); %first column is error angle, second is 1 if cs present in 50-200ms post sacade, 3rd is amplitude of error
        for trialNum = 1:numTrials
            thisTrialNum = theseCorProTrialNums(trialNum);
            
            thisTrial = theseStableTrials(thisTrialNum);
            thisTrialCode =thisTrial.conditionCode;
            %I think we wil need the saccade origin as well, for now just
            %using overall amplitude.
            xTarg = conditionLegend.test{thisTrialCode,2};
            yTarg = conditionLegend.test{thisTrialCode,3};
            
            xEnd =  thisTrial.saccadeAmplitudeX+thisTrial.saccadeOriginX;
            yEnd = thisTrial.saccadeAmplitudeY+thisTrial.saccadeOriginY;
           
            
            %angle of error is calculated by
            %  erAng = atan((xEnd-xTarg)/(yEnd-yTarg));
            erAng = atan2((yTarg-yEnd),(xTarg-xEnd));
            erMag = sqrt((xEnd-xTarg)^2 +(yEnd-yTarg)^2);
            %now find if thee is a CS
            saccEndTime = thisTrial.saccadeTime+thisTrial.saccadeDuration/1000;
            
            thisTrialCs = thisTrial.alignedSpikes{thisNeuronUnits(2)};
            numPostSaccCs = sum(thisTrialCs>saccEndTime+0.05 & thisTrialCs<saccEndTime+0.2);
            angsAndCs(trialNum,2) = numPostSaccCs;
            angsAndCs(trialNum,1) =  erAng;
            angsAndCs(trialNum,3) =  erMag;
        end
        
        
        hF = figure;
        hA(1) = subaxis(4,1,1,'SV',0.15); %dist of err angs
        hA(2) = subaxis(4,1,2,'SV',0.15); %dist of err angs with a cs
        hA(3) = subaxis(4,1,3,'SV',0.15); %dist of prob of cs on ang
        hA(4) = subaxis(4,1,4,'SV',0.15); %dist of err mags with a cs
        
        bincounts = histc(angsAndCs(:,1),angVec);
        bar(angVec,bincounts,'histc','parent',hA(1));
        
        %now get only cs trials
        angsAndCsWithCs =  angsAndCs(angsAndCs(:,2)==1,:);
        bincountsCs = histc( angsAndCsWithCs(:,1),angVec);
        bar(angVec,bincountsCs,'histc','parent',hA(2));
        %now get prob of cs per angle
       %make a poalr version of this showing prob of cs in each error
        %direction.
        bincounts= bincounts(:);
        bincountsCs = bincountsCs(:);
        csProbByAng = bincountsCs./bincounts;
         bar(angVec,csProbByAng,'histc','parent',hA(3));
        magVec = [0:0.25:5];
        bincountsCsMag = histc( angsAndCsWithCs(:,3),magVec);
        bar(magVec,bincountsCsMag,'histc','parent',hA(4));
        
        set(get(hA(1),'title'),'string','All trials Error Angle')
        set(get(hA(2),'title'),'string','Post-sacc CS trials Error Angle')
        set(get(hA(1),'xlabel'),'string','Error Angle (rad)')
        set(get(hA(2),'xlabel'),'string','Error Angle (rad)')
        set(get(hA(3),'title'),'string','Prob of CS')
        set(get(hA(3),'xlabel'),'string','Error Angle (rad)')
        set(get(hA(4),'xlabel'),'string','Error Magnitude')
        set(get(hA(4),'title'),'string','Post-sacc CS trials Error Magnitude')
        suptitle(strrep(neuronList(neuronNum).neuronName,'_','-'))
        
        export_fig(hF,'-pdf','-append',erAngPdfName)
        delete(hF)
        
        
        hF = figure;
        
        polar([0 nan],[1 nan])
        hold on;
        polar(angVec,csProbByAng')
        title(strrep(neuronList(neuronNum).neuronName,'_','-'))
        export_fig(hF,'-pdf','-append',erAngPdfName)
        delete(hF)
        
        
        allNeuronsCsProbByAng(neuronNum,:) = csProbByAng;
    end
    
    
    %now get each neurons and fin cs-on direction and replot all tuning curves
    %aligned to this. Separate fro verm and lat
   
    hF = figure;
    allNeuronsCsProbByAngAligned =nan(max(recordingsToAnalyse),9);
    allNeuronsCsOndDir = nan(max(recordingsToAnalyse),1);
    polar([0 nan],[1 nan])
    hold on;
    
    
    for neuronNum = vermNeurons
        
        [peakVal peakDir] = max(allNeuronsCsProbByAng(neuronNum,:));
        if numel(peakDir)~=1
            peakDir = peakDir(1); %TODO should check neighbours
        end
        
        
        
        allNeuronsCsOndDir(neuronNum) = peakDir;
        
        %so we now want to do a circshift to move peakDir to position 4 and
        %all results as well
        alignedResults = circshift(allNeuronsCsProbByAng(neuronNum,:),[0 alignmentTarget-peakDir]);
        polar(angVec,alignedResults)
        
        allNeuronsCsProbByAngAligned(neuronNum,:) = alignedResults;
    end
    title('Verm neurons prob cs by error direction aligned to max')
    export_fig(hF,'-pdf','-append',erAngPdfName)
    delete(hF)
    
    
    hF = figure;
    %used to initialise ax lims to 1
    polar([0 nan],[1 nan])
    hold on;
     
    for neuronNum = latNeurons
        
        [peakVal peakDir] = max(allNeuronsCsProbByAng(neuronNum,:));
        if numel(peakDir)~=1
            peakDir = peakDir(1); %TODO should check neighbours
        end
        
        
        
        
        allNeuronsCsOndDir(neuronNum) = peakDir;
        %so we now want to do a circshift to move peakDir to position 4 and
        %all results as well
        alignedResults = circshift(allNeuronsCsProbByAng(neuronNum,:),[0 alignmentTarget-peakDir]);
        polar(angVec,alignedResults)
        
        allNeuronsCsProbByAngAligned(neuronNum,:) = alignedResults;
    end
    title('Lateral neurons prob cs by error direction aligned to max')
    export_fig(hF,'-pdf','-append',erAngPdfName)
    delete(hF)
    
    
    save([resultsDir filesep endPointErrorsFileName],'allNeuronsCsOndDir','allNeuronsCsProbByAngAligned','allNeuronsCsProbByAng')
end

%% fixation saccade aligned rasters
doFixSacRast = 0;
fixSacRastPdfNAme = [resultsDir filesep 'fixSaccadesOnlyGoodEyeTraces.pdf'];
if doFixSacRast
    
    thisUnitNum  =1; %CS are marked on rasters anyway
    for neuronNum = recordingsToAnalyse; %do 16 though end
        
        
        %build figure filename
        
        %loop over all chosen rasters and create and export their figures
         rastSelected = 'corSacTrials'
            
            %this raster is sorted yb go cue time as default
            [hAx1 hAx2 hA2] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
                wholeNeuronResults(neuronNum).selectedTrials.(char(rastSelected)),smallTimeWindow,timeWindow,...
                'alignTo','fixSacStartTimeRelToTrialStart')
            title([neuronList(neuronNum).neuronName ' ' rastSelected '-Fixation Saccades'])
            export_fig(gcf,'-pdf','-append',fixSacRastPdfNAme)
            
            delete(gcf)
%             if ismember(neuronNum,exampleNeurons)
%                 %export to eps
%                 thisFigSaveName = [resultsDir filesep neuronList(neuronNum).neuronName '-' rastSelected{:} '.eps'];
%                 export_fig(gcf,'-eps',thisFigSaveName)
%             end
%             delete(gcf)
%             
%             thisNeuronThisRastPngName =  [resultsDir filesep neuronList(neuronNum).neuronName '-' (char(rastSelected)) unsortedRasterPngName]
%             
%             %now make a raster with no sroting and just keeping the original
%             %trail order to see any cahgnes over time
%             [hAx1 hAx2 hA2] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
%                 wholeNeuronResults(neuronNum).selectedTrials.(char(rastSelected)),smallTimeWindow,timeWindow,'sortBy','');
%             
%             
%             title([neuronList(neuronNum).neuronName ' ' rastSelected])
%             export_fig(gcf,'-png',thisNeuronThisRastPngName)
%             delete(gcf)
        
        
    end
    
end
%% fixation saccade end points %TODO fix

            xTarg =0;
            yTarg =0;
if doEndPointsFix
    
    %conditionLegend = load('condtionLegend.mat');
    allNeuronsCsProbByAngFix = nan(max(recordingsToAnalyse),9);
    for neuronNum = recordingsToAnalyse;
        %take all pro trials and get the saccade endpoint, compare to target location and get error direction.
        thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
        theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        %TODO could be across all trials
        numTrials = length(theseCorProTrialNums);
        angsAndCs = nan(numTrials,3); %first column is error angle, second is 1 if cs present in 50-200ms post sacade, 3rd is amplitude of error
        for trialNum = 1:numTrials
            thisTrialNum = theseCorProTrialNums(trialNum);
            
            thisTrial = theseStableTrials(thisTrialNum);
          %  thisTrialCode =thisTrial.conditionCode;
            
            if ~isnan(thisTrial.fixSacAmplitudeX)
            xEnd =  thisTrial.fixSacAmplitudeX+thisTrial.fixSacOriginX;
            yEnd = thisTrial.fixSacAmplitudeY+thisTrial.fixSacOriginY;
            
            %angle of error is calculated by
            %  erAng = atan((xEnd-xTarg)/(yEnd-yTarg));
           erAng = atan2((yTarg-yEnd),(xTarg-xEnd));
            erMag = sqrt((xEnd-xTarg)^2 +(yEnd-yTarg)^2);
            %now find if thee is a CS
            %saccEndTime = thisTrial.saccadeTime+thisTrial.saccadeDuration/1000;
            saccEndTime = thisTrial.fixSacStartTimeRelToTrialStart+thisTrial.fixSacDuration/1000;
            
            thisTrialCs = thisTrial.alignedSpikes{thisNeuronUnits(2)};
            numPostSaccCs = sum(thisTrialCs>saccEndTime+postSaccCsTimeWindow(1) & thisTrialCs<saccEndTime+postSaccCsTimeWindow(2));
            angsAndCs(trialNum,2) = numPostSaccCs;
            angsAndCs(trialNum,1) =  erAng;
            angsAndCs(trialNum,3) =  erMag;
            end
        end
        
        
        hF = figure;
        hA(1) = subaxis(4,1,1,'SV',0.15); %dist of err angs
        hA(2) = subaxis(4,1,2,'SV',0.15); %dist of err angs with a cs
        hA(3) = subaxis(4,1,3,'SV',0.15); %dist cs prob
        hA(4) = subaxis(4,1,4,'SV',0.15); %dist of err mags with a cs
        bincounts = histc(angsAndCs(:,1),angVec);
        bar(angVec,bincounts,'histc','parent',hA(1));
        
        %now get only cs trials
        angsAndCsWithCs =  angsAndCs(angsAndCs(:,2)==1,:);
        bincountsCs = histc( angsAndCsWithCs(:,1),angVec);
        bar(angVec,bincountsCs,'histc','parent',hA(2));
        
        
        bincounts= bincounts(:);
        bincountsCs = bincountsCs(:);
        csProbByAng = bincountsCs./bincounts;
         bar(angVec,csProbByAng,'histc','parent',hA(3));
        magVec = [0:0.25:5];
        bincountsCsMag = histc( angsAndCsWithCs(:,3),magVec);
        bar(magVec,bincountsCsMag,'histc','parent',hA(4));
        
        set(get(hA(1),'title'),'string','All trials Error Angle')
        set(get(hA(2),'title'),'string','Post-sacc CS trials Error Angle')
        set(get(hA(1),'xlabel'),'string','Error Angle (rad)')
        set(get(hA(2),'xlabel'),'string','Error Angle (rad)')
        set(get(hA(3),'title'),'string','Prob of CS')
        set(get(hA(3),'xlabel'),'string','Error Angle (rad)')
        set(get(hA(4),'xlabel'),'string','Error Magnitude')
        set(get(hA(4),'title'),'string','Post-sacc CS trials Error Magnitude')
        suptitle(['Fixation saccades' strrep(neuronList(neuronNum).neuronName,'_','-')])
        
        export_fig(hF,'-pdf','-append',erAngFixationPdfName)
        delete(hF)
        
        %make a poalr version of this showing prob of cs in each error
        %direction.
        
        hF = figure;
        
        polar([0 nan],[1 nan])
        hold on;
        polar(angVec,csProbByAng')
        title(strrep(neuronList(neuronNum).neuronName,'_','-'))
        export_fig(hF,'-pdf','-append',erAngFixationPdfName)
        delete(hF)
        
        
        allNeuronsCsProbByAng(neuronNum,:) = csProbByAng;
    end
    
    
    %now get each neurons and fin cs-on direction and replot all tuning curves
    %aligned to this. Separate fro verm and lat
    alignmentTarget = 5; %colun to put peak in
    hF = figure;
    allNeuronsCsProbByAngAligned =nan(max(recordingsToAnalyse),9);
    allNeuronsCsOndDir = nan(max(recordingsToAnalyse),1);
    polar([0 nan],[1 nan])
    hold on;
    for neuronNum = vermNeurons
        
        [peakVal peakDir] = max(allNeuronsCsProbByAng(neuronNum,:));
        if numel(peakDir)~=1
            peakDir = peakDir(1); %TODO should check neighbours
        end
        
        
        
        allNeuronsCsOndDir(neuronNum) = peakDir;
        
        %so we now want to do a circshift to move peakDir to position 4 and
        %all results as well
        alignedResults = circshift(allNeuronsCsProbByAng(neuronNum,:),[0 alignmentTarget-peakDir]);
        polar(angVec,alignedResults)
        
        allNeuronsCsProbByAngAligned(neuronNum,:) = alignedResults;
    end
    title('Verm neurons prob cs by error direction aligned to max')
    export_fig(hF,'-pdf','-append',erAngFixationPdfName)
    delete(hF)
    
    
    hF = figure;
    %used to initialise ax lims to 1
    polar([0 nan],[1 nan])
    hold on;
    for neuronNum = latNeurons
        
        [peakVal peakDir] = max(allNeuronsCsProbByAng(neuronNum,:));
        if numel(peakDir)~=1
            peakDir = peakDir(1); %TODO should check neighbours
        end
        
        
        
        
        allNeuronsCsOndDir(neuronNum) = peakDir;
        %so we now want to do a circshift to move peakDir to position 4 and
        %all results as well
        alignedResults = circshift(allNeuronsCsProbByAng(neuronNum,:),[0 alignmentTarget-peakDir]);
        polar(angVec,alignedResults)
        
        allNeuronsCsProbByAngAligned(neuronNum,:) = alignedResults;
    end
    title('Lateral neurons prob cs by error direction aligned to max')
    export_fig(hF,'-pdf','-append',erAngFixationPdfName)
    delete(hF)
    
    
    save([resultsDir filesep endPointErrorsFixationFileName],'allNeuronsCsOndDir','allNeuronsCsProbByAngAligned','allNeuronsCsProbByAng')
end

%% POPULATION BY DIRECTION (should eventually rotate to CS-on)
%population analysis settings
numDirections = 8;
numKinBins = 50;
axisTicks = [3250 4000 4750];
axTickLabs = {'-0.75','0','0.75'};
expRend = '-painters'; %for rendering heatmaps
 sortByString = 'saccadePeakVel';
%     validSaccadeLimits.reactionTime = [0.1 0.35]; %s %[min max] allowed
%     validSaccadeLimits.saccadeDuration = [30 120]; %ms
%     validSaccadeLimits.saccadeAmplitude = [3 17]; %deg
%     validSaccadeLimits.saccadePeakVel = [100 800]; %deg/s
validSaccadeLimits.reactionTime = [0.05 0.6]; %s %[min max] allowed
    validSaccadeLimits.saccadeDuration = [10 200]; %ms
    validSaccadeLimits.saccadeAmplitude = [1 20]; %deg
    validSaccadeLimits.saccadePeakVel = [100 800]; %deg/s
normFrDispRange = [-30 30;0 10];
frDispRange = [0 200;0 5;0 200;0 200;0 5]; %simple spikes then complex spikes for graph lims

if doPopResByDir
    %load end point error data
    epErs = load([resultsDir filesep endPointErrorsFileName]);
    
    
    numDirections = 8;
    validDirectionRows = find(~isnan(trialConditionCodes(:,1)))';
    
   
    stInit = cell(max(recordingsToAnalyse),8); %need space for all 8 directions
    popResByDir = struct('alignedTrials',stInit,'antiMu',stInit,'proMu',stInit);
    for neuronNum = recordingsToAnalyse;
        hF = figure;
        
        
        thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
        unitNum  =  thisNeuronUnits(1);
        
        %get all correct trial numbers
        theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
        theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        
        theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
        theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        %now loop over directions and get mean response
        
        %now we want only trials in the specified diretion
        for directionNum = 1:numDirections
            %plNum is plotNum and we should run to 9
            plNum = validDirectionRows(directionNum);
            hA(directionNum) = subaxis(3,3,plNum,'SpacingVert',0.2,'MR',0.05);
            
            thisProCode = trialConditionCodes(validDirectionRows(directionNum),1);
            thisAntiCode = trialConditionCodes(validDirectionRows(directionNum),2);
            
            tempLogCell = cellfun(@(x) x==thisProCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirProTrialNums = find([tempLogCell{:}]);
            thisDirCorProTrialNums = intersect(theseProTrialNums,thisDirProTrialNums);
            
            tempLogCell = cellfun(@(x) x==thisAntiCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirAntiTrialNums = find([tempLogCell{:}]);
            thisDirCorAntiTrialNums = intersect(theseAntiTrialNums,thisDirAntiTrialNums);
            
            %align all trils to the saccade time
            [alignedTrials.Anti antiMu] = alignAndExtractRates(theseStableTrials,thisDirCorAntiTrialNums,'saccadeTime','unitNum',unitNum);
            [alignedTrials.Pro proMu] = alignAndExtractRates(theseStableTrials,thisDirCorProTrialNums,'saccadeTime','unitNum',unitNum);
            
            %now remove any trials which are outside our valid saccade
            %parameters
            [cleanAlignedTrials.Anti rejectionCodesAnti percentTrialsRemoved.Anti] = trimBadTrials(alignedTrials.Anti,validSaccadeLimits);
            [cleanAlignedTrials.Pro rejectionCodesPro percentTrialsRemoved.Pro] = trimBadTrials(alignedTrials.Pro,validSaccadeLimits);
            
            %now copy clean trials over top of original for use in rest of
            %script
            alignedTrials = cleanAlignedTrials;
            proSdfs = vertcat(alignedTrials.Pro.sdf);
            antiSdfs = vertcat(alignedTrials.Anti.sdf);
            if isempty(antiSdfs)
                 antiMu = nan(1,length(-timeWindow(1):0.001:timeWindow(2)));
            else
                
                antiMu = nanmean(antiSdfs,1);
            end
            
            if isempty(proSdfs)
                 proMu = nan(1,length(-timeWindow(1):0.001:timeWindow(2)));
            else
                 proMu = nanmean(proSdfs,1);
                
            end
           
            
            %now get the mean on sem of sdf for pro and anti and plot on
            %top of each other, then export as png for the summary figures
            
            [hF] = createConfLimSdf(alignedTrials,[0.5 0.5],'hA',hA(directionNum));
            if directionNum~=1
                set(get(hA(directionNum),'ylabel'),'string','')
            end
            
            popResByDir(neuronNum,directionNum).alignedTrials = alignedTrials;
            popResByDir(neuronNum,directionNum).proMu = proMu;
            popResByDir(neuronNum,directionNum).antiMu = antiMu;
        end
        
        %now plot the end point error cs tuning as spider in middle
        
        set(hA(1:7),'xticklabel','')
        suptitle(strrep(neuronList(neuronNum).neuronName,'_','-'))
        
        hAp = subaxis(3,3,5,'SpacingVert',0.2,'MR',0.05);
        set(hAp,'visible','off')
        
        hFt = figure;
        polar(angVec, epErs.allNeuronsCsProbByAng(neuronNum,:))
        
        hT = copyobj(gca,hF);
        set(hT,'units','normalized','position',[get(hAp,'position') + [0 -0.1 0 0.1]])
        delete(hFt)
        export_fig(hF,'-pdf','-append',popByDirFigName)
        delete(hF)
    end


%now create population sdf heatmaps for each direction, ffirst all neurons
%then lat and verm separately
for directionNum = 1:numDirections
    
    numNeurons = max(recordingsToAnalyse);
    tempProTrials = cell(numNeurons,1);
%     tempProTrials = cell(max([vermNeurons,latNeurons]),1);
    tempAntiTrials = cell(numNeurons,1);
%     tempAntiTrials = cell(max([vermNeurons,latNeurons]),1);

    for neuronNum = recordingsToAnalyse
        
        if ~isempty(popResByDir(neuronNum,directionNum).alignedTrials)
            tempProTrials{neuronNum} = popResByDir(neuronNum,directionNum).alignedTrials.Pro;
            
            tempAntiTrials{neuronNum} = popResByDir(neuronNum,directionNum).alignedTrials.Anti;
            
        end
        
    end
     sacVelrange = [350 450]; 
    allProTrials = [tempProTrials{:}];
    allAntiTrials = [tempAntiTrials{:}];
    [hA] = createSdfHeatMap(allProTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
    set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
   
        figure(get(hA,'parent'))
        suptitle(['All Neurons Pro, direction ' num2str(directionNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popByDirFigName)
        delete(gcf)
        [hA] = createSdfHeatMap(allAntiTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
       set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
        figure(get(hA,'parent'))
        suptitle(['All Neurons Anti, direction ' num2str(directionNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popByDirFigName)
        delete(gcf)
        
        %we also want the responses from each neuron in this direction,
        %coloured by burst/pause classification
        recsUsedClass = allNeuronClassification(recordingsToAnalyse); %get only the recordings used
        burstNeuronInds = find(strcmpi(recsUsedClass,'burst'));
        pauseNeuronInds = find(strcmpi(recsUsedClass,'pause'));
        nsNeuronInds = find(strcmpi(recsUsedClass,'not sig'));
        recsLoc = recList(4,recordingsToAnalyse);
        vermNeuronInds = find(strcmpi(recsLoc,'vermis'));
        latNeuronInds = find(strcmpi(recsLoc,'lateral'));
        burstNeuronIndsVerm = intersect(burstNeuronInds,vermNeuronInds);
        pauseNeuronIndsVerm = intersect(pauseNeuronInds,vermNeuronInds);
        nsNeuronIndsVerm = intersect(nsNeuronInds,vermNeuronInds);
        burstNeuronIndsLat = intersect(burstNeuronInds,latNeuronInds);
        pauseNeuronIndsLat = intersect(pauseNeuronInds,latNeuronInds);
        nsNeuronIndsLat = intersect(nsNeuronInds,latNeuronInds);
        
        hF = figure;
        hA(1) = subplot(2,2,1,'parent',hF);
        hA(2) = subplot(2,2,2,'parent',hF);
        hA(3) = subplot(2,2,3,'parent',hF);
        hA(4) = subplot(2,2,4,'parent',hF);
        tV  = -timeWindow(1):0.001:timeWindow(2);
        proMus = vertcat(popResByDir(recordingsToAnalyse,directionNum).proMu);
        antiMus = vertcat(popResByDir(recordingsToAnalyse,directionNum).antiMu);
        %these ned to be normalised soyou can see waht's going on
        proBaselineRates = mean(proMus(:,3800:3900),2);
        antiBaselineRates = mean(antiMus(:,3800:3900),2);
        normPros = bsxfun(@minus,proMus,proBaselineRates); %TODO switch to eric's max normalisation
        normAntis = bsxfun(@minus,antiMus,antiBaselineRates);
        
        if ~isempty(burstNeuronIndsVerm)
            line(tV,normPros(burstNeuronIndsVerm,:),'parent',hA(1),'color','b');
            line(tV,normAntis(burstNeuronIndsVerm,:),'parent',hA(2),'color','b');
        end
        if ~isempty(pauseNeuronIndsVerm)
            line(tV,normPros(pauseNeuronIndsVerm,:),'parent',hA(1),'color','r');
            line(tV,normAntis(pauseNeuronIndsVerm,:),'parent',hA(2),'color','r');
        end
        if ~isempty(nsNeuronIndsVerm)
            line(tV,normPros(nsNeuronIndsVerm,:),'parent',hA(1),'color','k');
            line(tV,normAntis(nsNeuronIndsVerm,:),'parent',hA(2),'color','k');
        end
        
        
        
        if ~isempty(burstNeuronIndsLat)
            line(tV,normPros(burstNeuronIndsLat,:),'parent',hA(3),'color','b');
            line(tV,normAntis(burstNeuronIndsLat,:),'parent',hA(4),'color','b');
        end
        if ~isempty(pauseNeuronIndsLat)
            line(tV,normPros(pauseNeuronIndsLat,:),'parent',hA(3),'color','r');
            line(tV,normAntis(pauseNeuronIndsLat,:),'parent',hA(4),'color','r');
        end
        if ~isempty(nsNeuronIndsLat)
            line(tV,normPros(nsNeuronIndsLat,:),'parent',hA(3),'color','k');
            line(tV,normAntis(nsNeuronIndsLat,:),'parent',hA(4),'color','k');
        end
        set(hA(:),'xlim',[-0.5 0.5])
        set(get(hA(1),'ylabel'),'string','vermis')
        set(get(hA(3),'ylabel'),'string','lateral')
        set(get(hA(3),'xlabel'),'string','pro')
        set(get(hA(4),'xlabel'),'string','anti')
        suptitle(['Direction: ' num2str(directionNum)])
        
         export_fig(gcf,'-pdf',expRend,'-append',popByDirFigName)
        delete(gcf)
        
        
        tempProTrials = cell(numNeurons,1);
    tempAntiTrials = cell(numNeurons,1);
        %now we want to do the same with only saccades in a small parameter
        %window
       
        %for all neurons remove any trials outside this range
        for neuronNum = vermNeurons
            
            if ~isempty(popResByDir(neuronNum,directionNum).alignedTrials)
                tempTrials =  popResByDir(neuronNum,directionNum).alignedTrials.Pro;
                inRangeTrials = find([tempTrials.saccadePeakVel]>sacVelrange(1) &  [tempTrials.saccadePeakVel]<sacVelrange(2));
                tempProTrials{neuronNum} = tempTrials(inRangeTrials);
                
                tempTrials = popResByDir(neuronNum,directionNum).alignedTrials.Anti;
                inRangeTrials = find([tempTrials.saccadePeakVel]>sacVelrange(1) &  [tempTrials.saccadePeakVel]<sacVelrange(2));
                
                tempAntiTrials{neuronNum} = tempTrials(inRangeTrials);
                
            end
            
        end
        
         selProTrials = [tempProTrials{:}];
    selAntiTrials = [tempAntiTrials{:}];
    if ~isempty(selProTrials)
    %now we want for this direction to get the population response, to do
    %this we do a ksdensity with kernel of 2.5ms
    proSpikes = vertcat(selProTrials.alignedSpikes);
    proSpikes = proSpikes(:,1); %get only ss
    proSpikes = vertcat(proSpikes{:});
    proPopSdf = ksdensity(proSpikes,tV,'width',0.0025)*length(proSpikes)/length(selProTrials);
    else
        proPopSdf = nan(size(tV));
        
    end
     if ~isempty(selAntiTrials)
     antiSpikes = vertcat(selAntiTrials.alignedSpikes);
    antiSpikes = antiSpikes(:,1); %get only ss
    antiSpikes = vertcat(antiSpikes{:});
    antiPopSdf = ksdensity(antiSpikes,tV,'width',0.0025)*length(antiSpikes)/length(selAntiTrials);
     else
         antiPopSdf = nan(size(tV));
     end
    hF = figure;
    hA = axes('parent',hF); 
    line(tV,proPopSdf,'parent',hA,'color','g')
    line(tV,antiPopSdf,'parent',hA,'color','r')
    xlim([-0.5 0.5])
    title(['Vermis. Direction ' num2str(directionNum) ', numPro:' num2str(size(selProTrials,2)) ', numAnti:' num2str(size(selAntiTrials,2)) ', Sac Vel:' num2str(sacVelrange(1)) ':' num2str(sacVelrange(2))])
     export_fig(gcf,'-pdf',expRend,'-append',popByDirFigName)
        delete(gcf)
        
        
        %now exactly the same for lateral neurons
            tempProTrials = cell(numNeurons,1);
    tempAntiTrials = cell(numNeurons,1);
        %now we want to do the same with only saccades in a small parameter
        %window
        sacVelrange = [350 450]; 
        %for all neurons remove any trials outside this range
        for neuronNum = latNeurons
            
            if ~isempty(popResByDir(neuronNum,directionNum).alignedTrials)
                tempTrials =  popResByDir(neuronNum,directionNum).alignedTrials.Pro;
                inRangeTrials = find([tempTrials.saccadePeakVel]>sacVelrange(1) &  [tempTrials.saccadePeakVel]<sacVelrange(2));
                tempProTrials{neuronNum} = tempTrials(inRangeTrials);
                
                tempTrials = popResByDir(neuronNum,directionNum).alignedTrials.Anti;
                inRangeTrials = find([tempTrials.saccadePeakVel]>sacVelrange(1) &  [tempTrials.saccadePeakVel]<sacVelrange(2));
                
                tempAntiTrials{neuronNum} = tempTrials(inRangeTrials);
                
            end
            
        end
        
         selProTrials = [tempProTrials{:}];
    selAntiTrials = [tempAntiTrials{:}];
    if ~isempty(selProTrials)
    %now we want for this direction to get the population response, to do
    %this we do a ksdensity with kernel of 2.5ms
    proSpikes = vertcat(selProTrials.alignedSpikes);
    proSpikes = proSpikes(:,1); %get only ss
    proSpikes = vertcat(proSpikes{:});
    proPopSdf = ksdensity(proSpikes,tV,'width',0.0025)*length(proSpikes)/length(selProTrials);
    else
        proPopSdf = nan(size(tV));
        
    end
     if ~isempty(selAntiTrials)
     antiSpikes = vertcat(selAntiTrials.alignedSpikes);
    antiSpikes = antiSpikes(:,1); %get only ss
    antiSpikes = vertcat(antiSpikes{:});
    antiPopSdf = ksdensity(antiSpikes,tV,'width',0.0025)*length(antiSpikes)/length(selAntiTrials);
     else
         antiPopSdf = nan(size(tV));
     end
    hF = figure;
    hA = axes('parent',hF); 
    line(tV,proPopSdf,'parent',hA,'color','g')
    line(tV,antiPopSdf,'parent',hA,'color','r')
    xlim([-0.5 0.5])
    title(['Lateral. Direction ' num2str(directionNum) ', numPro:' num2str(size(selProTrials,2)) ', numAnti:' num2str(size(selAntiTrials,2)) ', Sac Vel:' num2str(sacVelrange(1)) ':' num2str(sacVelrange(2))])
     export_fig(gcf,'-pdf',expRend,'-append',popByDirFigName)
        delete(gcf)
end

end
%% POP RES BY DIR BUT ROTATED TO CS-on
 %alignemnt isn't simple circ shift, first set up a matrix with plot
        %numbers and the corresponsing angle.
        plAngCode = [1 0.75*pi; 2 0.5*pi; 3 0.25*pi; 4 pi; 5 0; 6 -0.75*pi; 7 -0.5*pi; 8 -0.25*pi];
if doPopResByDir
    %load end point error data
    epErs = load([resultsDir filesep endPointErrorsFileName]);
    
    
    numDirections = 8;
    validDirectionRows = find(~isnan(trialConditionCodes(:,1)))';
    
   
    stInit = cell(max(recordingsToAnalyse),8); %need space for all 8 directions
    popResByDirRot = struct('alignedTrials',stInit,'antiMu',stInit,'proMu',stInit);
    for neuronNum = recordingsToAnalyse;
        hF = figure;
        
        
        thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
        unitNum  =  thisNeuronUnits(1);
        
        %get all correct trial numbers
        theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
        theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        
        theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
        theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        %now loop over directions and get mean response
        %get this neurons cs on
        peakDir = epErs.allNeuronsCsOndDir(neuronNum);
         alignedCs = circshift(epErs.allNeuronsCsProbByAng(neuronNum,:),[0 alignmentTarget-peakDir]);
       % thisNeuronCsAlignedDirs = circshift(1:numDirections,[0 -1*(alignmentTarget-peakDir)]); %opposite direction of switch as polar runs counterlcok this runs clock
        
       
        %now find the angle changed by cs alignemnt
        
        alignmentCorrection = -1*(alignmentTarget-peakDir)*pi/4;
        thisAngCode = plAngCode(:,2)- alignmentCorrection;
        thisAngCode = wrapToPi(thisAngCode);
        thisAngCode(thisAngCode==pi)=-pi;
        thisNeuronCsAlignedDirs = nan(8,1);
        for directionNum = 1:numDirections
            [val thisNeuronCsAlignedDirs(directionNum)] = min(abs(thisAngCode-plAngCode(directionNum,2)));
        end
        %now we want only trials in the specified diretion
        for directionNum = 1:numDirections
            %plNum is plotNum and we should run to 9
            plNum = validDirectionRows(directionNum);
            hA(directionNum) = subaxis(3,3,plNum,'SpacingVert',0.2,'MR',0.05);
            %updated to get the correct codes for aligned directions
            thisProCode = trialConditionCodes(validDirectionRows(thisNeuronCsAlignedDirs(directionNum)),1);
            thisAntiCode = trialConditionCodes(validDirectionRows(thisNeuronCsAlignedDirs(directionNum)),2);
            
            tempLogCell = cellfun(@(x) x==thisProCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirProTrialNums = find([tempLogCell{:}]);
            thisDirCorProTrialNums = intersect(theseProTrialNums,thisDirProTrialNums);
            
            tempLogCell = cellfun(@(x) x==thisAntiCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirAntiTrialNums = find([tempLogCell{:}]);
            thisDirCorAntiTrialNums = intersect(theseAntiTrialNums,thisDirAntiTrialNums);
            
            %align all trils to the saccade time
            [alignedTrials.Anti antiMu] = alignAndExtractRates(theseStableTrials,thisDirCorAntiTrialNums,'saccadeTime','unitNum',unitNum);
            [alignedTrials.Pro proMu] = alignAndExtractRates(theseStableTrials,thisDirCorProTrialNums,'saccadeTime','unitNum',unitNum);
            
            %now remove any trials which are outside our valid saccade
            %parameters
            [cleanAlignedTrials.Anti rejectionCodesAnti percentTrialsRemoved.Anti] = trimBadTrials(alignedTrials.Anti,validSaccadeLimits);
            [cleanAlignedTrials.Pro rejectionCodesPro percentTrialsRemoved.Pro] = trimBadTrials(alignedTrials.Pro,validSaccadeLimits);
            
            %now copy clean trials over top of original for use in rest of
            %script
            alignedTrials = cleanAlignedTrials;
            proSdfs = vertcat(alignedTrials.Pro.sdf);
            antiSdfs = vertcat(alignedTrials.Anti.sdf);
            if isempty(antiSdfs)
                 antiMu = nan(1,length(-timeWindow(1):0.001:timeWindow(2)));
            else
                
                antiMu = nanmean(antiSdfs,1);
            end
            
            if isempty(proSdfs)
                 proMu = nan(1,length(-timeWindow(1):0.001:timeWindow(2)));
            else
                 proMu = nanmean(proSdfs,1);
                
            end
           
            
            %now get the mean on sem of sdf for pro and anti and plot on
            %top of each other, then export as png for the summary figures
            
            [hF] = createConfLimSdf(alignedTrials,[0.5 0.5],'hA',hA(directionNum));
            if directionNum~=1
                set(get(hA(directionNum),'ylabel'),'string','')
            end
            
            popResByDirRot(neuronNum,directionNum).alignedTrials = alignedTrials;
            popResByDirRot(neuronNum,directionNum).proMu = proMu;
            popResByDirRot(neuronNum,directionNum).antiMu = antiMu;
        end
        
        %now plot the end point error cs tuning as spider in middle
        
        set(hA(1:7),'xticklabel','')
        suptitle(strrep(neuronList(neuronNum).neuronName,'_','-'))
        
        hAp = subaxis(3,3,5,'SpacingVert',0.2,'MR',0.05);
        set(hAp,'visible','off')
        
        hFt = figure;
        polar(angVec,alignedCs)
        
        hT = copyobj(gca,hF);
        set(hT,'units','normalized','position',[get(hAp,'position') + [0 -0.1 0 0.1]])
        delete(hFt)
        export_fig(hF,'-pdf','-append',popByDirFigNameRot)
        delete(hF)
    end


%now create population sdf heatmaps for each direction, ffirst all neurons
%then lat and verm separately
for directionNum = 1:numDirections
    
    numNeurons = max(recordingsToAnalyse);
    tempProTrials = cell(numNeurons,1);
    tempAntiTrials = cell(numNeurons,1);
    
    for neuronNum = recordingsToAnalyse
        
        if ~isempty(popResByDir(neuronNum,directionNum).alignedTrials)
            tempProTrials{neuronNum} = popResByDir(neuronNum,directionNum).alignedTrials.Pro;
            
            tempAntiTrials{neuronNum} = popResByDir(neuronNum,directionNum).alignedTrials.Anti;
            
        end
        
    end
     sacVelrange = [350 450]; 
    allProTrials = [tempProTrials{:}];
    allAntiTrials = [tempAntiTrials{:}];
    [hA] = createSdfHeatMap(allProTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
    set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
   
        figure(get(hA,'parent'))
        suptitle(['All Neurons Pro, direction ' num2str(directionNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popByDirFigNameRot)
        delete(gcf)
        [hA] = createSdfHeatMap(allAntiTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
       set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
        figure(get(hA,'parent'))
        suptitle(['All Neurons Anti, direction ' num2str(directionNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popByDirFigNameRot)
        delete(gcf)
        
        %we also want the responses from each neuron in this direction,
        %coloured by burst/pause classification
        recsUsedClass = allNeuronClassification(recordingsToAnalyse); %get only the recordings used
        burstNeuronInds = find(strcmpi(recsUsedClass,'burst'));
        pauseNeuronInds = find(strcmpi(recsUsedClass,'pause'));
        nsNeuronInds = find(strcmpi(recsUsedClass,'not sig'));
        recsLoc = recList(4,recordingsToAnalyse);
        vermNeuronInds = find(strcmpi(recsLoc,'vermis'));
        latNeuronInds = find(strcmpi(recsLoc,'lateral'));
        burstNeuronIndsVerm = intersect(burstNeuronInds,vermNeuronInds);
        pauseNeuronIndsVerm = intersect(pauseNeuronInds,vermNeuronInds);
        nsNeuronIndsVerm = intersect(nsNeuronInds,vermNeuronInds);
        burstNeuronIndsLat = intersect(burstNeuronInds,latNeuronInds);
        pauseNeuronIndsLat = intersect(pauseNeuronInds,latNeuronInds);
        nsNeuronIndsLat = intersect(nsNeuronInds,latNeuronInds);
        
        hF = figure;
        hA(1) = subplot(2,2,1,'parent',hF);
        hA(2) = subplot(2,2,2,'parent',hF);
        hA(3) = subplot(2,2,3,'parent',hF);
        hA(4) = subplot(2,2,4,'parent',hF);
        tV  = -timeWindow(1):0.001:timeWindow(2);
        proMus = vertcat(popResByDir(recordingsToAnalyse,directionNum).proMu);
        antiMus = vertcat(popResByDir(recordingsToAnalyse,directionNum).antiMu);
        %these ned to be normalised soyou can see waht's going on
        proBaselineRates = mean(proMus(:,3800:3900),2);
        antiBaselineRates = mean(antiMus(:,3800:3900),2);
        normPros = bsxfun(@minus,proMus,proBaselineRates); %TODO switch to eric's max normalisation
        normAntis = bsxfun(@minus,antiMus,antiBaselineRates);
        
        if ~isempty(burstNeuronIndsVerm)
            line(tV,normPros(burstNeuronIndsVerm,:),'parent',hA(1),'color','b');
            line(tV,normAntis(burstNeuronIndsVerm,:),'parent',hA(2),'color','b');
        end
        if ~isempty(pauseNeuronIndsVerm)
            line(tV,normPros(pauseNeuronIndsVerm,:),'parent',hA(1),'color','r');
            line(tV,normAntis(pauseNeuronIndsVerm,:),'parent',hA(2),'color','r');
        end
        if ~isempty(nsNeuronIndsVerm)
            line(tV,normPros(nsNeuronIndsVerm,:),'parent',hA(1),'color','k');
            line(tV,normAntis(nsNeuronIndsVerm,:),'parent',hA(2),'color','k');
        end
        
        
        
        if ~isempty(burstNeuronIndsLat)
            line(tV,normPros(burstNeuronIndsLat,:),'parent',hA(3),'color','b');
            line(tV,normAntis(burstNeuronIndsLat,:),'parent',hA(4),'color','b');
        end
        if ~isempty(pauseNeuronIndsLat)
            line(tV,normPros(pauseNeuronIndsLat,:),'parent',hA(3),'color','r');
            line(tV,normAntis(pauseNeuronIndsLat,:),'parent',hA(4),'color','r');
        end
        if ~isempty(nsNeuronIndsLat)
            line(tV,normPros(nsNeuronIndsLat,:),'parent',hA(3),'color','k');
            line(tV,normAntis(nsNeuronIndsLat,:),'parent',hA(4),'color','k');
        end
        set(hA(:),'xlim',[-0.5 0.5])
        set(get(hA(1),'ylabel'),'string','vermis')
        set(get(hA(3),'ylabel'),'string','lateral')
        set(get(hA(3),'xlabel'),'string','pro')
        set(get(hA(4),'xlabel'),'string','anti')
        suptitle(['Direction: ' num2str(directionNum)])
        
         export_fig(gcf,'-pdf',expRend,'-append',popByDirFigNameRot)
        delete(gcf)
        
        
        tempProTrials = cell(numNeurons,1);
    tempAntiTrials = cell(numNeurons,1);
        %now we want to do the same with only saccades in a small parameter
        %window
       
        %for all neurons remove any trials outside this range
        for neuronNum = vermNeurons
            
            if ~isempty(popResByDir(neuronNum,directionNum).alignedTrials)
                tempTrials =  popResByDir(neuronNum,directionNum).alignedTrials.Pro;
                inRangeTrials = find([tempTrials.saccadePeakVel]>sacVelrange(1) &  [tempTrials.saccadePeakVel]<sacVelrange(2));
                tempProTrials{neuronNum} = tempTrials(inRangeTrials);
                
                tempTrials = popResByDir(neuronNum,directionNum).alignedTrials.Anti;
                inRangeTrials = find([tempTrials.saccadePeakVel]>sacVelrange(1) &  [tempTrials.saccadePeakVel]<sacVelrange(2));
                
                tempAntiTrials{neuronNum} = tempTrials(inRangeTrials);
                
            end
            
        end
        
         selProTrials = [tempProTrials{:}];
    selAntiTrials = [tempAntiTrials{:}];
    if ~isempty(selProTrials)
    %now we want for this direction to get the population response, to do
    %this we do a ksdensity with kernel of 2.5ms
    proSpikes = vertcat(selProTrials.alignedSpikes);
    proSpikes = proSpikes(:,1); %get only ss
    proSpikes = vertcat(proSpikes{:});
    proPopSdf = ksdensity(proSpikes,tV,'width',0.0025)*length(proSpikes)/length(selProTrials);
    else
        proPopSdf = nan(size(tV));
        
    end
     if ~isempty(selAntiTrials)
     antiSpikes = vertcat(selAntiTrials.alignedSpikes);
    antiSpikes = antiSpikes(:,1); %get only ss
    antiSpikes = vertcat(antiSpikes{:});
    antiPopSdf = ksdensity(antiSpikes,tV,'width',0.0025)*length(antiSpikes)/length(selAntiTrials);
     else
         antiPopSdf = nan(size(tV));
     end
    hF = figure;
    hA = axes('parent',hF); 
    line(tV,proPopSdf,'parent',hA,'color','g')
    line(tV,antiPopSdf,'parent',hA,'color','r')
    xlim([-0.5 0.5])
    title(['Vermis. Direction ' num2str(directionNum) ', numPro:' num2str(size(selProTrials,2)) ', numAnti:' num2str(size(selAntiTrials,2)) ', Sac Vel:' num2str(sacVelrange(1)) ':' num2str(sacVelrange(2))])
     export_fig(gcf,'-pdf',expRend,'-append',popByDirFigNameRot)
        delete(gcf)
        
        
        %now exactly the same for lateral neurons
            tempProTrials = cell(numNeurons,1);
    tempAntiTrials = cell(numNeurons,1);
        %now we want to do the same with only saccades in a small parameter
        %window
        sacVelrange = [350 450]; 
        %for all neurons remove any trials outside this range
        for neuronNum = latNeurons
            
            if ~isempty(popResByDir(neuronNum,directionNum).alignedTrials)
                tempTrials =  popResByDir(neuronNum,directionNum).alignedTrials.Pro;
                inRangeTrials = find([tempTrials.saccadePeakVel]>sacVelrange(1) &  [tempTrials.saccadePeakVel]<sacVelrange(2));
                tempProTrials{neuronNum} = tempTrials(inRangeTrials);
                
                tempTrials = popResByDir(neuronNum,directionNum).alignedTrials.Anti;
                inRangeTrials = find([tempTrials.saccadePeakVel]>sacVelrange(1) &  [tempTrials.saccadePeakVel]<sacVelrange(2));
                
                tempAntiTrials{neuronNum} = tempTrials(inRangeTrials);
                
            end
            
        end
        
         selProTrials = [tempProTrials{:}];
    selAntiTrials = [tempAntiTrials{:}];
    if ~isempty(selProTrials)
    %now we want for this direction to get the population response, to do
    %this we do a ksdensity with kernel of 2.5ms
    proSpikes = vertcat(selProTrials.alignedSpikes);
    proSpikes = proSpikes(:,1); %get only ss
    proSpikes = vertcat(proSpikes{:});
    proPopSdf = ksdensity(proSpikes,tV,'width',0.0025)*length(proSpikes)/length(selProTrials);
    else
        proPopSdf = nan(size(tV));
        
    end
     if ~isempty(selAntiTrials)
     antiSpikes = vertcat(selAntiTrials.alignedSpikes);
    antiSpikes = antiSpikes(:,1); %get only ss
    antiSpikes = vertcat(antiSpikes{:});
    antiPopSdf = ksdensity(antiSpikes,tV,'width',0.0025)*length(antiSpikes)/length(selAntiTrials);
     else
         antiPopSdf = nan(size(tV));
     end
    hF = figure;
    hA = axes('parent',hF); 
    line(tV,proPopSdf,'parent',hA,'color','g')
    line(tV,antiPopSdf,'parent',hA,'color','r')
    xlim([-0.5 0.5])
    title(['Lateral. Direction ' num2str(directionNum) ', numPro:' num2str(size(selProTrials,2)) ', numAnti:' num2str(size(selAntiTrials,2)) ', Sac Vel:' num2str(sacVelrange(1)) ':' num2str(sacVelrange(2))])
     export_fig(gcf,'-pdf',expRend,'-append',popByDirFigNameRot)
        delete(gcf)
end
end
%% CORRELATION TO KINEMATICS
%the idea is for every neuron to get a time series which is the correlation
%between all trials (pro only) and the kinematics. For example if saccade
%velocity then look in a 5ms bin in all trials aligned to saccade onset and
%then get r and p value for that bin, then move the bin in time and repeat.
if doCorKin

corSigLevel = 0.01;
saccadeInd = 4000;
binSize = 0.005;
binSizeSamples = binSize*1000;
binningLims = [0.5 0.5]; %pre and post saccade time in s
binningLimsSamples = [saccadeInd-binningLims(1)*1000 saccadeInd+binningLims(2)*1000]; %in trial relative samples
numBins = sum(binningLims)/binSize;
baselineInds = [2750 3250];
%{'saccadePeakVel','saccadeAmplitude','reactionTime','saccadeDuration','saccadeAmplitudeX','saccadeAmplitudeY'};
for sacKin = {'saccadeAmplitudeX','saccadeAmplitudeY'};
allNeuronsKinCors = nan(max(recordingsToAnalyse),numBins);
allNeuronsKinCorPvals = nan(max(recordingsToAnalyse),numBins);
for neuronNum = recordingsToAnalyse;
    
    thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
    unitNum = thisNeuronUnits(1);
    theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
    
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
    
    [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseCorProTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
    theseSdfs = vertcat(alignedTrials.Pro.sdf);
     theseBaseRates = mean(theseSdfs(:,baselineInds(1):baselineInds(2)),2);
    %normalise to baseline trial by trial
    %theseNormSdfs = bsxfun(@ldivide,theseBaseRates,theseSdfs); %other norm
    theseNormSdfs = bsxfun(@minus,theseSdfs,theseBaseRates); 
    theseKinematics = [alignedTrials.Pro.(char(sacKin))];
    theseKinematics = theseKinematics(:);
    numTrials = length(theseKinematics);
    theseCorrResultsPvals = nan(numBins,1);
    theseCorrResults = nan(numBins,1);
     for xBinNum = 1:numBins
         xBinInds = [binningLimsSamples(1)+binSizeSamples*(xBinNum-1)+1 binningLimsSamples(1)+binSizeSamples*(xBinNum)];
      
         
          theseSdfValues = nanmean(theseNormSdfs(:,xBinInds(1):xBinInds(2)),2);
          
          %do gibbs outlier removal 
        x_date = cellfun(@(x) num2str(x),num2cell(1:numTrials),'uniformoutput',false);
        badTrialsNanSdf = find(isnan(theseSdfValues));
        if ~isempty(badTrialsNanSdf)
            goodTrialsNanSdf = setdiff(1:numTrials,badTrialsNanSdf);
        else
             goodTrialsNanSdf = 1:numTrials;
        end
        [Gmax_loop Gmin_loop G2_loop lambda1_loop lambda2_loop outlier outlier_num] = grubbs(theseSdfValues, x_date');
       if ~isempty(outlier_num)
           badTrials = outlier_num(:,1);
       else
           badTrials = [];
       end
       goodTrialsOut = setdiff(1:numTrials,badTrials);
       goodTrials  = intersect(goodTrialsNanSdf,goodTrialsOut);
       trialsToUse = goodTrials;
        [theseCorrResults(xBinNum) theseCorrResultsPvals(xBinNum)] = corr(theseSdfValues(trialsToUse) ,theseKinematics(trialsToUse));
   
     end
     timeVec = linspace(-binningLims(1),binningLims(2), numBins);
     hF = figure;
     hA(1) = axes('parent',hF)
     hA(2) = axes('parent',hF,'position',get(hA(1),'position'),'yaxislocation','right','color','none','ycolor','r')
     line(timeVec,theseCorrResults,'parent',hA(1),'color','b')
     line(timeVec,theseCorrResultsPvals,'parent',hA(2),'color','r')
     line([timeVec(1) timeVec(numBins)],[corSigLevel corSigLevel],'parent',hA(2),'color','r','linestyle','-.')
     set(get(hA(1),'ylabel'),'string','Rho')
     set(get(hA(2),'ylabel'),'string','p')
     set(get(hA(1),'xlabel'),'string','time to saccade onset (s)')
     suptitle([neuronList(neuronNum).neuronName '-' (char(sacKin))])
export_fig([resultsDir filesep corKinFigNameStart '-' char(sacKin) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
    allNeuronsKinCors(neuronNum,:) = theseCorrResults;
    allNeuronsKinCorPvals(neuronNum,:) = theseCorrResultsPvals;
end

    [val inds] = sort(epErs.allNeuronsCsOndDir(recordingsToAnalyse));
%now a summary figure showing a heatplot of all
hF = figure;
hA = axes('parent',hF);
colormap(jet)
theseCorrResults = allNeuronsKinCors(recordingsToAnalyse,:);
imagesc(theseCorrResults,[-0.5 0.5])
set(get(hA,'ylabel'),'string','Neuron Number')
set(get(hA,'xlabel'),'string','time to saccade onset (s)')
line([1 numBins],[length(vermNeurons) length(vermNeurons)]+0.5,'parent',hA,'color','k','linewidth',2)
set(hA,'xtick',[1 numBins/2 numBins],'xticklabel',{'-0.5','0','0.5'})
hC = colorbar;
set(get(hC,'title'),'string','R')

export_fig([resultsDir filesep corKinFigNameStart '-' char(sacKin) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    %now exactly the same thing but arrange neurons by CS-on Dir
   
    theseSortedCorrResults = theseCorrResults(inds,:);
    hF = figure;
colormap(jet)
hA = axes('parent',hF);

imagesc(theseSortedCorrResults,[-0.5 0.5])
set(get(hA,'ylabel'),'string','Neuron CS ON dir')
set(get(hA,'xlabel'),'string','time to saccade onset (s)')
line([1 numBins],[length(vermNeurons) length(vermNeurons)]+0.5,'parent',hA,'color','k','linewidth',2)
set(hA,'xtick',[1 numBins/2 numBins],'xticklabel',{'-0.5','0','0.5'})
hC = colorbar;
set(get(hC,'title'),'string','p')
export_fig([resultsDir filesep corKinFigNameStart '-' char(sacKin) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
  %now thresholded p vals
hF = figure;
colormap([1 1 1; 0 0 0])  
hA = axes('parent',hF);
theseCorrResults = allNeuronsKinCorPvals(recordingsToAnalyse,:);
imagesc(theseCorrResults,[0 0.05])
set(get(hA,'ylabel'),'string','Neuron Number')
set(get(hA,'xlabel'),'string','time to saccade onset (s)')
line([1 numBins],[length(vermNeurons) length(vermNeurons)]+0.5,'parent',hA,'color','k','linewidth',2)
set(hA,'xtick',[1 numBins/2 numBins],'xticklabel',{'-0.5','0','0.5'})
hC = colorbar;
set(get(hC,'title'),'string','p')
export_fig([resultsDir filesep corKinFigNameStart '-' char(sacKin) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
       %now exactly the same thing but arrange neurons by CS-on Dir
    
    
    theseSortedCorrResults = theseCorrResults(inds,:);
    hF = figure;
colormap([1 1 1; 0 0 0])  
hA = axes('parent',hF);

imagesc(theseSortedCorrResults,[0 0.05])
set(get(hA,'ylabel'),'string','Neuron CS ON dir')
set(get(hA,'xlabel'),'string','time to saccade onset (s)')
line([1 numBins],[length(vermNeurons) length(vermNeurons)]+0.5,'parent',hA,'color','k','linewidth',2)
set(hA,'xtick',[1 numBins/2 numBins],'xticklabel',{'-0.5','0','0.5'})
hC = colorbar;
set(get(hC,'title'),'string','p')
export_fig([resultsDir filesep corKinFigNameStart '-' char(sacKin) '.pdf'], '-pdf','-append', hF);
    delete(hF)
end

end





%% CS TRIGGERED EYE TRACES WHOLENEURON
%should get the complex spike pause for all files concatenated.

if wholeNeurCsSta

for neuronNum = recordingsToAnalyse;
   
    %load matfile into memory
     thisMatFileName = [processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd];
     mLoadedData = matfile(thisMatFileName); %create mat file
   %build filename for analysis, not needed for CS pause
        thisAnalysisFileName = [resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd];
        
        thisFigFileName = [resultsDir filesep neuronList(neuronNum).neuronName csStaFigNameEnd];
        
        numFiles = size(mLoadedData.editedFileList,1);
        load(thisAnalysisFileName);
      %  stablePeriodTrials = cell(1,numFiles);
        
        thisNeuronAlignedStas = cell(numFiles,3); %x first colummy second, 3rd linV
        %loop over files and delete trials falling otuside marked times
        for fileNum = 1:numFiles
            
            %get file info
            thisFileInfo = mLoadedData.editedFileList(fileNum,1);
            
            %get start and end of stable section in this file
            if isinf(thisFileInfo.startTime)
                thisStartTime = -inf;
            else
                thisStartTime = thisFileInfo.startTime;
            end
            thisStopTime = thisFileInfo.stopTime;
            
            
            thisFileStableSpikes = cellfun(@(x) x(x>thisStartTime & x<thisStopTime),analysisResults(fileNum).alignedSpikeTimesCell,'uniformoutput',false);
            calData = mLoadedData.calibEye(fileNum,1); %get eye data
            saccFileName = [processDir filesep thisFileInfo.fileName(1:end-4) '-saccDetectNewest.mat'];
           load(saccFileName)
           
            %now we need to get the eye trace triggerd off every cs
            if length(thisFileInfo.unitNumbers)>1
                %hF = figure;
%                 [hRast thisNeuronAlignedRasters{fileNum}] = createRaster(thisFileStableSpikes{thisFileInfo.unitNumbers(1)},...
%                     thisFileStableSpikes{thisFileInfo.unitNumbers(2)},[0.1 0.1],'spikeMarkSize',3); %TODO also make this plot other complex spikes
%                  title(thisFileInfo.fileName(1:end-4))
%                 %export_fig(figureFileName, '-pdf', '-append', gcf);
%                 delete(gcf)
                %TODO could get output from Marios to get velocty and cleaned traces etc.
%                 thisNeuronAlignedStas{fileNum,1} = extractEpochs(calData.eyeXcal,thisFileStableSpikes{thisFileInfo.unitNumbers(2)},[4 4],calData.eyeFsCal,'displayFlag',1);
%             set(get(gca,'title'),'string',[thisFileInfo.fileName 'X pos'])
%             xlim([-1 1])
%             export_fig(gcf,figureFileName,'-append','-pdf')
%             delete(gcf)
             if~ isempty(thisFileStableSpikes{thisFileInfo.unitNumbers(2)})
             thisNeuronAlignedStas{fileNum,1} = extractEpochs(saccadeDetectOut.traces.position.raw(:,1),thisFileStableSpikes{thisFileInfo.unitNumbers(2)},[4 4],calData.eyeFsCal,'displayFlag',1);
            set(get(gca,'title'),'string',[thisFileInfo.fileName 'X pos  mar'])
            xlim([-1 1])
            export_fig(gcf,thisFigFileName,'-append','-pdf')
            delete(gcf)
               
             thisNeuronAlignedStas{fileNum,2} = extractEpochs(saccadeDetectOut.traces.position.raw(:,2),thisFileStableSpikes{thisFileInfo.unitNumbers(2)},[4 4],calData.eyeFsCal,'displayFlag',1);
            set(get(gca,'title'),'string',[thisFileInfo.fileName 'Y pos  mar'])
            xlim([-1 1])
            export_fig(gcf,thisFigFileName,'-append','-pdf')
            delete(gcf)
            
            thisNeuronAlignedStas{fileNum,3} = extractEpochs(saccadeDetectOut.traces.velocity.v,thisFileStableSpikes{thisFileInfo.unitNumbers(2)},[4 4],calData.eyeFsCal,'displayFlag',1);
            set(get(gca,'title'),'string',[thisFileInfo.fileName 'Lin V  mar'])
            xlim([-1 1])
            export_fig(gcf,thisFigFileName,'-append','-pdf')
            delete(gcf)
             end
%             thisNeuronAlignedStas{fileNum,2} = extractEpochs(calData.eyeXcal,thisFileStableSpikes{thisFileInfo.unitNumbers(2)},[4 4],calData.eyeFsCal,'displayFlag',1);
%             set(get(gca,'title'),'string',[thisFileInfo.fileName 'Y pos'])
%              xlim([-1 1])
%             export_fig(gcf,figureFileName,'-append','-pdf')
%             delete(gcf)
            end
%             %list all trial starts falling within these
%             trialsToKeep = cellfun(@(x) x>thisStartTime & x<thisStopTime,{analysisResults(fileNum).newTrialStructure.trialStart},'uniformoutput',false);
%             stablePeriodTrials{fileNum} = analysisResults(fileNum).newTrialStructure([trialsToKeep{:}]);
%             
            %stablePeriodTrials{fileNum} =  burstAnalysisResults(fileNum).burstResults;
        end
        
%         thisNeuronAlignedRastersCombined = vertcat(thisNeuronAlignedRasters{:});
%         if ~isempty(thisNeuronAlignedRastersCombined)
%         %now plot the resulting combined raster
%         [hA] = createCombinedRaster(thisNeuronAlignedRastersCombined,[0.1 0.1]);
%         set(get(hA,'title'),'string',neuronList(neuronNum).neuronName,'interpreter','none')
%         thisFig = get(hA,'parent');
%         
%         
%         else
%             
%             thisFig = figure;
%             hA = axes('parent',thisFig,'xlim',[0 1],'ylim',[0 1]);
%             hT = text(0.5,0.5,'No CS!','parent',hA);
%             set(get(hA,'title'),'string',neuronList(neuronNum).neuronName,'interpreter','none')
%        
%         end
%         export_fig(allCsPausePdfName, '-pdf', '-append', thisFig);
%         %remove title before savinga s fig so don't have to deal with it
%         %when placing in summary fig
%         delete(get(hA,'title'))
%         %hgsave(thisFig,thisFigFileName)
%         %chagned to png as easier to add to summary fig
%         export_fig([resultsDir filesep neuronList(neuronNum).neuronName csPauseFigNameEnd], '-png', thisFig);
%         delete(thisFig)
end

end
%% CS PAUSE, WHOLE NEURON
%should get the complex spike pause for all files concatenated.

if wholeNeurCsPause

for neuronNum = recordingsToAnalyse;
   
    %load matfile into memory
     thisMatFileName = [processDir filesep neuronList(neuronNum).neuronName dataFileNameEnd];
     mLoadedData = matfile(thisMatFileName); %create mat file
   %build filename for analysis, not needed for CS pause
        thisAnalysisFileName = [resultsDir filesep neuronList(neuronNum).neuronName analysisFileNameEnd];
        
        thisFigFileName = [resultsDir filesep neuronList(neuronNum).neuronName csPauseFigNameEnd];
        
        numFiles = size(mLoadedData.editedFileList,1);
        load(thisAnalysisFileName);
      %  stablePeriodTrials = cell(1,numFiles);
        
        thisNeuronAlignedRasters = cell(numFiles,1);
        %loop over files and delete trials falling otuside marked times
        for fileNum = 1:numFiles
            
            %get file info
            thisFileInfo = mLoadedData.editedFileList(fileNum,1);
            
            %get start and end of stable section in this file
            if isinf(thisFileInfo.startTime)
                thisStartTime = -inf;
            else
                thisStartTime = thisFileInfo.startTime;
            end
            thisStopTime = thisFileInfo.stopTime;
            
            
            thisFileStableSpikes = cellfun(@(x) x(x>thisStartTime & x<thisStopTime),analysisResults(fileNum).alignedSpikeTimesCell,'uniformoutput',false);
            
            %now we need to get the cs pauase for this file
            if length(thisFileInfo.unitNumbers)>1
                %hF = figure;
                [hRast thisNeuronAlignedRasters{fileNum}] = createRaster(thisFileStableSpikes{thisFileInfo.unitNumbers(1)},...
                    thisFileStableSpikes{thisFileInfo.unitNumbers(2)},[0.1 0.1],'spikeMarkSize',3); %TODO also make this plot other complex spikes
                 title(thisFileInfo.fileName(1:end-4))
                %export_fig(figureFileName, '-pdf', '-append', gcf);
                delete(gcf)
                
                
            end
%             %list all trial starts falling within these
%             trialsToKeep = cellfun(@(x) x>thisStartTime & x<thisStopTime,{analysisResults(fileNum).newTrialStructure.trialStart},'uniformoutput',false);
%             stablePeriodTrials{fileNum} = analysisResults(fileNum).newTrialStructure([trialsToKeep{:}]);
%             
            %stablePeriodTrials{fileNum} =  burstAnalysisResults(fileNum).burstResults;
        end
        
        thisNeuronAlignedRastersCombined = vertcat(thisNeuronAlignedRasters{:});
        if ~isempty(thisNeuronAlignedRastersCombined)
        %now plot the resulting combined raster
        [hA] = createCombinedRaster(thisNeuronAlignedRastersCombined,[0.1 0.1]);
        set(get(hA,'title'),'string',neuronList(neuronNum).neuronName,'interpreter','none')
        thisFig = get(hA,'parent');
        
        
        else
            
            thisFig = figure;
            hA = axes('parent',thisFig,'xlim',[0 1],'ylim',[0 1]);
            hT = text(0.5,0.5,'No CS!','parent',hA);
            set(get(hA,'title'),'string',neuronList(neuronNum).neuronName,'interpreter','none')
       
        end
        export_fig(allCsPausePdfName, '-pdf', '-append', thisFig);
        %remove title before savinga s fig so don't have to deal with it
        %when placing in summary fig
        delete(get(hA,'title'))
        %hgsave(thisFig,thisFigFileName)
        %chagned to png as easier to add to summary fig
        export_fig([resultsDir filesep neuronList(neuronNum).neuronName csPauseFigNameEnd], '-png', thisFig);
        delete(thisFig)
end

end

%% Trial by trial correlations on population
if doPopCorMatrix
    a =1;
    
    
end
%% Trial BY TRAIL CORRELATION MATRIX a la MICHIEL
%for every neuron we loop over every trial and get the simple spike sdf
%and the eye traces. then normalise the trials eye trace to it's own
%baseline (/mean). TODO check with micheil this is what he did.
%now over all trials in this neuron get the R^2 correlation between the sfd
%value (bin x) and the eye value (bin y). Now plot the r2 on a colour
%square (imagesec)
%for now the correlation will be with the x position value of the eye in
%only correct pro trials
if doCorMatrix
corFigNameStart = 'corMatrix-';
corSigLevel = 0.01;
saccadeInd = 4000;
binSize = 0.005;
binSizeSamples = binSize*1000;
binningLims = [0.5 0.5]; %pre and post saccade time in s
binningLimsSamples = [saccadeInd-binningLims(1)*1000 saccadeInd+binningLims(2)*1000]; %in trial relative samples
numBins = sum(binningLims)/binSize;
%baseline [-1.25 -0.75] relative to saccade. saccade is at 4000 samples

for eyeVariableToUse = {'eyePositionX','eyePositionY','eyeVelocityX','eyeVelocityY','eyeVelocityV','eyeAngle'};

allNeuronsCorrResultsPvals = cell(max(recordingsToAnalyse),1);
allNeuronsCorrResults = cell(max(recordingsToAnalyse),1);
baselineInds = [2750 3250];
for neuronNum = recordingsToAnalyse
    thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
    unitNum = thisNeuronUnits(1);
    theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
    
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
    
    [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseCorProTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
    theseSdfs = vertcat(alignedTrials.Pro.sdf); %at both traces
    theseEyeTraces = vertcat(alignedTrials.Pro.(char(eyeVariableToUse)));
    %get baseline spike rates
    theseBaseRates = mean(theseSdfs(:,baselineInds(1):baselineInds(2)),2);
    %normalise to baseline trial by trial
    %theseNormSdfs = bsxfun(@ldivide,theseBaseRates,theseSdfs); %TODO swtich norm
    theseNormSdfs = bsxfun(@minus,theseSdfs,theseBaseRates); %TODO swtich norm
    numTrials = size(alignedTrials.Pro,2);
    thisSingleResult = {nan(numTrials,2)}; %first sdf secod eye
    theseResults = repmat(thisSingleResult,numBins,numBins); %cell to hold all results
    theseCorrResults = nan(numBins,numBins); %matrix to hold r values
    theseCorrResultsPvals = nan(numBins,numBins); %matrix to hold p values
    for xBinNum = 1:numBins
        %in the xbins we select part of the sdf to compare to
        %various time points of eye trace
        xBinInds = [binningLimsSamples(1)+binSizeSamples*(xBinNum-1)+1 binningLimsSamples(1)+binSizeSamples*(xBinNum)];
        %do gibbs outlier removal 
        x_date = cellfun(@(x) num2str(x),num2cell(1:numTrials),'uniformoutput',false);
        theseSdfValues = nanmean(theseNormSdfs(:,xBinInds(1):xBinInds(2)),2);
         %we need to ignore any trials witha  nan in the sdf
            badTrialsNanSdf = find(isnan(theseSdfValues));
        if ~isempty(badTrialsNanSdf)
            goodTrialsNanSdf = setdiff(1:numTrials,badTrialsNanSdf);
        else
             goodTrialsNanSdf = 1:numTrials;
        end
        [Gmax_loop Gmin_loop G2_loop lambda1_loop lambda2_loop outlier outlier_num] = grubbs(theseSdfValues, x_date');
       if ~isempty(outlier_num)
           badTrials = outlier_num(:,1);
       else
           badTrials = [];
       end
       goodTrialsOut = setdiff(1:numTrials,badTrials);
       goodTrials  = intersect(goodTrialsNanSdf,goodTrialsOut); 
        for yBinNum = 1:numBins
            %ybins select the part of the ey trace to use
            yBinInds = [binningLimsSamples(1)+binSizeSamples*(yBinNum-1)+1 binningLimsSamples(1)+binSizeSamples*(yBinNum)];
            theseEyeValues =  nanmean(theseEyeTraces(:,yBinInds(1):yBinInds(2)),2);
            
            %we need to ignore any trials witha  nan in the eye trace
            goodTrialsNan = find(~isnan(theseEyeValues));
            
            trialsToUse = intersect(goodTrialsNan,goodTrials);
            %now save actual trial by trials values as well as the R^2
            theseResults{xBinNum,yBinNum} = [theseSdfValues theseEyeValues];
            [theseCorrResults(xBinNum,yBinNum) theseCorrResultsPvals(xBinNum,yBinNum)] = corr(theseSdfValues(trialsToUse) ,theseEyeValues(trialsToUse));
        end
    end
    allNeuronsCorrResults{neuronNum} = theseCorrResults; %save results matrix
   allNeuronsCorrResultsPvals{neuronNum} = theseCorrResultsPvals;
    hF = figure;
   colormap(jet)
    hE = axes('parent',hF,'units','normalized','position',[0.05 0.3 0.2 0.6]);
    hS = axes('parent',hF,'units','normalized','position',[0.3 0.05 0.48 0.2]);
     hA = axes('parent',hF,'units','normalized','position',[0.3 0.3 0.6 0.6]);
    imagesc(theseCorrResults,[-0.5 0.5])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',[neuronList(neuronNum).neuronName ' ' char(eyeVariableToUse)])
    line([0 400],[200 200],'color','k','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','k','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    
    timeVec = linspace(binningLimsSamples(1),binningLimsSamples(2),diff(binningLimsSamples));
    line(theseEyeTraces(:,binningLimsSamples(1):binningLimsSamples(2)-1),timeVec ,'parent',hE)
    line(timeVec ,theseNormSdfs(:,binningLimsSamples(1):binningLimsSamples(2)-1),'parent',hS)
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    %also drwa pvalues
     hF = figure;
    hA = axes('parent',hF);
    imagesc(theseCorrResultsPvals,[0 1])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',[neuronList(neuronNum).neuronName ' ' char(eyeVariableToUse) ' p-vals'])
    line([0 400],[200 200],'color','k','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','k','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
    %also draw only significant p's
     hF = figure;
    hA = axes('parent',hF);
    sigPvals = theseCorrResultsPvals< corSigLevel;
    colormap([ 1 1 1; 0 0 0])
    imagesc(sigPvals,[0 1])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',[neuronList(neuronNum).neuronName ' ' char(eyeVariableToUse) ' if sig'])
    line([0 400],[200 200],'color','r','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','r','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
end

%TOOD try abs value of correlation? or pvalue?
%now draw average of all, also split average into lat and verm
fullCorrRes = cat(3,allNeuronsCorrResults{recordingsToAnalyse});
theseCorrResults = nanmean(abs(fullCorrRes),3);

 hF = figure;
  colormap(jet)
    hA = axes('parent',hF);
    imagesc(theseCorrResults,[-0.5 0.5])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',['All Neurons ' char(eyeVariableToUse)])
    line([0 400],[200 200],'color','k','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','k','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
    
    fullCorrRes = cat(3,allNeuronsCorrResults{vermNeurons});
theseCorrResults = nanmean(abs(fullCorrRes),3);

 hF = figure;
    hA = axes('parent',hF);
    imagesc(theseCorrResults,[-0.5 0.5])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',['Vermis Neurons ' char(eyeVariableToUse)])
    line([0 400],[200 200],'color','k','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','k','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
     fullCorrRes = cat(3,allNeuronsCorrResults{latNeurons});
theseCorrResults = nanmean(abs(fullCorrRes),3);

 hF = figure;
    hA = axes('parent',hF);
    imagesc(theseCorrResults,[-0.5 0.5])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',['Lateral Neurons ' char(eyeVariableToUse)])
    line([0 400],[200 200],'color','k','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','k','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
    %now draw same three graphs for p-values
    fullCorrRes = cat(3,allNeuronsCorrResultsPvals{recordingsToAnalyse});
theseCorrResults = nanmean(fullCorrRes,3);

 hF = figure;
    hA = axes('parent',hF);
    imagesc(theseCorrResults,[0 1])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',['All Neurons ' char(eyeVariableToUse) ' p-vals'])
    line([0 400],[200 200],'color','k','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','k','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
    
    fullCorrRes = cat(3,allNeuronsCorrResultsPvals{vermNeurons});
theseCorrResults = nanmean(fullCorrRes,3);

 hF = figure;
    hA = axes('parent',hF);
    imagesc(theseCorrResults,[0 1])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',['Vermis Neurons ' char(eyeVariableToUse) ' p-vals'])
    line([0 400],[200 200],'color','k','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','k','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
    
     fullCorrRes = cat(3,allNeuronsCorrResultsPvals{latNeurons});
theseCorrResults = nanmean(fullCorrRes,3);

 hF = figure;
    hA = axes('parent',hF);
    imagesc(theseCorrResults,[0 1])
    set(hA,'YDir','normal')
    colorbar
    set(get(hA,'title'),'string',['Lateral Neurons ' char(eyeVariableToUse) ' p-vals'])
    line([0 400],[200 200],'color','k','linestyle','-.','parent',hA)
    line([200 200],[0 400],'color','k','linestyle','-.','parent',hA)
    set(get(hA,'ylabel'),'string',eyeVariableToUse)
    set(get(hA,'xlabel'),'string','SS norm sdf')
    export_fig([resultsDir filesep corFigNameStart char(eyeVariableToUse) '.pdf'], '-pdf','-append', hF);
    delete(hF)
  
    
    
    allVarCorRes.(char(eyeVariableToUse)) = allNeuronsCorrResults;
allVarCorResPvals.(char(eyeVariableToUse)) = allNeuronsCorrResultsPvals;
    
end
save([resultsDir filesep 'corResTest2.mat'],'allVarCorRes','allVarCorResPvals')
end
%% SPIKE PSDS
doSpikeSpectra = 0;
if doSpikeSpectra
indTrialPsdFlag =1;
psdFs = 2000;
% Set power for segment length as 10.
% Thus T=2^10 = 1024.
% 1024 points @ 1000/sec sampling gives segment length of 1.024 sec
% Frequency resolution is inverse of this 0.977 Hz.
%nextpow2(2000) = 11 for us
seg_pwr=9;
trialLength = 8; %length of full trial in seconds
opt_str='';    % Clear options %TODO read and understand
%spike frequency analysis using neurospec
for neuronNum = recordingsToAnalyse(1:3);
    theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
    
    numTrials = size(theseStableTrials,2);
    stInit = cell(numTrials,1);
    neuronPsds = struct('freqAxis1',stInit,'wholeTrain',stInit,'freqAxis2',stInit,'period1',stInit,'period2',stInit)
    %now loop over trials and get the frequency content of each 
    for trialNum = 1:numTrials
        
        %we have 4 seconds before and after the start so first lets use all
        %of it
        %TODO need to convert spikes to a vector of integer times in
        %samples
        if ~isempty(theseStableTrials(trialNum).alignedSpikes)
        thisTrialSpikes = theseStableTrials(trialNum).alignedSpikes(1);
        
        spikesForTesting = round(psdFs*thisTrialSpikes{:}); %+1 to ensuse no 0 spike
        %however alignedSpikes start with negative numbers so align all to
        %first spike %TODO need to add a delay to stop loss of first
        %spike??
        spikesForTesting = spikesForTesting+abs(min(spikesForTesting))+1;
        
        sp1 = spikesForTesting;
        sp2 = spikesForTesting;
        
        [f4,t4,cl4]=sp2_m1(0,sp1,sp2,trialLength,psdFs,seg_pwr,opt_str);
        %cl4.what='Independent Spike Trains';
%         figure
        freq = 200;     % Parameters for plotting
%         lag_tot=100;
%         lag_neg=50;
%         ch_max=0.5;
%         psp2(f4,t4,cl4,freq,lag_tot,lag_neg,ch_max)
        
       
        neuronPsds(trialNum).freqAxis1 = f4(:,1);
        neuronPsds(trialNum).wholeTrain = f4(:,2);
        
        trialLengthPoints = 550;
        %this time split into 2 sections
        sp3 = spikesForTesting(spikesForTesting<8000&spikesForTesting>(8000-trialLengthPoints));
        sp4 = spikesForTesting(spikesForTesting>=8000&spikesForTesting<=8000+trialLengthPoints);
        trialLength = trialLengthPoints/2000;
        
        sp3 = sp3-abs(sp3(1))+1;
        sp4 = sp4-abs(sp4(1))+1;
        [f5,t5,cl5]=sp2_m1(0,sp3,sp4,trialLength,psdFs,seg_pwr,opt_str);
        neuronPsds(trialNum).freqAxis2 = f5(:,1);
        neuronPsds(trialNum).period1 = f5(:,2);
        neuronPsds(trialNum).period2 = f5(:,3);
         
        if indTrialPsdFlag
        hF = figure;
        subplot(1,3,1)
        psp_fa1(f4,cl4,freq); %power spectrum plot subfunction of above function
        subplot(1,3,2)
         psp_fa1(f5,cl5,freq); %power spectrum plot subfunction of above function
         subplot(1,3,3)
         psp_fb1(f5,cl5,freq); %power spectrum plot subfunction of above function
         
          %set(get(gca,'title'),'string',['Trial ' num2str(trialNum)]);
          suptitle([neuronList(neuronNum).neuronName '-Trial ' num2str(trialNum)])
          %TODO export fig
          delete(hF)
        end
        end 
    end
    
    %get frequency axis of first trial as all should eb equal %TODO need
    %error check here and for empty first trial erorr
    f1 = neuronPsds(1).freqAxis1;
    allPsds = [neuronPsds.wholeTrain];
    muPsd = mean(allPsds,2);
    hF2 = figure;
    subplot(1,3,1)
   plot(f1,muPsd);
   xlim([0 200])
    f2 = neuronPsds(1).freqAxis1;
    allPsds = [neuronPsds.period1];
    muPsd = mean(allPsds,2);
    subplot(1,3,2)
     plot(f2,muPsd);
     xlim([0 200])
     allPsds = [neuronPsds.period2];
    muPsd = mean(allPsds,2);
    subplot(1,3,3)
     plot(f2,muPsd);
     xlim([0 200])
     
     suptitle([neuronList(neuronNum).neuronName '-AllTrials'])
    
end
end
%% BURST PLOTS
% if doBursts
% if burstSumPlotFlag
%     
%     for neuronNum = recordingsToAnalyse;
%         %get all of teh correct trials with a valid saccade
%         theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
%         theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
%         allFrs = {theseStableTrials(theseTrialNums).overlapFr}
%         %now get extremes to set limits of spider of plots
%         %allFrs = {allNeuronResults(neuronNum).allStableBurstTrials(corSacTrialNums).overlapFr};
%         spiderLims = [min([allFrs{:}]) max([allFrs{:}])]*1.5;
%         %do for pro and anti sperately
%         hF2 = figure;
%         muBurstPeakRateByDir = nan(9,2); %to keep the saccade period firing rae for each direction
%         semBurstPeakRateByDir = nan(9,2); %for SEMs
%         burstPeakPercent = nan(9,2); %percentage of trials with a dtetected burst in the saccade period
%         ratesToPlot = cell(2,1); %to hold rates at each direction for both pro and anti for tuning curve
%         errToPlot = cell(2,1); %to hold SEMs for conf limits on polar plot
%         for dirLabelNum = 1:2
%             
%             
%             for plNum = find(~isnan(trialConditionCodes(:,dirLabelNum)))'
%                 %find all trials fort his direction witha  valid saccade
%                 thisConditionCode = trialConditionCodes(plNum,dirLabelNum);
%                 tempLogCell = cellfun(@(x) x==thisConditionCode,{theseStableTrials.conditionCode},'uniformoutput',false);
%                 theseDirTrialNums = find([tempLogCell{:}]);
%                 corDirTrialNums = intersect(theseTrialNums,theseDirTrialNums);
%                 
%                 
%                 [hAx1 hAx2 hAx22] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
%                     corDirTrialNums,smallTimeWindow,timeWindow,'burstMarkingFlag',1)
%                 hFtemp = get(hAx1,'parent');
%                 set(get(hAx1,'title'),'string',['Cond ' num2str(thisConditionCode)])
%                 export_fig(hFtemp,'-pdf','-append',burstSumFigName)
%                 delete(hFtemp)
%                 
%                 
%                 %get mean of peak firing rate and numbre of non-nan trials to
%                 %calcualte SEMs
%                 burstPeakFiringRates = {theseStableTrials(corDirTrialNums).overlapFr};
%                 numBurstTrials = sum(cellfun(@(x) any(x),burstPeakFiringRates));
%                 numTrials = size(burstPeakFiringRates,2);
%                 burstPeakPercent(plNum,dirLabelNum) = 100*(numBurstTrials)/(numTrials);
%                 theseBurstFrs = [burstPeakFiringRates{:}];
%                 muBurstPeakRateByDir(plNum,dirLabelNum) = nanmean(theseBurstFrs(:)); %TODO check how a nan ets through
%                 semBurstPeakRateByDir(plNum,dirLabelNum) = nanstd(theseBurstFrs(:))/sqrt(numBurstTrials);
%                 %now on the outisde mark the percentage of trials witha  valid
%                 %burst
%                 set(0, 'currentfigure', hF2);
%                 hAs = subaxis(3,3,plNum)
%                 text(0.5,0.2+dirLabelNum*0.2,num2str(burstPeakPercent(plNum,dirLabelNum)),'parent',hAs,'color',ccProAnti(dirLabelNum,:))
%                 set(hAs,'visible','off')
%             end
%             
%             
%             theseRatesToPlot = muBurstPeakRateByDir([1 2 3 4 6 7 8 9],dirLabelNum);
%             theseErrsToPlot = semBurstPeakRateByDir([1 2 3 4 6 7 8 9],dirLabelNum);
%             
%             
%             %spider starts at 3 o'clock and goes anti-clockwise
%             %plNum starts at 10 oclock and goes clockwise
%             %so flipud
%             ratesToPlot{dirLabelNum} = flipud(circshift(theseRatesToPlot,3));
%             errToPlot{dirLabelNum} = flipud(circshift(theseErrsToPlot,3));
%             %
%             %                      [f, ca, o] = createSpider(ratesToPlot ,'',[0 spiderLims(2)],repmat({''},8,1),'');
%             %             hT = copyobj(ca,hF2);
%             %             set(hT,'units','normalized','position',[get(hA,'position') + [0 -0.1 0 0.1]])
%         end
%         %Create tuning plot with coloured areas representing confidnce limits
%         hA = subaxis(3,3,5,'ML',0,'MR',0,'MT',0,'MB',0)
%         anglesToTest = linspace(0,2*pi,9);
%         anglesToTest = anglesToTest(1:8); %8 evenly spaced going counterclockwise
%         
%         
%         [hA] = createConfLimSpider(ratesToPlot,errToPlot,anglesToTest,hA);
%         maxVal = max(max([ratesToPlot{:}]));
%         set(hA,'color','none','xlim',[-maxVal maxVal],'ylim',[-maxVal maxVal])
%         
%         
%         export_fig(hF2,'-pdf','-append',burstSumFigName)
%         
%         %now crate tuning plot with signinficance test using circular
%         %statistics
%         hFsig = figure;
%         hAsig = axes('parent',hFsig);
%         anglesToTest = linspace(0,2*pi,9);
%         anglesToTest = anglesToTest(1:8)
%         [hAsig] = createTuningPlot(anglesToTest',ratesToPlot,hAsig)
%         
%         export_fig(hFsig,'-pdf','-append',burstSumFigName)
%         
%         delete(hFsig)
%         delete(hF2)
%         
%         %    circStats.circMean = circ_mean(ori,ratesToPlot',2);
%         %         %   circStats.circRtest = circ_rtest(ori,ratesToPlot',dori);
%         %         %
%         %
%     end
%     
% end
% %TODO fix so that it saves a pdf per file or teh seperate conditions not
% %for entire set of recordings, should also save to summary the tuning plots
% %with titles.
% %TODO circular stats not working. Fix
% end
%% REWARD AND ITI
% %TODO in saccade trimming get the next saccade after the important one as
% this interferes with the reward period
% for neuronNum = recordingsToAnalyse;
%     %ALL SACCADES BIG WINDOW
%     [hAx1 hAx2 hA2] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
%         wholeNeuronResults(neuronNum).selectedTrials.corSacTrials,[0.3 0.3],timeWindow,'alignTo','relativeRewardTime','markEvent','saccadeTime','sortBy','saccadeTime');%,...
%        %  'sortBy','goCueTime','alignTo','bit6time','markEvent','goCueTime')
%     title([neuronList(neuronNum).neuronName ' all saccades, align to rewardTime'])
%     export_fig(gcf,'-pdf','-append','rewardRasters2.pdf')
%     
%     
%     
%     
% end
%% CORRELATIONS
%do correlations between main burst features: peak, duration, start time,
%                                                           end time
%                   and saccade features: peakVel, duration, velocity

if doCor
    
allCorrResults = struct('kuniRes',[],'kuniResultsStats',[]);

for neuronNum = recordingsToAnalyse;
    %build figure filename
    thisNeuronSumFigName = [resultsDir filesep neuronList(neuronNum).neuronName indNeuronSumFigEnd];
    
    %get the trials and the numbers of stable ones
    theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
    
    theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
    theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
    ccLog = zeros(size(theseStableTrials,2),1);
    ccLog(theseAntiTrialNums) = 2;
    ccLog(theseProTrialNums) = 1;
    %need to make sure that only one burst is included for correlation
    %purposes
    [trialBurstStruct] = selOneBurst(theseStableTrials);
    theseSelBursts = [trialBurstStruct.selBurst];
    
    
    burstStatsToTest = {'relStartTime','duration','relEndTime','overlapFr'};
    statsToTest = {'saccadeAmplitude','saccadeDuration','saccadePeakVel','reactionTime'}
    hFc = figure('renderer','opengl');
    
    numBurstStats = size(burstStatsToTest,2)
    numStats = size(statsToTest,2)
    for burstStatsToTestNum = 1:numBurstStats;
        thisBurstStat = burstStatsToTest{burstStatsToTestNum};
        for statsToTestNum = 1:numStats
            thisStat = statsToTest{statsToTestNum};
            testData1 = {theseStableTrials.(thisStat)};
            %testData2 = {theseStableTrials.overlapFr};
            testData2 = {theseSelBursts.(thisBurstStat)};
            
            %now check for which trials have valid statistics (independent first,
            %then intersection)
            
            tempLogCell = cellfun(@(x) any(x),testData1,'uniformoutput',false);
            
            theseValidData1 = find([tempLogCell{:}]);
            tempLogCell = cellfun(@(x) any(x),testData2,'uniformoutput',false);
            
            theseValidData2 = find([tempLogCell{:}]);
            
            theseValidDataBoth = intersect(theseValidData1,theseValidData2);
            %find all correct trials (btoth pro and anti)
            theseCorTrials = find(ccLog~=0);
            
            theseValidTrials = intersect(theseCorTrials,theseValidDataBoth);
            
            
            hAxc = subaxis(numBurstStats,numStats,sub2ind([numStats numBurstStats],statsToTestNum,burstStatsToTestNum),'S',0.1);
            
            %now create a scatter graph using the vecot rof groups from ealier
            %to control colours
            scatter([testData1{theseValidTrials}],[testData2{theseValidTrials}],...
                10,ccProAnti(ccLog(theseValidTrials),:),'parent',hAxc)
            
            validProTrials = intersect(theseProTrialNums,theseValidTrials);
            validAntiTrials = intersect(theseAntiTrialNums,theseValidTrials);
            
            if ~isempty(validProTrials) && ~isempty(validAntiTrials)
            proVals1 = [testData1{validProTrials}];
            proVals2 = [testData2{validProTrials}];
            antiVals1 = [testData1{validAntiTrials}];
            antiVals2 = [testData2{validAntiTrials}];
            
            antiCoeffs = polyfit(antiVals1,antiVals2,1)
            fittedAntiX = linspace(min(antiVals1), max(antiVals1), 200);
            fittedAntiY = polyval(antiCoeffs, fittedAntiX);
            
            line(fittedAntiX,fittedAntiY,'color',ccProAnti(2,:),'parent',hAxc,...
            'linewidth',2);
            
            proCoeffs = polyfit(proVals1,proVals2,1)
            fittedProX = linspace(min(proVals1), max(proVals1), 200);
            fittedProY = polyval(proCoeffs, fittedProX);
            
            line(fittedProX,fittedProY,'color',ccProAnti(1,:),'parent',hAxc,...
                'linewidth',2)
            %lsline %TODO replace with smarter fitting
            
            %only label appropraite ones
            if statsToTestNum==1
                ylabel(thisBurstStat)
            end
            if burstStatsToTestNum==4
                xlabel(thisStat)
                
            end
            %also label with R value of correlation
            %        RHO = corr([[testData1{theseValidDataBoth}]' [testData2{theseValidDataBoth}]']);
            %        title(['R = ' num2str(RHO(1,2))])
            %get R for pro and anti trials seperately
            rPro = corr([proVals1' proVals2']);
            rAnti = corr([antiVals1' antiVals2']);
            corAnti =roundTo(rAnti(1,2),0.01)
            corPro = roundTo(rPro(1,2),0.01)
            curYlim = get(hAxc,'ylim')
            curXlim = get(hAxc,'xlim')
            text(curXlim(1),curYlim(2)*1.1,['R=' num2str(corPro)],...
                'color',ccProAnti(1,:),'fontsize',8)
            text(curXlim(1)+diff(curXlim)/2,curYlim(2)*1.1,...
                ['R=' num2str(corAnti)],'color',ccProAnti(2,:),'fontsize',8)
            %        title(['R = ' num2str(RHO(1,2))])
            end
        end
        
    end
    suptitle(neuronList(neuronNum).neuronName)
    %TODO export as png as well as resolution is poor
    export_fig(hFc,'-pdf','-append','-opengl',[resultsDir filesep corrPdfName])
    delete(hFc)
end
end
%% LFP
% 
% for neuronNum = recordingsToAnalyse
%     %get the trials and the numbers of stable ones
%     theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
%     theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
%    [res] = specFromTrial(theseStableTrials(theseTrialNums),timeWindow) 
%    set(get(gca,'title'),'string',[neuronList(neuronNum).neuronName 'all'])
%     export_fig(gcf,'-pdf','-append','-opengl','ericSpecgramTesting.pdf')
%     
%      theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
%     theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
%    [res] = specFromTrial(theseStableTrials(theseTrialNums),timeWindow) 
%   set(get(gca,'title'),'string',[neuronList(neuronNum).neuronName 'anti'])
%     export_fig(gcf,'-pdf','-append','-opengl','ericSpecgramTesting.pdf')
%     
%      theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
%     theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
%    [res] = specFromTrial(theseStableTrials(theseTrialNums),timeWindow) 
%   set(get(gca,'title'),'string',[neuronList(neuronNum).neuronName 'pro'])
%     export_fig(gcf,'-pdf','-append','-opengl','ericSpecgramTesting.pdf')
% end
%%
%% CV2 POPULATION RESPONSE


 

sortByString = 'saccadeAmplitude';
frDispRange = [0 200;0 5;0 200;0 200;0 5];
intTimeVec = -4:0.001:4; %ms spacing

unitNum = 1;

if doPopCV2
    
    
    
    
    %if exist(popResFileName,'file')~=2
    stInit = cell(max(recordingsToAnalyse),2);
    popResCV2 = struct('alignedTrials',stInit,'antiMu',stInit,'proMu',stInit);
    for neuronNum = recordingsToAnalyse;
        thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
        %unitCount =1;
        % for unitNum  = thisNeuronUnits;
        %get the trials and the numbers of stable ones
        theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
        theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        
        theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
        theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        
        
        [alignedTrials.Anti antiMu] = alignAndExtractRates(theseStableTrials,theseAntiTrialNums,'saccadeTime','unitNum',thisNeuronUnits(unitNum));
        [alignedTrials.Pro proMu] = alignAndExtractRates(theseStableTrials,theseProTrialNums,'saccadeTime','unitNum',thisNeuronUnits(unitNum));
        
        %after getting the aligned trials just need to extract the
        %cv2
        
        
        
        hFs = figure; %to hold both pro and anti heatplots
        [cv2mat hFC hAxH alignedTrials.Pro] = extractCV2fromAlignedTrials(alignedTrials.Pro,intTimeVec,thisNeuronUnits(unitNum),sortByString);
        delete(hFC)
        popResCV2(neuronNum,1).proMu = nanmean(cv2mat);
        
        hAxS1 = subplot(1,2,1,'parent',hFs);
        hAxS12 = copyobj(get(hAxH,'children'),hAxS1);
        delete(get(hAxH,'parent'))
        [cv2mat hFC hAxH alignedTrials.Anti] = extractCV2fromAlignedTrials(alignedTrials.Anti,intTimeVec,thisNeuronUnits(unitNum),sortByString);
        delete(hFC)
        
        popResCV2(neuronNum,1).antiMu = nanmean(cv2mat);
        
        set(hAxS1,'ylim',[0 size(cv2mat,1)])
        hAxS2 = subplot(1,2,2,'parent',hFs);
        hAxS22 = copyobj(get(hAxH,'children'),hAxS2);
        delete(get(hAxH,'parent'));
        set([hAxS1 hAxS2],'xlim',[axisTicks(1) axisTicks(3)],'clim',[0 2],'xtick',axisTicks,'xticklabel',axTickLabs)
        set(get(hAxS1,'title'),'string','Pro Trials CV2')
        set(get(hAxS2,'title'),'string','Anti Trials CV2')
        suptitle(recList{1,neuronNum})
        set(hAxS2,'ylim',[0 size(cv2mat,1)])
        try
            export_fig(hFs,'-pdf',expRend,'-append',cvPopResPdfName)
        catch
            export_fig(hFs,'-pdf',expRend,cvPopResPdfName)
        end
        delete(hFs)
        
        
        close all
    end
    
    %save([resultsDir filesep popResFileName],'popRes','-v7.3')
    
    %else
    %   load([resultsDir filesep popResFileName]);
    %end
    %
    %
    tV = intTimeVec; %get the itnerpolated time vector for the cv2 traces
    hFd = figure; %for showing the difference between pro and anti means
    hAd(1) = subplot(1,2,1,'parent',hFd);
    hAd(2) = subplot(1,2,2,'parent',hFd);
    hF = figure;
    hAxs =nan(1,4);
    
    %first plot the mean response to all pro trials for vermis, then
    %lateral, on a seperate figure also plot the difference between the two
    hAxs(1) = subplot(2,2,1,'parent',hF);
    hAxs(2) = subplot(2,2,2,'parent',hF);
    
    proMus = vertcat(popResCV2(vermNeurons).proMu);
    antiMus =  vertcat(popResCV2(vermNeurons).antiMu);
    pMinusA = proMus-antiMus;
    
    line(tV',antiMus,'parent',hAxs(2))
    
    line(tV',proMus,'parent',hAxs(1))
    line(tV',pMinusA,'parent',hAd(1))
    
    %now also plot mean of this as thick black line
    line(tV',nanmean(antiMus),'parent',hAxs(2),'color','k','linewidth',2)
    
    line(tV',nanmean(proMus),'parent',hAxs(1),'color','k','linewidth',2)
    line(tV',nanmean(pMinusA),'parent',hAd(1),'color','k','linewidth',2)
    set(get(hAxs(2),'title'),'string','Anti')
    set(get(hAxs(1),'title'),'string','Pro')
    set(get(hAxs(2),'ylabel'),'string','Vermis')
    set(get(hAd(1),'title'),'string','Vermis')
    
    
    hAxs(3) = subplot(2,2,3,'parent',hF);
    hAxs(4) = subplot(2,2,4,'parent',hF);
    
    proMus = vertcat(popResCV2(latNeurons).proMu);
    antiMus =  vertcat(popResCV2(latNeurons).antiMu);
    pMinusA = proMus-antiMus;
   
    line(tV',antiMus,'parent',hAxs(4))
    
    line(tV',proMus,'parent',hAxs(3))
    line(tV',pMinusA,'parent',hAd(2))
    
    %now also plot mean of this as thick black line
    line(tV',nanmean(antiMus),'parent',hAxs(4),'color','k','linewidth',2)
    
    line(tV',nanmean(proMus),'parent',hAxs(3),'color','k','linewidth',2)
    line(tV',nanmean(pMinusA),'parent',hAd(2),'color','k','linewidth',2)
    
    set(get(hAxs(4),'title'),'string','Anti')
    set(get(hAxs(3),'title'),'string','Pro')
    set(get(hAxs(4),'ylabel'),'string','Lateral')
    set(get(hAd(2),'title'),'string','Lateral')
    temp = cell2mat(get(hAxs,'ylim'));
    set(hAxs,'ylim',[0 1.5],'xlim',[-1 1])
    
    
    suptitle(['CV2 Unit ' num2str(unitNum)])
    export_fig(hF,'-pdf','-opengl','-append',cvPopResPdfName)
    delete(hF)
    set(hAd,'xlim',[-1 1],'ylim',[-0.5 0.5])
    suptitle(['P-A CV2 Unit ' num2str(unitNum)])
    export_fig(hFd,'-pdf','-opengl','-append',cvPopResPdfName)
    delete(hFd)
    
    
    
    
    %we will now do exactly the same thing but for normlaised data
    %isntead
    
    %  normData = normToBaseInds(data,baseInd) %us this function to
    %  normalise
    
    
    hFd = figure; %for showing the difference between pro and anti means
    hAd(1) = subplot(2,1,1,'parent',hFd);
    hAd(2) = subplot(2,1,2,'parent',hFd);
    hF = figure;
    hAxs =nan(1,4);
    
    %first plot the mean response to all pro trials for vermis, then
    %lateral, on a seperate figure also plot the difference between the two
    hAxs(1) = subplot(2,2,1,'parent',hF);
    hAxs(2) = subplot(2,2,2,'parent',hF);
    
    proMus = vertcat(popResCV2(vermNeurons).proMu);
    antiMus =  vertcat(popResCV2(vermNeurons).antiMu);
    
    
    normProMus = normToBaseInds(proMus,baseInd);
    normAntiNums = normToBaseInds(antiMus,baseInd);
    pMinusA = normProMus-normAntiNums;
    line(tV',normAntiNums,'parent',hAxs(2))
    
    line(tV',normProMus,'parent',hAxs(1))
    line(tV',pMinusA,'parent',hAd(1))
    vermPvA = pMinusA;
    %now also plot mean of this as thick black line
    line(tV',nanmean(normAntiNums),'parent',hAxs(2),'color','k','linewidth',2)
    
    line(tV',nanmean(normProMus),'parent',hAxs(1),'color','k','linewidth',2)
    line(tV',nanmean(pMinusA),'parent',hAd(1),'color','k','linewidth',2)
    set(get(hAxs(2),'title'),'string','Anti')
    set(get(hAxs(1),'title'),'string','Pro')
    set(get(hAxs(2),'ylabel'),'string','Vermis')
    set(get(hAd(1),'title'),'string','Vermis')
    
    
    hAxs(3) = subplot(2,2,3,'parent',hF);
    hAxs(4) = subplot(2,2,4,'parent',hF);
    
    proMus = vertcat(popResCV2(latNeurons).proMu);
    antiMus =  vertcat(popResCV2(latNeurons).antiMu);
    normProMus = normToBaseInds(proMus,baseInd);
    normAntiNums = normToBaseInds(antiMus,baseInd);
    pMinusA = normProMus-normAntiNums;
     latPvA = pMinusA;
    
    line(tV',normAntiNums,'parent',hAxs(4))
    
    line(tV',normProMus,'parent',hAxs(3))
    line(tV',pMinusA,'parent',hAd(2))
    
    %now also plot mean of this as thick black line
    line(tV',nanmean(normAntiNums),'parent',hAxs(4),'color','k','linewidth',2)
    
    line(tV',nanmean(normProMus),'parent',hAxs(3),'color','k','linewidth',2)
    line(tV',nanmean(pMinusA),'parent',hAd(2),'color','k','linewidth',2)
    
    set(get(hAxs(4),'title'),'string','Anti')
    set(get(hAxs(3),'title'),'string','Pro')
    set(get(hAxs(4),'ylabel'),'string','Lateral')
    set(get(hAd(2),'title'),'string','Lateral')
    temp = cell2mat(get(hAxs,'ylim'));
    set(hAxs,'ylim',[-30 30],'xlim',[-1 1])
    set(hAxs,'ytick',[-30 0 30])
    
    suptitle(['Norm CV2 Unit ' num2str(unitNum)])
    export_fig(hF,'-pdf','-opengl','-append',cvPopResPdfName)
    delete(hF)
    set(hAd,'xlim',[-1 1],'ylim',[-15 15],'ytick',[-15 0 15])
    suptitle(['Norm P-A CV2 Unit ' num2str(unitNum)])
    export_fig(hFd,'-pdf','-opengl','-append',cvPopResPdfName)
    delete(hFd)
    
    
    %TIMINGS OF SIGNIFICANT POINTS
    %now neuron by neuron we should find the first time that each is abs>4,
    %also find the maximum
    sigStatsVermPvA = nan(size(vermPvA,1),2);
    sigStatsLatPvA = nan(size(latPvA,1),2);
   numStds = 5; %number of stds away from baseline is counted as significant
    for neuronNum = 1:size(vermPvA,1)
        
         firstSigPoint = find(abs(vermPvA(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsVermPvA(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(vermPvA(neuronNum,3500:5500)));
       
         sigStatsVermPvA(neuronNum,2)=maxSigPoint-500;
        end
    end
    %now for lateral
    for neuronNum = 1:size(latPvA,1)
        
         firstSigPoint = find(abs(latPvA(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsLatPvA(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(latPvA(neuronNum,3500:5500)));
       
         sigStatsLatPvA(neuronNum,2)=maxSigPoint-500;
        end
    end
    
    ccBar = repmat([0 0 1;1 0 0],2,1); %alternating r and b for bars
    vermMu = nanmean(sigStatsVermPvA);
    vermSem = nanstd(sigStatsVermPvA)./sqrt(size(vermPvA,1));
    latMu = nanmean(sigStatsLatPvA);
    latSem = nanstd(sigStatsLatPvA)./sqrt(size(latPvA,1));
    xIn = [vermMu' latMu']./1000; %convert to s
    erIn = [vermSem' latSem']./1000;
    hF = figure;
    hAx = axes('parent',hF);
   
   [hPatches] = createGroupedBarGraph(xIn,erIn,hAx,ccBar,0.3);
   suptitle('Time of CV2 differences PvA')
   set(hAx,'xtick',[1 2],'xticklabel',{'First Time','Maximum Time'})
   set(get(hAx,'ylabel'),'string','Time Relative to Saccade (s)')
     export_fig(hF,'-pdf','-opengl','-append',cvPopResPdfName)
    delete(hF)
  
    
    %now dot he same but instead of looking at pro vs anti differences
    %let's just get the times of modualtion from normalised baseline
    %this is for pro trials
    proMus = vertcat(popResCV2(latNeurons).proMu);
    normLatMus = normToBaseInds(proMus,baseInd);
    proMus = vertcat(popResCV2(vermNeurons).proMu);
    normVermMus = normToBaseInds(proMus,baseInd);
    
    
    sigStatsVerm = nan(size(normVermMus,1),2);
    sigStatsLat = nan(size(normLatMus,1),2);
   %numStds = 7; %number of stds away from baseline is counted as significant
    for neuronNum = 1:size(normVermMus,1)
        
         firstSigPoint = find(abs(normVermMus(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsVerm(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(normVermMus(neuronNum,3500:5500)));
       
         sigStatsVerm(neuronNum,2)=maxSigPoint-500;
        end
    end
    %now for lateral
    for neuronNum = 1:size(normLatMus,1)
        
         firstSigPoint = find(abs(normLatMus(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsLat(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(normLatMus(neuronNum,3500:5500)));
       
         sigStatsLat(neuronNum,2)=maxSigPoint-500;
        end
    end
    
    ccBar = repmat([0 0 1;1 0 0],2,1); %alternating r and b for bars
    vermMu = nanmean(sigStatsVerm);
    vermSem = nanstd(sigStatsVerm)./sqrt(size(vermPvA,1));
    latMu = nanmean(sigStatsLat);
    latSem = nanstd(sigStatsLat)./sqrt(size(latPvA,1));
    xIn = [vermMu' latMu']./1000; %convert to s
    erIn = [vermSem' latSem']./1000;
    hF = figure;
    hAx = axes('parent',hF);
   
   [hPatches] = createGroupedBarGraph(xIn,erIn,hAx,ccBar,0.3);
   suptitle('Time of CV2 differences From baseline, Pro Trials')
   set(hAx,'xtick',[1 2],'xticklabel',{'First Time','Maximum Time'})
   set(get(hAx,'ylabel'),'string','Time Relative to Saccade (s)')
     export_fig(hF,'-pdf','-opengl','-append',cvPopResPdfName)
    delete(hF)
     cv2PopTimeStats.sigStatsVermPro = sigStatsVerm;
    cv2PopTimeStats.sigStatsLatPro = sigStatsLat;
    
     %now dot he same but instead of looking at pro vs anti differences
    %let's just get the times of modualtion from normalised baseline
    %this time for anti trials
    antiMus = vertcat(popResCV2(latNeurons).antiMu);
    normLatMus = normToBaseInds(antiMus,baseInd);
    antiMus = vertcat(popResCV2(vermNeurons).antiMu);
    normVermMus = normToBaseInds(proMus,baseInd);
    
    
    sigStatsVerm = nan(size(normVermMus,1),2);
    sigStatsLat = nan(size(normLatMus,1),2);
   %numStds = 7; %number of stds away from baseline is counted as significant
    for neuronNum = 1:size(normVermMus,1)
        
         firstSigPoint = find(abs(normVermMus(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsVerm(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(normVermMus(neuronNum,3500:5500)));
       
         sigStatsVerm(neuronNum,2)=maxSigPoint-500;
        end
    end
    %now for lateral
    for neuronNum = 1:size(normLatMus,1)
        
         firstSigPoint = find(abs(normLatMus(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsLat(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(normLatMus(neuronNum,3500:5500)));
       
         sigStatsLat(neuronNum,2)=maxSigPoint-500;
        end
    end
    
    ccBar = repmat([0 0 1;1 0 0],2,1); %alternating r and b for bars
    vermMu = nanmean(sigStatsVerm);
    vermSem = nanstd(sigStatsVerm)./sqrt(size(vermPvA,1));
    latMu = nanmean(sigStatsLat);
    latSem = nanstd(sigStatsLat)./sqrt(size(latPvA,1));
    xIn = [vermMu' latMu']./1000; %convert to s
    erIn = [vermSem' latSem']./1000;
    hF = figure;
    hAx = axes('parent',hF);
   
   [hPatches] = createGroupedBarGraph(xIn,erIn,hAx,ccBar,0.3);
   suptitle('Time of CV2 differences From baseline, Anti Trials')
   set(hAx,'xtick',[1 2],'xticklabel',{'First Time','Maximum Time'})
   set(get(hAx,'ylabel'),'string','Time Relative to Saccade (s)')
     export_fig(hF,'-pdf','-opengl','-append',cvPopResPdfName)
    delete(hF)
    
   cv2PopTimeStats.sigStatsVermAnti = sigStatsVerm;
    cv2PopTimeStats.sigStatsLatAnti = sigStatsLat;
    cv2PopTimeStats.sigStatsVermPvA = sigStatsVermPvA;
    cv2PopTimeStats.sigStatsLatPvA = sigStatsLatPvA;
    save([resultsDir filesep 'cv2PopTimeStats2015.mat'],'cv2PopTimeStats')
end

%% DEFUNCT SUMMARy FIGURES
%     %below is old method for creating sumamry figures using loop
%      %now create figure showing the average CV2 trace over time for each
%      %neuron in 2x2 subplot with pro and anti and verm and lateral seperate
%      tV = intTimeVec;
%    
%     unitNum =1;
%     hFd = figure; %for showing the difference between pro and anti means
%     hAd(1) = subplot(1,2,1,'parent',hFd);
%     hAd(2) = subplot(1,2,2,'parent',hFd);
%     hF = figure;
%    hAxs =nan(1,4);
% 
%        hAxs(1) = subplot(2,2,1,'parent',hF);
%        hAxs(2) = subplot(2,2,2,'parent',hF);
%        nCount = 1; %for counting colours
%        for neuronNum = recordingsToAnalyse(vermNeurons)
%            if ~isempty(popResCV2(neuronNum,unitNum).antiMu) && ~isempty(popResCV2(neuronNum,unitNum).proMu)
%            line(tV',popResCV2(neuronNum,unitNum).antiMu,'parent',hAxs(2),'color',cc(mod(nCount,12)+1,:))
%           
%            line(tV',popResCV2(neuronNum,unitNum).proMu,'parent',hAxs(1),'color',cc(mod(nCount,12)+1,:))
%           %now get difference and plot on other graph
%           pMinusA = popResCV2(neuronNum,unitNum).proMu-popResCV2(neuronNum,unitNum).antiMu;
%            line(tV',pMinusA,'parent',hAd(1),'color',cc(mod(nCount,12)+1,:))
%            
%             nCount = nCount+1;
%            end
%        end
%        %vertcat(popResCV2(latNeurons).proMu) will give mean response of pro
%        set(get(hAxs(2),'title'),'string','Anti')
%        set(get(hAxs(1),'title'),'string','Pro')
%         set(get(hAxs(2),'ylabel'),'string','Vermis')
%         set(get(hAd(1),'title'),'string','Vermis')
%         
%        hAxs(3) = subplot(2,2,3,'parent',hF);
%        hAxs(4) = subplot(2,2,4,'parent',hF);
%        nCount = 1;
%        for neuronNum = recordingsToAnalyse(latNeurons)
%            if ~isempty(popResCV2(neuronNum,unitNum).antiMu) && ~isempty(popResCV2(neuronNum,unitNum).proMu)
%            
%            line(tV',popResCV2(neuronNum,unitNum).antiMu,'parent',hAxs(4),'color',cc(mod(nCount,12)+1,:))
%          
%            line(tV',popResCV2(neuronNum,unitNum).proMu,'parent',hAxs(3),'color',cc(mod(nCount,12)+1,:))
%             pMinusA = popResCV2(neuronNum,unitNum).proMu-popResCV2(neuronNum,unitNum).antiMu;
%            line(tV',pMinusA,'parent',hAd(2),'color',cc(mod(nCount,12)+1,:))
%             nCount = nCount+1;
%             end
%        end
%        set(get(hAd(2),'title'),'string','Lateral')
%        set(get(hAxs(4),'title'),'string','Anti')
%        set(get(hAxs(3),'title'),'string','Pro')
%         set(get(hAxs(4),'ylabel'),'string','Lateral')
%   % end
%     
%     temp = cell2mat(get(hAxs,'ylim'));
%     set(hAxs,'ylim',[0 2],'xlim',[-1 1])
%     
%     set(hAd,'xlim',[-1 1],'ylim',[-0.6 0.6])
%     suptitle(['CV2 Unit ' num2str(unitNum)])
%      export_fig(hF,'-pdf','-opengl','-append',cvPopResPdfName)
%         delete(hF)
%         
        
        
        %now exactly the same but this time normlise all to baseline firing rate.
%baseline here will be taken as ?[-0.4 0] to trialStart in kunimatsu so
%[-1.25 -0.75] relative to saccade. 
%now do a 2x2 subplot with Verm left, Lat right, Pro Top Anti Bottom
% baseInd = [3000 3250];
%  nFunc = @(x) ((x-mean(x(baseInd(1):baseInd(2))))./std(x(baseInd(1):baseInd(2))));
% 
% normCV2DispRange = [-30 30;0 10];
% 
%  hFd = figure; %for showing the difference between pro and anti means
%     hAd(1) = subplot(1,2,1,'parent',hFd);
%     hAd(2) = subplot(1,2,2,'parent',hFd);
%    hF = figure;
%    hAxs =nan(1,4)
% 
%        hAxs(1) = subplot(2,2,1,'parent',hF);
%        hAxs(2) = subplot(2,2,2,'parent',hF);
%        nCount = 1;
%        for neuronNum = recordingsToAnalyse(vermNeurons)
%            if ~isempty(popResCV2(neuronNum,unitNum).antiMu) && ~isempty(popResCV2(neuronNum,unitNum).proMu)
%                
%                
%                
%                line(tV',nFunc(popResCV2(neuronNum,unitNum).antiMu),'parent',hAxs(2),'color',cc(mod(nCount,12)+1,:))
%                
%                line(tV',nFunc(popResCV2(neuronNum,unitNum).proMu),'parent',hAxs(1),'color',cc(mod(nCount,12)+1,:))
%                pMinusA = nFunc(popResCV2(neuronNum,unitNum).proMu)-nFunc(popResCV2(neuronNum,unitNum).antiMu);
%                line(tV',pMinusA,'parent',hAd(1),'color',cc(mod(nCount,12)+1,:))
%                
%                nCount = nCount+1;
%            end
%        end
%        set(get(hAxs(2),'title'),'string','Anti')
%        set(get(hAxs(1),'title'),'string','Pro')
%        set(get(hAxs(2),'ylabel'),'string','Vermis')
%        
%        
%        hAxs(3) = subplot(2,2,3,'parent',hF);
%        hAxs(4) = subplot(2,2,4,'parent',hF);
%        nCount = 1;
%        for neuronNum = recordingsToAnalyse(latNeurons)
%            if ~isempty(popResCV2(neuronNum,unitNum).antiMu) && ~isempty(popResCV2(neuronNum,unitNum).proMu)
%                line(tV',nFunc(popResCV2(neuronNum,unitNum).antiMu),'parent',hAxs(3),'color',cc(mod(nCount,12)+1,:))
%                
%                line(tV',nFunc(popResCV2(neuronNum,unitNum).proMu),'parent',hAxs(4),'color',cc(mod(nCount,12)+1,:))
%                pMinusA = nFunc(popResCV2(neuronNum,unitNum).proMu)-nFunc(popResCV2(neuronNum,unitNum).antiMu);
%                line(tV',pMinusA,'parent',hAd(1),'color',cc(mod(nCount,12)+1,:))
%                nCount = nCount+1;
%            end
%        end
%        
%           set(get(hAxs(4),'title'),'string','Anti')
%        set(get(hAxs(3),'title'),'string','Pro')
%         set(get(hAxs(4),'ylabel'),'string','Lateral')
%   % end
%     
%     temp = cell2mat(get(hAxs,'ylim'));
%     set(hAxs,'ylim',normCV2DispRange(unitNum,:))
%     set(hAxs,'xlim',[-1 1])
%     suptitle(['CV2 Unit ' num2str(unitNum)])
%     export_fig(hF,'-pdf',expRend,'-append',cvPopResPdfName)
%        delete(hF) 
%% POPULATION RESPONSE
%load end point error data
    epErs = load([resultsDir filesep endPointErrorsFileName]);
 
numKinBins = 50;
axisTicks = [3250 4000 4750];
axTickLabs = {'-0.75','0','0.75'};
expRend = '-painters'; %for rendering heatmaps
validSaccadeLimits.reactionTime = [0.05 0.6]; %s %[min max] allowed
    validSaccadeLimits.saccadeDuration = [10 200]; %ms
    validSaccadeLimits.saccadeAmplitude = [1 20]; %deg
    validSaccadeLimits.saccadePeakVel = [350 450]; %deg/s
normFrDispRange = [-30 30;0 10];
frDispRange = [0 200;0 5;0 200;0 200;0 5]; %simple spikes then complex spikes for graph lims

if doPop
    sortByString = 'saccadeAngle';
   % sortByString = 'saccadeAmplitude';
   % sortByString = 'saccadePeakVel';
    %sortByString = 'saccadeDuration';
    %sortByString = 'reactionTime'; %TODO place loop here
    %if exist(popResFileName,'file')~=2
    stInit = cell(max(recordingsToAnalyse),2);
    popRes = struct('alignedTrials',stInit,'antiMu',stInit,'proMu',stInit);
    for neuronNum = recordingsToAnalyse;
        thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
        unitCount =1;
        for unitNum  = thisNeuronUnits;
            %get the trials and the numbers of stable ones
            theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
            theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
            
            theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
            theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
            
            %align all trils to the saccade time
            [alignedTrials.Anti antiMu] = alignAndExtractRates(theseStableTrials,theseAntiTrialNums,'saccadeTime','unitNum',unitNum);
            [alignedTrials.Pro proMu] = alignAndExtractRates(theseStableTrials,theseProTrialNums,'saccadeTime','unitNum',unitNum);
            
            %now remove any trials which are outside our valid saccade
            %parameters
            [cleanAlignedTrials.Anti rejectionCodesAnti percentTrialsRemoved.Anti] = trimBadTrials(alignedTrials.Anti,validSaccadeLimits);
            [cleanAlignedTrials.Pro rejectionCodesPro percentTrialsRemoved.Pro] = trimBadTrials(alignedTrials.Pro,validSaccadeLimits);
            
            %now copy clean trials over top of original for use in rest of
            %script
            alignedTrials = cleanAlignedTrials;
            
            %now get the mean on sem of sdf for pro and anti and plot on
            %top of each other, then export as png for the summary figures
            [hF] = createConfLimSdf(alignedTrials,[0.5 0.5]);
            %now export as png for loading into summary figs
            thisPngName = [resultsDir filesep neuronList(neuronNum).neuronName char(muSdfPngNameEnds{unitCount})];
             export_fig(thisPngName, '-png', hF);
            delete(hF)
            
            %get csOn dir to label plot
            peakDir = epErs.allNeuronsCsOndDir(neuronNum);
 csOnAngAng = plAngCode(peakDir,2);
            
            %now go over pro and anti and create seperate heat maps
            for dirLabelNum = 1:2;
                thisLabel = byDirLabel{dirLabelNum,2};
                [hA] = createSdfHeatMap(alignedTrials.(thisLabel),sortByString,[],'frRange',frDispRange(unitNum,:));
                set(get(hA,'title'),'string',[neuronList(neuronNum).neuronName ' unit ' num2str(unitNum) ' ' thisLabel ' trials by ' sortByString ', CS-on:' num2str(csOnAngAng)])
                %set(hA,'xlim',[3500 4500],'xtick',[3500 4000 4500],'xticklabel',{'-0.5','0','0.5'})
       set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
       
                try
                    export_fig(gcf,'-pdf','-append',expRend,popResPdfName)
                catch
                    export_fig(gcf,'-pdf',expRend,popResPdfName)
                end
            end
            close all
            popRes(neuronNum,unitCount).alignedTrials = alignedTrials;
            popRes(neuronNum,unitCount).antiMu = antiMu;
            popRes(neuronNum,unitCount).proMu = proMu;
            unitCount= unitCount+1;
        end
    end
    
   % save([resultsDir filesep sortByString popResFileName],'popRes','-v7.3')
    
    %else
    %   load([resultsDir filesep popResFileName]);
    %end
    %
    %
    %need to stick all pro trials across all neuorns
    
    
    for unitNum = 1:2
        numNeurons = max(recordingsToAnalyse);
        tempProTrials = cell(numNeurons,1);
        tempAntiTrials = cell(numNeurons,1);
        %                 for neuronNum = recordingsToAnalyse
        %                     %remove the lfpEoch filed as not present in all analysed
        %                     %recordings
        %                     if ~isempty(popRes(neuronNum,unitNum).alignedTrials)
        %                         tempRes = popRes(neuronNum,unitNum).alignedTrials.Pro;
        %                         if isfield(tempRes,'lfpEpoch')
        %                             tempProTrials{neuronNum} = rmfield(tempRes,'lfpEpoch');
        %
        %                         else
        %                             tempProTrials{neuronNum} =tempRes;
        %                         end
        %                         tempRes = popRes(neuronNum,unitNum).alignedTrials.Anti;
        %                         if isfield(tempRes,'lfpEpoch')
        %                             tempAntiTrials{neuronNum} = rmfield(tempRes,'lfpEpoch');
        %                         else
        %                             tempAntiTrials{neuronNum} = tempRes;
        %                         end
        %
        %                     end
        %
        %                 end
        for neuronNum = recordingsToAnalyse
            
            if ~isempty(popRes(neuronNum,unitNum).alignedTrials)
                tempProTrials{neuronNum} = popRes(neuronNum,unitNum).alignedTrials.Pro;
                
                tempAntiTrials{neuronNum} = popRes(neuronNum,unitNum).alignedTrials.Anti;
                
            end
            
        end
     
        
        allProTrials = [tempProTrials{:}];
        allAntiTrials = [tempAntiTrials{:}];
        [hA] = createSdfHeatMap(allProTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
        %set(hA,'xlim',[3500 4500],'xtick',[3500 4000 4500],'xticklabel',{'-0.5','0','0.5'})
       set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
        figure(get(hA,'parent'))
        suptitle(['All Neurons Pro, unit ' num2str(unitNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popResPdfName)
        delete(gcf)
        [hA] = createSdfHeatMap(allAntiTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
       set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
        figure(get(hA,'parent'))
        suptitle(['All Neurons Anti, unit ' num2str(unitNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popResPdfName)
        delete(gcf)
        %align each trial to saccade time and get psth with bin of 1ms smoothed with gaussian of 10ms
        
        
        %now do exactly the sma ebut split inot lateral and vermal neurons
        tempVermProTrials = tempProTrials(vermNeurons);
        tempVermAntiTrials = tempAntiTrials(vermNeurons);
        allVermProTrials = [tempVermProTrials{:}];
        allVermAntiTrials = [tempVermAntiTrials{:}];
        
        [hA] = createSdfHeatMap(allVermProTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
        figure(get(hA,'parent'))
       set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
        suptitle(['Verm Neurons Pro, unit ' num2str(unitNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popResPdfName)
        delete(gcf)
        [hA] = createSdfHeatMap(allVermAntiTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
        set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
        figure(get(hA,'parent'))
        suptitle(['Verm Neurons Anti, unit ' num2str(unitNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popResPdfName)
        delete(gcf)
        
        
        if ~isempty(latNeurons)
        tempLatProTrials = tempProTrials(latNeurons);
        tempLatAntiTrials = tempAntiTrials(latNeurons);
        allLatProTrials = [tempLatProTrials{:}];
        allLatAntiTrials = [tempLatAntiTrials{:}];
        
        [hA] = createSdfHeatMap(allLatProTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
        set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
        figure(get(hA,'parent'))
        suptitle(['Lat Neurons Pro, unit ' num2str(unitNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popResPdfName)
        delete(gcf)
        [hA] = createSdfHeatMap(allLatAntiTrials,sortByString,[],'numBins',numKinBins,'frRange',frDispRange(unitNum,:));
        set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
        figure(get(hA,'parent'))
        suptitle(['Lat Neurons Anti, unit ' num2str(unitNum)])
        export_fig(gcf,'-pdf',expRend,'-append',popResPdfName)
        delete(gcf)
        end
   end     
   
    %  numIters=10;
    %  for randIter = 1:numIters
    %      randNeurons = randperm(length(recordingsToAnalyse),10);
    %     if randIter==numIters
    %         %randNeurons = latNeurons;
    %         randNeurons = randperm(19,10)+20;
    %     end
    %      tempRanProTrials = tempProTrials(recordingsToAnalyse(randNeurons));
    %      tempRanAntiTrials = tempAntiTrials(recordingsToAnalyse(randNeurons));
    %      allRanProTrials = [tempRanProTrials{:}];
    %      allRanAntiTrials = [tempRanAntiTrials{:}];
    %
    %      [hA] = createSdfHeatMap(allRanProTrials,sortByString,[],'numBins',100,'frRange',[0 200]);
    %      figure(get(hA,'parent'))
    %      suptitle(['iter ' num2str(randIter) ' Neurons Pro'])
    %      export_fig(gcf,'-pdf','-opengl','-append',popResPdfName)
    %      delete(gcf)
    %      [hA] = createSdfHeatMap(allRanAntiTrials,sortByString,[],'numBins',100,'frRange',[0 200]);
    %      figure(get(hA,'parent'))
    %      suptitle(['iter ' num2str(randIter) ' Neurons Anti'])
    %      export_fig(gcf,'-pdf','-opengl','-append',popResPdfName)
    %      delete(gcf)
    %
    %  end
    
    
    
    
    %now do a 2x2 subplot with Verm left, Lat right, Pro Top Anti Bottom
    bothUnitVermPvA = cell(1,2);
    bothUnitLatPvA = cell(1,2);
    for unitNum = 1:2 %doing unit 2 overwreites vermPvA
        %method below was replaced by the non loop version further below
        %         hF = figure;
        %         hAxs =nan(1,4)
        %         %for dirLabelNum = 1:2
        %         hAxs(1) = subplot(2,2,1);
        %         hAxs(2) = subplot(2,2,2);
        %         for neuronNum = recordingsToAnalyse(vermNeurons)
        %             if ~isempty(popRes(neuronNum,unitNum).alignedTrials)
        %                 tV = linspace(0,4,4001);
        %                 line(tV',popRes(neuronNum,unitNum).antiMu,'parent',hAxs(2),'color',cc(ceil(rand(1)*12),:))
        %                 title('Anti')
        %                 ylabel('Vermis')
        %                 line(tV',popRes(neuronNum,unitNum).proMu,'parent',hAxs(1),'color',cc(ceil(rand(1)*12),:))
        %                 title('Pro')
        %             end
        %         end
        %
        %         hAxs(3) = subplot(2,2,3);
        %         hAxs(4) = subplot(2,2,4);
        %         for neuronNum = recordingsToAnalyse(latNeurons)
        %             if ~isempty(popRes(neuronNum,unitNum).alignedTrials)
        %                 tV = linspace(0,4,4001);
        %                 line(tV',popRes(neuronNum,unitNum).antiMu,'parent',hAxs(3),'color',cc(ceil(rand(1)*12),:))
        %                 title('Anti')
        %                 ylabel('Lateral')
        %                 line(tV',popRes(neuronNum,unitNum).proMu,'parent',hAxs(4),'color',cc(ceil(rand(1)*12),:))
        %                 title('Pro')
        %             end
        %         end
        %         % end
        %
        %         temp = cell2mat(get(hAxs,'ylim'));
        %         set(hAxs,'ylim',frDispRange(unitNum,:))
        %         set(hAxs,'xlim',[1 3])
        %         suptitle(['Unit ' num2str(unitNum)])
        %         export_fig(hF,'-pdf','-opengl','-append',popResPdfName)
        %         delete(hF)
        
        
        
        tV = intTimeVec; %get the itnerpolated time vector for the cv2 traces
        hFd = figure; %for showing the difference between pro and anti means
        hAd(1) = subplot(1,2,1,'parent',hFd);
        hAd(2) = subplot(1,2,2,'parent',hFd);
        hF = figure;
        hAxs =nan(1,4);
        
        %first plot the mean response to all pro trials for vermis, then
        %lateral, on a seperate figure also plot the difference between the two
        hAxs(1) = subplot(2,2,1,'parent',hF);
        hAxs(2) = subplot(2,2,2,'parent',hF);
        
        proMus = [popRes(vermNeurons,unitNum).proMu]';
        antiMus =  [popRes(vermNeurons,unitNum).antiMu]';
        pMinusA = proMus-antiMus;
        
        line(tV',antiMus,'parent',hAxs(2))
        
        line(tV',proMus,'parent',hAxs(1))
        line(tV',pMinusA,'parent',hAd(1))
        
        %now also plot mean of this as thick black line
        line(tV',nanmean(antiMus),'parent',hAxs(2),'color','k','linewidth',2)
        
        line(tV',nanmean(proMus),'parent',hAxs(1),'color','k','linewidth',2)
        line(tV',nanmean(pMinusA),'parent',hAd(1),'color','k','linewidth',2)
        set(get(hAxs(2),'title'),'string','Anti')
        set(get(hAxs(1),'title'),'string','Pro')
        set(get(hAxs(2),'ylabel'),'string','Vermis')
        set(get(hAd(1),'title'),'string','Vermis')
        
        
        hAxs(3) = subplot(2,2,3,'parent',hF);
        hAxs(4) = subplot(2,2,4,'parent',hF);
        
        if ~isempty(latNeurons)
        proMus = [popRes(latNeurons,unitNum).proMu]';
        antiMus =  [popRes(latNeurons,unitNum).antiMu]';
        pMinusA = proMus-antiMus;
        
        line(tV',antiMus,'parent',hAxs(4))
        
        line(tV',proMus,'parent',hAxs(3))
        line(tV',pMinusA,'parent',hAd(2))
        
        %now also plot mean of this as thick black line
        line(tV',nanmean(antiMus),'parent',hAxs(4),'color','k','linewidth',2)
        
        line(tV',nanmean(proMus),'parent',hAxs(3),'color','k','linewidth',2)
        line(tV',nanmean(pMinusA),'parent',hAd(2),'color','k','linewidth',2)
        
        set(get(hAxs(4),'title'),'string','Anti')
        set(get(hAxs(3),'title'),'string','Pro')
        set(get(hAxs(4),'ylabel'),'string','Lateral')
        set(get(hAd(2),'title'),'string','Lateral')
        temp = cell2mat(get(hAxs,'ylim'));
         end
        set(hAxs,'ylim',frDispRange(unitNum,:),'xlim',[-1 1])
        
        
        suptitle(['Unit ' num2str(unitNum)])
        export_fig(hF,'-pdf','-opengl','-append',popResPdfName)
        delete(hF)
        set(hAd,'xlim',[-1 1],'ylim',[-8 8])
        suptitle(['P-A Unit ' num2str(unitNum)])
        
       
        export_fig(hFd,'-pdf','-opengl','-append',popResPdfName)
        delete(hFd)
        
        
        %normalised to abseline version
        hFd = figure; %for showing the difference between pro and anti means
        hAd(1) = subplot(1,2,1,'parent',hFd);
        hAd(2) = subplot(1,2,2,'parent',hFd);
        hF = figure;
        hAxs =nan(1,4);
        
        %first plot the mean response to all pro trials for vermis, then
        %lateral, on a seperate figure also plot the difference between the two
        hAxs(1) = subplot(2,2,1,'parent',hF);
        hAxs(2) = subplot(2,2,2,'parent',hF);
        
        proMus = [popRes(vermNeurons,unitNum).proMu]';
        antiMus =  [popRes(vermNeurons,unitNum).antiMu]';
        
        
        normProMus = normToBaseInds(proMus,baseInd);
        normAntiNums = normToBaseInds(antiMus,baseInd);
        pMinusA = normProMus-normAntiNums;
        vermPvA = pMinusA;
        line(tV',normAntiNums,'parent',hAxs(2))
        
        line(tV',normProMus,'parent',hAxs(1))
        line(tV',pMinusA,'parent',hAd(1))
        
        %now also plot mean of this as thick black line
        line(tV',nanmean(normAntiNums),'parent',hAxs(2),'color','k','linewidth',2)
        
        line(tV',nanmean(normProMus),'parent',hAxs(1),'color','k','linewidth',2)
        line(tV',nanmean(pMinusA),'parent',hAd(1),'color','k','linewidth',2)
        set(get(hAxs(2),'title'),'string','Anti')
        set(get(hAxs(1),'title'),'string','Pro')
        set(get(hAxs(2),'ylabel'),'string','Vermis')
        set(get(hAd(1),'title'),'string','Vermis')
        
        
        hAxs(3) = subplot(2,2,3,'parent',hF);
        hAxs(4) = subplot(2,2,4,'parent',hF);
        if ~isempty(latNeurons)
        proMus = [popRes(latNeurons,unitNum).proMu]';
        antiMus =  [popRes(latNeurons,unitNum).antiMu]';
        normProMus = normToBaseInds(proMus,baseInd);
        normAntiNums = normToBaseInds(antiMus,baseInd);
        pMinusA = normProMus-normAntiNums;
        latPvA = pMinusA;
        
        line(tV',normAntiNums,'parent',hAxs(4))
        
        line(tV',normProMus,'parent',hAxs(3))
        line(tV',pMinusA,'parent',hAd(2))
        
        %now also plot mean of this as thick black line
        line(tV',nanmean(normAntiNums),'parent',hAxs(4),'color','k','linewidth',2)
        
        line(tV',nanmean(normProMus),'parent',hAxs(3),'color','k','linewidth',2)
        line(tV',nanmean(pMinusA),'parent',hAd(2),'color','k','linewidth',2)
        
        set(get(hAxs(4),'title'),'string','Anti')
        set(get(hAxs(3),'title'),'string','Pro')
        set(get(hAxs(4),'ylabel'),'string','Lateral')
        set(get(hAd(2),'title'),'string','Lateral')
        temp = cell2mat(get(hAxs,'ylim'));
        end
        set(hAxs,'ylim',normFrDispRange(unitNum,:),'xlim',[-1 1])
        
        
        suptitle(['Norm Unit ' num2str(unitNum)])
        export_fig(hF,'-pdf','-opengl','-append',popResPdfName)
        delete(hF)
        set(hAd,'xlim',[-1 1],'ylim',[-30 30],'ytick',[-30 0 30])
        suptitle(['Norm P-A Unit ' num2str(unitNum)])
        export_fig(hFd,'-pdf','-opengl','-append',popResPdfName)
        delete(hFd)
        
        bothUnitVermPvA{unitNum} = vermPvA;
        
        bothUnitLatPvA{unitNum} = latPvA;
    end
    
    
    unitNum = 1; %only for SS so far
    %TIMINGS OF SIGNIFICANT POINTS
    %now neuron by neuron we should find the first time that each is abs>4,
    %also find the maximum
    vermPvA = bothUnitVermPvA{unitNum};
    latPvA = bothUnitLatPvA{unitNum};
    sigStatsVermPvA = nan(size(vermPvA,1),2);
    sigStatsLatPvA = nan(size(latPvA,1),2);
   numStds = 5; %number of stds away from baseline is counted as significant
    for neuronNum = 1:size(vermPvA,1)
        
         firstSigPoint = find(abs(vermPvA(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsVermPvA(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(vermPvA(neuronNum,3500:5500)));
       
         sigStatsVermPvA(neuronNum,2)=maxSigPoint-500;
        end
    end
    %now for lateral
    for neuronNum = 1:size(latPvA,1)
        
         firstSigPoint = find(abs(latPvA(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsLatPvA(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(latPvA(neuronNum,3500:5500)));
       
         sigStatsLatPvA(neuronNum,2)=maxSigPoint-500;
        end
    end
    
    ccBar = repmat([0 0 1;1 0 0],2,1); %alternating r and b for bars
    vermMu = nanmean(sigStatsVermPvA);
    vermSem = nanstd(sigStatsVermPvA)./sqrt(size(vermPvA,1));
    latMu = nanmean(sigStatsLatPvA);
    latSem = nanstd(sigStatsLatPvA)./sqrt(size(latPvA,1));
    xIn = [vermMu' latMu']./1000; %convert to s
    erIn = [vermSem' latSem']./1000;
    hF = figure;
    hAx = axes('parent',hF);
   
   [hPatches] = createGroupedBarGraph(xIn,erIn,hAx,ccBar,0.3);
   suptitle('Time of Rate differences PvA')
   set(hAx,'xtick',[1 2],'xticklabel',{'First Time','Maximum Time'})
   set(get(hAx,'ylabel'),'string','Time Relative to Saccade (s)')
     export_fig(hF,'-pdf','-opengl','-append',popResPdfName)
    delete(hF)
    

    %now dot he same but instead of looking at pro vs anti differences
    %let's just get the times of modualtion from normalised baseline
    proMus = [popRes(latNeurons,unitNum).proMu]';
    normLatMus = normToBaseInds(proMus,baseInd);
    proMus = [popRes(vermNeurons,unitNum).proMu]';
  
    normVermMus = normToBaseInds(proMus,baseInd);
    
    
    sigStatsVerm = nan(size(normVermMus,1),2);
    sigStatsLat = nan(size(latPvA,1),2);
   %numStds = 7; %number of stds away from baseline is counted as significant
    for neuronNum = 1:size(normVermMus,1)
        
         firstSigPoint = find(abs(normVermMus(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsVerm(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(normVermMus(neuronNum,3500:5500)));
       
         sigStatsVerm(neuronNum,2)=maxSigPoint-500;
        end
    end
    %now for lateral
    for neuronNum = 1:size(normLatMus,1)
        
         firstSigPoint = find(abs(normLatMus(neuronNum,3500:5500))>numStds,1,'first') -500;
        if~isempty(firstSigPoint)
         sigStatsLat(neuronNum,1) = firstSigPoint;
         [val maxSigPoint] = max(abs(normLatMus(neuronNum,3500:5500)));
       
         sigStatsLat(neuronNum,2)=maxSigPoint-500;
        end
    end
    
    ccBar = repmat([0 0 1;1 0 0],2,1); %alternating r and b for bars
    vermMu = nanmean(sigStatsVerm);
    vermSem = nanstd(sigStatsVerm)./sqrt(size(vermPvA,1));
    latMu = nanmean(sigStatsLat);
    latSem = nanstd(sigStatsLat)./sqrt(size(latPvA,1));
    xIn = [vermMu' latMu']./1000; %convert to s
    erIn = [vermSem' latSem']./1000;
    hF = figure;
    hAx = axes('parent',hF);
   
   [hPatches] = createGroupedBarGraph(xIn,erIn,hAx,ccBar,0.3);
   suptitle('Time of Rate differences From baseline')
   set(hAx,'xtick',[1 2],'xticklabel',{'First Time','Maximum Time'})
   set(get(hAx,'ylabel'),'string','Time Relative to Saccade (s)')
     export_fig(hF,'-pdf','-opengl','-append',popResPdfName)
    delete(hF)
    
    
    
     ratePopTimeStats.sigStatsVerm = sigStatsVerm;
    ratePopTimeStats.sigStatsLat = sigStatsLat;
    ratePopTimeStats.sigStatsVermPvA = sigStatsVermPvA;
    ratePopTimeStats.sigStatsLatPvA = sigStatsLatPvA;
    save([resultsDir filesep 'ratePopTimeStats2015final.mat'],'ratePopTimeStats')
    
    %now lets do a correlation between the times that CV2 and rate
    %modualted from baseline to see if they represent the same thing.
    
    %this doesn't seem to be used ever and causes an error with the complex
    %spikes with rate having one less entry than cv2
%     [r,p] = corrcoef(ratePopTimeStats.sigStatsVerm(:,1),cv2PopTimeStats.sigStatsVerm(:,1),'rows','complete')
%     [r,p] = corrcoef(ratePopTimeStats.sigStatsVerm(:,2),cv2PopTimeStats.sigStatsVerm(:,2),'rows','complete')
%     [r,p] = corrcoef(ratePopTimeStats.sigStatsVermPvA(:,1),cv2PopTimeStats.sigStatsVermPvA(:,1),'rows','complete')
%     [r,p] = corrcoef(ratePopTimeStats.sigStatsVermPvA(:,2),cv2PopTimeStats.sigStatsVermPvA(:,2),'rows','complete')
%    
    %
end
  a =1;
     
%% KUNIMATSU STYLE ANALYSIS
sigLevel= 0.01;
ccBar = repmat([0 0 1;1 0 0],3,1); %alternating r and b for bars
%TODO check why kunimatsu analysis produces 0's, should be nan?!?
distLims = [0 200; 0 5];
axLims = [0 150; 0 5]; %fitst for ss seond for cs rates
axLimsCV2 = [0 2];
%for spider plot angles
anglesToTest = linspace(0,2*pi,9);
anglesToTest = anglesToTest(1:8); %8 evenly spaced going counterclockwise
inAngleOrder = [6 3 2 1 4 7 8 9]; %order of trialConditionCodes(:,dirLabelNum) in terms of starting at angle 0 (horizontal right) and counting coutner clockwise
           
if exist(allKuniResFile,'file')~=2
    %semRates = cell(max(recordingsToAnalyse),2,3); %first colum mean pro rate, second anti
    stInit = cell(max(recordingsToAnalyse),2);
    allNeuronResults = struct('kuniRes',stInit,'kuniResultsStats',stInit,...
        'muRates',stInit,'baselineRates',stInit,'muCV2s',stInit);
    
    for unitNum = 1:2; %TODO can be used to do kunimatsu analysis on complex spikes as well
        %allNeuronResults uses 1 for spikes pikes 2 for complex, thisUnitNum is
        %required to get real spikes from trail structure whichincldues units 4
        %and 5
        muRates = cell(max(recordingsToAnalyse),3,3); %first colum mean pro rate, second anti, third if sig diff between p and a, 3rd dim is time point to be tested
        muCV2s = muRates;
        baselineRates = cell(max(recordingsToAnalyse),2); %first colum is simple spikes econd is CS
        baselineCV2s = baselineRates;
        ratePvals = baselineRates; %to hold the p values for each period's modualtion, first colum pro, second anti
        cv2Pvals = baselineRates;
        for neuronNum = recordingsToAnalyse;
            %get unit num for this unit for thsi recording
            if unitNum<=length(neuronList(neuronNum).fileList{1,3})
                thisUnitNum = neuronList(neuronNum).fileList{1,3}(unitNum);
                
                %build figure filename
                thisNeuronSumFigName = [resultsDir filesep neuronList(neuronNum).neuronName indNeuronSumFigEnd];
                %loop over all chosen rasters and create and export their figures
                for rastSelected = rastToPlot
                    
                    [kuniResults kuniResultsStats hAx] = kunimatsuAnalysis(wholeNeuronResults(neuronNum).allStableTrials,...
                        wholeNeuronResults(neuronNum).selectedTrials.(char(rastSelected)),timeWindow,'',thisUnitNum,'sigLevel',sigLevel);
                    set(hAx{1},'ylim',distLims(unitNum,:))
                    title([neuronList(neuronNum).neuronName ' ' rastSelected ' unit ' num2str(thisUnitNum)])
                    
                    export_fig(gcf,'-pdf','-append',thisNeuronSumFigName)
                    delete(gcf)
                    
                    allNeuronResults(neuronNum,unitNum).kuniResults.(char(rastSelected)) = kuniResults;
                    allNeuronResults(neuronNum,unitNum).kuniResultsStats.(char(rastSelected)) = kuniResultsStats;
                end
                baselineRates{neuronNum,unitNum} = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corSacTrials.baselineRates; %was previously only for pro trials but why not use all
                baselineCV2s{neuronNum,unitNum} = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corSacTrials.baselineCV2s;
                %keep a copy of teh pValues for each period being tested
                %against baseline
                ratePvals{neuronNum,1} = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corSacTrials.pValues;
                %ratePvals{neuronNum,2} = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.pValues;
                cv2Pvals{neuronNum,1} = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corSacTrials.pValuesCV2;
                %cv2Pvals{neuronNum,2} = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.pValuesCV2;
                
                
                %now get mean rates i matrix form
                muRates{neuronNum,1,1} = mean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.saccRates);
                muRates{neuronNum,2,1} = mean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.saccRates);
                muRates{neuronNum,1,2} = mean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.dirInstRates);
                muRates{neuronNum,2,2} = mean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.dirInstRates);
                muRates{neuronNum,1,3} = mean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.typeInstRates);
                muRates{neuronNum,2,3} = mean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.typeInstRates);
                
                
                muCV2s{neuronNum,1,1} = nanmean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.saccCV2s,2); %ph changed as nanmean now seems to need dimension even if singleton involved
                muCV2s{neuronNum,2,1} = nanmean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.saccCV2s,2);
                muCV2s{neuronNum,1,2} = nanmean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.dirInstCV2s,2);
                muCV2s{neuronNum,2,2} = nanmean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.dirInstCV2s,2);
                muCV2s{neuronNum,1,3} = nanmean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.typeInstCV2s,2);
                muCV2s{neuronNum,2,3} = nanmean(allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.typeInstCV2s,2);
                
                
                %TODO do a tuning plot of teh byDirection kunimatsu stats
                %this should also be done for CV2 to compare
                %now do by direction for both pro and anti seperately
                muRateByDir = nan(9,2); %to keep the saccade period firing rae for each direction
                semRateByDir = nan(9,2);
                muCV2ByDir = nan(9,2);
                semCV2ByDir = nan(9,2);
                for dirLabelNum = 1:2
                    
                    
                    for plNum = find(~isnan(trialConditionCodes(:,dirLabelNum)))'
                        
                        
                        
                        %first get all relevant trials
                        thisConditionCode = trialConditionCodes(plNum,dirLabelNum)
                        tempLogCell = cellfun(@(x) x==thisConditionCode,{wholeNeuronResults(neuronNum).allStableTrials.conditionCode},'uniformoutput',false);
                        theseDirTrialNums = find([tempLogCell{:}]);
                        corDirTrialNums = intersect(wholeNeuronResults(neuronNum).selectedTrials.corSacTrials,theseDirTrialNums);
                        
                        
                        
                        
                        %         %now get the firing rate in the saccade period using the kunimatsu
                        %         %anyslsis function
                        if ~isempty(corDirTrialNums)
                        [kuniResDir kuniResDirStats] = kunimatsuAnalysis(wholeNeuronResults(neuronNum).allStableTrials,...
                            corDirTrialNums,timeWindow,'',thisUnitNum,'sigLevel',sigLevel);
                        muRateByDir(plNum,dirLabelNum) = kuniResDirStats.muRates(4);
                        semRateByDir(plNum,dirLabelNum) = kuniResDirStats.semRates(4);
                        
                        muCV2ByDir(plNum,dirLabelNum) = kuniResDirStats.muRates(9);
                        semCV2ByDir(plNum,dirLabelNum) = kuniResDirStats.semRates(9);
                        
                        delete(gcf)
                        else
                            muRateByDir(plNum,dirLabelNum) = NaN;
                             semRateByDir(plNum,dirLabelNum)= NaN;
                             muCV2ByDir(plNum,dirLabelNum)= NaN;
                             semCV2ByDir(plNum,dirLabelNum) = NaN;
                        end
                    end
                end
                allNeuronResults(neuronNum,unitNum).kuniResults.muCV2ByDir = muCV2ByDir;
                allNeuronResults(neuronNum,unitNum).kuniResults.muRateByDir = muRateByDir;
                allNeuronResults(neuronNum,unitNum).kuniResults.semRateByDir = semRateByDir;
                allNeuronResults(neuronNum,unitNum).kuniResults.semCV2ByDir = semCV2ByDir;
%                 rateSpidLim = nanmax(muRateByDir(:));
%                 cv2SpidLim = nanmax(muCV2ByDir(:));
%                 for dirLabelNum = 1:2
%                     %now add spider plot to middle showing saccade period firing rate for
%                     %each direction
%                     theseMuRates = allNeuronResults(neuronNum,unitNum).kuniResults.muRateByDir(:,dirLabelNum);
%                     %hA = subaxis(3,3,5,'SpacingVert',0.2,'MR',0.05);
%                     
%                     ratesToPlot = theseMuRates(~isnan(theseMuRates));
%                     %spider starts at 3 o'clock and goes anti-clockwise
%                     %plNum starts at 10 oclock and goes clockwise
%                     %so flipud
%                     %TODO repalce with:
%                     %[hAx] = createConfLimSpider(datVal,datErr,anglesToTest,hAx)
%                     ratesToPlot = flipud(circshift(ratesToPlot,3));
%                     [f, ca, o] = createSpider(ratesToPlot ,'',[0 rateSpidLim],repmat({''},8,1),'');
%                     
%                     suptitle([neuronList(neuronNum).neuronName ' ' byDirLabel{dirLabelNum,2} ' unit ' num2str(thisUnitNum) ' Rate'])
%                     export_fig(f,'-pdf','-append',thisNeuronSumFigName)
%                     delete(f)
%                     
%                     
%                     
%                     theseMuCV2s = allNeuronResults(neuronNum,unitNum).kuniResults.muCV2ByDir(:,dirLabelNum);
%                     
%                     %ratesToPlot = theseMuCV2s(~isnan(theseMuCV2s)); %this line was designed to ignore the 5th entry which is always nan but with cv2 there are sometimes nan so isntead will harddcode index 5
%                     ratesToPlot = theseMuCV2s([1 2 3 4 6 7 8 9]);
%                    
%                         ratesToPlot = flipud(circshift(ratesToPlot,3));
%                         [f, ca, o] = createSpider(ratesToPlot ,'',[0 cv2SpidLim],repmat({''},8,1),'');
%                  
%                     suptitle([neuronList(neuronNum).neuronName ' ' byDirLabel{dirLabelNum,2} ' unit ' num2str(thisUnitNum) ' CV2'])
%                     export_fig(f,'-pdf','-append',thisNeuronSumFigName)
%                     delete(f)
%                 end
                
                
                %here put the conf lim spider approach whcih shows tuning
                %for both pro and anti on the same graph with SEM
                %confidence limits as shaded patches, this repalces the
                %seperate spider plot apporach above.
                datErr = {[semRateByDir(inAngleOrder,1)] [semRateByDir(inAngleOrder,2)]};
                datVal= {[muRateByDir(inAngleOrder,1)] [muRateByDir(inAngleOrder,2)]};
                hFs = figure;
                hAs = axes('parent',hFs);
                [hAx] = createConfLimSpider(datVal,datErr,anglesToTest,hAs,axLims(unitNum,:));
                set(get(hAx,'title'),'string',[neuronList(neuronNum).neuronName  ' Rate for unit ' num2str(thisUnitNum)])
                
                export_fig(hFs,'-pdf','-append',allNeuronSumFigName) %keep a copy in the main pdf and also export to png
               thisPngName = [resultsDir filesep neuronList(neuronNum).neuronName  ' Rate for unit ' num2str(thisUnitNum) 'tuning.png'];
               export_fig(hFs,'-png',thisPngName) 
               delete(hFs)
                datErr = {[semCV2ByDir(inAngleOrder,1)] [semCV2ByDir(inAngleOrder,2)]};
                datVal= {[muCV2ByDir(inAngleOrder,1)] [muCV2ByDir(inAngleOrder,2)]};
                hFs = figure;
                hAs = axes('parent',hFs);
                [hAx] = createConfLimSpider(datVal,datErr,anglesToTest,hAs,axLimsCV2);
                set(get(hAx,'title'),'string',[neuronList(neuronNum).neuronName  ' CV2 for unit ' num2str(thisUnitNum)])
                
                export_fig(hFs,'-pdf','-append',allNeuronSumFigName)
                 thisPngName = [resultsDir filesep neuronList(neuronNum).neuronName  ' CV2 for unit ' num2str(thisUnitNum) 'tuning.png'];
               export_fig(hFs,'-png',thisPngName) 
                delete(hFs)
                
                
                theseProTrials = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials;
                theseAntiTrials = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials;
                %allNeuronResults(neuronNum).proAntiSigSacc = ttest2(theseProTrials.saccRates,theseAntiTrials.saccRates)
                %allNeuronResults(neuronNum).proAntiSigTi = ttest2(theseProTrials.typeInstRates,theseAntiTrials.typeInstRates)
                %allNeuronResults(neuronNum).proAntiSigDi = ttest2(theseProTrials.dirInstRates,theseAntiTrials.dirInstRates)
                muRates{neuronNum,3,1} = ranksum(theseProTrials.saccRates,theseAntiTrials.saccRates);
                muRates{neuronNum,3,2} = ranksum(theseProTrials.dirInstRates,theseAntiTrials.dirInstRates);
                muRates{neuronNum,3,3} = ranksum(theseProTrials.typeInstRates,theseAntiTrials.typeInstRates);
                
                
                muCV2s(neuronNum,3,:) = {nan};
                
                if sum(~isnan(theseProTrials.saccCV2s))>1 && sum(~isnan(theseAntiTrials.saccCV2s))>1
                    muCV2s{neuronNum,3,1} = ranksum(theseProTrials.saccCV2s,theseAntiTrials.saccCV2s);
                end
                if sum(~isnan(theseProTrials.dirInstCV2s))>1 && sum(~isnan(theseAntiTrials.dirInstCV2s))>1
                    muCV2s{neuronNum,3,2} = ranksum(theseProTrials.dirInstCV2s,theseAntiTrials.dirInstCV2s);
                end
                if sum(~isnan(theseProTrials.typeInstCV2s))>1 && sum(~isnan(theseAntiTrials.typeInstCV2s))>1
                    muCV2s{neuronNum,3,3} = ranksum(theseProTrials.typeInstCV2s,theseAntiTrials.typeInstCV2s);
                end
                allNeuronResults(neuronNum,unitNum).muRates = muRates;
                allNeuronResults(neuronNum,unitNum).muCV2s = muCV2s;
                %allNeuronResults(neuronNum,thisUnitNum).baselineRates
            end
        end
        
        
        %for this unit we should get a count of the number of significantly
        %modualted neurons in each period and tehn the intersections
        %firs need lonear inds to access struct
        %         thisUnitVermNeurons = max(recordingsToAnalyse)*(unitNum-1)+vermNeurons
        %         thisUnitLatNeurons = max(recordingsToAnalyse)*(unitNum-1)+latNeurons
        
        %         thisUnitVermNeurons = sub2ind(size(allNeuronResults),vermNeurons,unitNum*ones(1,length(vermNeurons)));
        %         thisUnitLatNeurons = sub2ind(size(allNeuronResults),latNeurons,unitNum*ones(1,length(latNeurons)));
        %           ratePvals{neuronNum,1} = allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.pValues;
        
        selSigLevel = sigLevel;
        
        vermPs = vertcat(ratePvals{vermNeurons,1}); %all saccades now
        percChangeVerm = 100*sum(vermPs<selSigLevel)./length(vermNeurons);
        latPs = vertcat(ratePvals{latNeurons,1});
        if ~isempty(latPs)
        percChangeLat = 100*sum(latPs<selSigLevel)./length(latNeurons);
        else
            percChangeLat = nan(1,3);
        end
        
        hFb = figure;
        hAb = axes('parent',hFb);
        [hPatches] = createGroupedBarGraph([percChangeVerm; percChangeLat]',zeros(3,2),hAb,ccBar,0.3);
        set(hAb,'xtick',[1:3],'xticklabel',{'typeInst','dirInst','sacc'})
        ylim([0 100])
        ylabel('Percent Modulated')
        title(['Firing Rate Modulated, Unit ' num2str(unitNum)])
        export_fig(hFb,'-pdf','-append',allNeuronSumFigName)
        delete(hFb)
        
        vermPs = vertcat(cv2Pvals{vermNeurons,1}); %all saccades now
        percChangeVerm = 100*sum(vermPs<selSigLevel)./length(vermNeurons);
        if ~isempty(latNeurons)
        latPs = vertcat(cv2Pvals{latNeurons,1});
        percChangeLat = 100*sum(latPs<selSigLevel)./length(latNeurons);
        else
            percChangeLat = nan(1,3);
        end
        
        hFb = figure;
        hAb = axes('parent',hFb);
        [hPatches] = createGroupedBarGraph([percChangeVerm; percChangeLat]',zeros(3,2),hAb,ccBar,0.3);
        set(hAb,'xtick',[1:3],'xticklabel',{'typeInst','dirInst','sacc'})
        ylim([0 100])
        ylabel('Percentage Modulated')
        title(['CV2 Modulated, Unit ' num2str(unitNum)])
        export_fig(hFb,'-pdf','-append',allNeuronSumFigName)
        delete(hFb)
        
        
        
        %initialise summary figures
        %on these the mean firing rate of each neuron during sacade period is represented by a point on the scatter
        %x axis pro trial and y axis anti trials
        %coulored based on significant differences
        %now need to split into lateral and vermis neurons
        
        %label each poitn with neuronNum
         scatPointLabels = cellstr(num2str([1:size(recList,2)]')); %using actual labels is too big for graph so just numbers
    %TODO this should be fixed to recNumbers not just 1:n
        %TODO also get percentages of neurons in each period which are
        %different between pro and anti
         dispX = diff(axLims(unitNum,:))/20; %displacement so not oevr top of dots
      dispY = 0;
        
        testPeriods = {'Saccade Period','Type Instruction','Direction Instruction'};
        sigPercents = nan(3,2); %(testPeriod, [verm lat])hold the percentage of all neurons that are significantly different between the pro and anti trials
        %now loop over the 3 testing periods and porduce a figure for each
        for testPeriodNum = 1:3
            hFsacc = figure;
            hAsacc = axes('parent',hFsacc);
            hF2 = figure; %for lat and verm subplots
            
            emptyRecs = find(cellfun('isempty',muRates))
            muRates(emptyRecs) = {nan};
            
            proSaccadeMeanRates = [muRates{:,1,testPeriodNum}];
            antiSaccadeMeanRates = [muRates{:,2,testPeriodNum}];
            proAntiSig = [muRates{:,3,testPeriodNum}];
            
            set(hAsacc,'xlim',axLims(unitNum,:),'ylim',axLims(unitNum,:),'Nextplot','add')
            line(axLims(unitNum,:),axLims(unitNum,:),'color','k','parent',hAsacc); %LOI
            scatter(proSaccadeMeanRates,antiSaccadeMeanRates,20,'k','parent',hAsacc);
            %now scatter signifcant points with filled circles
            sigSacs = find(proAntiSig<selSigLevel);
            scatter(proSaccadeMeanRates(sigSacs),antiSaccadeMeanRates(sigSacs),20,'k','filled','parent',hAsacc)
            %now place the name of each cell aobve the point
           %TODO fix text(proSaccadeMeanRates+dispX,antiSaccadeMeanRates+dispY,scatPointLabels,'parent',hAsacc)
           
            set(get(hAsacc,'title'),'string',['All Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
            set(get(hAsacc,'xlabel'),'string','Pro Rate (Hz)')
            set(get(hAsacc,'ylabel'),'string','Anti Rate (Hz)')
            export_fig(hFsacc,'-pdf','-append',allNeuronSumFigName)
            delete(hFsacc)
            
            %now split into vermis and lateral cerebellar recordings
            sigVermSaccs = intersect(sigSacs,vermNeurons);
            sigLatSaccs = intersect(sigSacs,latNeurons);
            
            sigPercents(testPeriodNum,1) = 100*length(sigVermSaccs)/length(vermNeurons);
            sigPercents(testPeriodNum,2) = 100*length(sigLatSaccs)/length(latNeurons);
            
            hAx1 = subplot(1,2,1,'parent',hF2);
            set(hAx1,'xlim',axLims(unitNum,:),'ylim',axLims(unitNum,:),'Nextplot','add')
            line(axLims(unitNum,:),axLims(unitNum,:),'color','k','parent',hAx1); %LOI
            scatter(proSaccadeMeanRates(vermNeurons),antiSaccadeMeanRates(vermNeurons),20,'k','parent',hAx1);
            %now scatter signifcant points with filled circles
             scatter(proSaccadeMeanRates(sigVermSaccs),antiSaccadeMeanRates(sigVermSaccs),20,'k','filled','parent',hAx1);
              %now place the name of each cell aobve the point
            %TODO fix text(proSaccadeMeanRates(vermNeurons)+dispX,antiSaccadeMeanRates(vermNeurons)+dispY,scatPointLabels(vermNeurons),'parent',hAx1)
           
            set(get(hAx1,'title'),'string',['Vermis Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
            set(get(hAx1,'xlabel'),'string','Pro Rate (Hz)')
            set(get(hAx1,'ylabel'),'string','Anti Rate (Hz)')
            %export_fig(hF2,'-pdf','-append',allNeuronSumFigName)
            
            
            hAx2 = subplot(1,2,2,'parent',hF2);
            set(hAx2,'xlim',axLims(unitNum,:),'ylim',axLims(unitNum,:),'Nextplot','add')
            line(axLims(unitNum,:),axLims(unitNum,:),'color','k','parent',hAx2); %LOI
            scatter(proSaccadeMeanRates(latNeurons),antiSaccadeMeanRates(latNeurons),20,'k','parent',hAx2);
            %now scatter signifcant points with filled circles
              scatter(proSaccadeMeanRates(sigLatSaccs),antiSaccadeMeanRates(sigLatSaccs),20,'k','filled','parent',hAx2)
             %now place the name of each cell aobve the point
       %TODO fix     text(proSaccadeMeanRates(latNeurons)+dispX,antiSaccadeMeanRates(latNeurons)+dispY,scatPointLabels(latNeurons),'parent',hAx2)
           
            set(get(hAx2,'title'),'string',['Lateral Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
            set(get(hAx2,'xlabel'),'string','Pro Rate (Hz)')
            set(get(hAx2,'ylabel'),'string','Anti Rate (Hz)')
            export_fig(hF2,'-pdf','-append',allNeuronSumFigName)
            
            delete(hF2)
            
        end
        
        
       
       %sigPercents needs reordering to match the labels as is currently in
       %test period order
       sigPercents2 = sigPercents([2 3 1],:);
        hFb = figure;
        hAb = axes('parent',hFb);
        [hPatches] = createGroupedBarGraph(sigPercents2,zeros(3,2),hAb,ccBar,0.3);
        set(hAb,'xtick',[1:3],'xticklabel',{'typeInst','dirInst','sacc'})
        ylim([0 100])
        ylabel('Percent Modulated')
        title(['Firing Rate Different PvA, Unit ' num2str(unitNum)])
        export_fig(hFb,'-pdf','-append',allNeuronSumFigName)
        delete(hFb)
        
        
        %TODO replace with scatter with bidirectional error bars
        %repeated for CV2 values
        %testPeriods = {'Saccade Period','Type Instruction','Direction Instruction'};
          sigPercents = nan(3,2); %(testPeriod, [verm lat])hold the percentage of all neurons that are significantly different between the pro and anti trials
       dispX = 0.01; %displacement so not oevr top of dots
      dispY = 0;
        %now loop over the 3 testing periods and porduce a figure for each
        for testPeriodNum = 1:3
            hFsacc = figure;
            hAsacc = axes('parent',hFsacc);
            hF2 = figure; %for lat and verm subplots
            
            emptyRecs = find(cellfun('isempty',muCV2s))
            muCV2s(emptyRecs) = {nan};
            
            proSaccadeMeanRates = [muCV2s{:,1,testPeriodNum}];
            antiSaccadeMeanRates = [muCV2s{:,2,testPeriodNum}];
            proAntiSig = [muCV2s{:,3,testPeriodNum}];
            
            
            
            set(hAsacc,'xlim',axLimsCV2,'ylim',axLimsCV2,'Nextplot','add')
            line(axLimsCV2,axLimsCV2,'color','k','parent',hAsacc); %LOI
            scatter(proSaccadeMeanRates,antiSaccadeMeanRates,20,'k','parent',hAsacc);
            %now scatter signifcant points with filled circles
            sigSacs = find(proAntiSig<selSigLevel);
            scatter(proSaccadeMeanRates(sigSacs),antiSaccadeMeanRates(sigSacs),20,'k','filled','parent',hAsacc)
            
            %now place the name of each cell aobve the point
         %TODO fix   text(proSaccadeMeanRates+dispX,antiSaccadeMeanRates+dispY,scatPointLabels,'parent',hAsacc)
            set(get(hAsacc,'title'),'string',['All Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
            set(get(hAsacc,'xlabel'),'string','Pro CV2')
            set(get(hAsacc,'ylabel'),'string','Anti CV2')
            export_fig(hFsacc,'-pdf','-append',allNeuronSumFigName)
            
            delete(hFsacc)
            %now split into vermis and lateral cerebellar recordings
            sigVermSaccs = intersect(sigSacs,vermNeurons);
            sigLatSaccs = intersect(sigSacs,latNeurons);
            
            sigPercents(testPeriodNum,1) = 100*length(sigVermSaccs)/length(vermNeurons);
            sigPercents(testPeriodNum,2) = 100*length(sigLatSaccs)/length(latNeurons);
            
            hAx1 = subplot(1,2,1,'parent',hF2);
            set(hAx1,'xlim',axLimsCV2,'ylim',axLimsCV2,'Nextplot','add')
            line(axLimsCV2,axLimsCV2,'color','k','parent',hAx1); %LOI
            scatter(proSaccadeMeanRates(vermNeurons),antiSaccadeMeanRates(vermNeurons),20,'k','parent',hAx1);
            %now scatter signifcant points with filled circles
          
            scatter(proSaccadeMeanRates(sigVermSaccs),antiSaccadeMeanRates(sigVermSaccs),20,'k','filled','parent',hAx1);
             %now place the name of each cell aobve the point
           %TODO fix text(proSaccadeMeanRates(vermNeurons)+dispX,antiSaccadeMeanRates(vermNeurons)+dispY,scatPointLabels(vermNeurons),'parent',hAx1)
           
            set(get(hAx1,'title'),'string',['Vermis Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
            set(get(hAx1,'xlabel'),'string','Pro CV2')
            set(get(hAx1,'ylabel'),'string','Anti CV2')
            %export_fig(hF2,'-pdf','-append',allNeuronSumFigName)
            
            
            hAx2 = subplot(1,2,2,'parent',hF2);
            set(hAx2,'xlim',axLimsCV2,'ylim',axLimsCV2,'Nextplot','add')
            line(axLimsCV2,axLimsCV2,'color','k','parent',hAx2); %LOI
            scatter(proSaccadeMeanRates(latNeurons),antiSaccadeMeanRates(latNeurons),20,'k','parent',hAx2);
            %now scatter signifcant points with filled circles
             scatter(proSaccadeMeanRates(sigLatSaccs),antiSaccadeMeanRates(sigLatSaccs),20,'k','filled','parent',hAx2)
            %now place the name of each cell aobve the point
         %TODO fix   text(proSaccadeMeanRates(latNeurons)+dispX,antiSaccadeMeanRates(latNeurons)+dispY,scatPointLabels(latNeurons),'parent',hAx2)
           
            set(get(hAx2,'title'),'string',['Lateral Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
            set(get(hAx2,'xlabel'),'string','Pro CV2')
            set(get(hAx2,'ylabel'),'string','Anti CV2')
            export_fig(hF2,'-pdf','-append',allNeuronSumFigName)
            
            delete(hF2)
            %TODO for each time period tested also produce a tuned spider
        end
        
        
          hFb = figure;
        hAb = axes('parent',hFb);
         sigPercents2 = sigPercents([2 3 1],:);
        [hPatches] = createGroupedBarGraph(sigPercents2,zeros(3,2),hAb,ccBar,0.3);
        set(hAb,'xtick',[1:3],'xticklabel',{'typeInst','dirInst','sacc'})
        ylim([0 100])
        ylabel('Percent Modulated')
        title(['CV2 Different PvA, Unit ' num2str(unitNum)])
        export_fig(hFb,'-pdf','-append',allNeuronSumFigName)
        delete(hFb)
        
        
        
        %compare lateral and vermis baseline rates and cv2s
        %now do 1 final plot that is a line with errorbars showing the mean and
        %sem bseline firing rate for each neuron and then also divide into
        %vermis and lateral
        hF = figure;
       % baseRates = {baselineRates{recordingsToAnalyse,unitNum}};
        baseRates = {baselineRates{:,unitNum}};
        baseRates(cellfun('isempty',baseRates))={nan};
        
        muBaseRate = cellfun(@(x) mean(x),baseRates);
        semBaseRate = cellfun(@(x) std(x)/sqrt(length(x)),baseRates);
        
        hA = axes('parent',hF);
        errorbar(vermNeurons,muBaseRate(vermNeurons),semBaseRate(vermNeurons),'color','b')
        hold on
        errorbar(latNeurons,muBaseRate(latNeurons),semBaseRate(latNeurons),'color','r')
        %now also plot line at overall mean
        line([0 200],[nanmean(muBaseRate(:)) nanmean(muBaseRate(:))],'color','k','xliminclude','off')
        line([0 200],[nanmean(muBaseRate(vermNeurons),2) nanmean(muBaseRate(vermNeurons),2)],'color','b','xliminclude','off')
        line([0 200],[nanmean(muBaseRate(latNeurons),2) nanmean(muBaseRate(latNeurons),2)],'color','r','xliminclude','off')
        %TODO also show mean of lat and verm sepearately and split the errorbar
        %call int he same way
        set(hA,'xtick',[1:length(recordingsToAnalyse)],...
            'xticklabel',{neuronList(recordingsToAnalyse).neuronName});
        xticklabel_rotate
        hL = createLegend([0 0 1;1 0 0],{'Vermis';'Lateral'},{'-','-'},hA,[0.7 0.7 0.2 0.2]);
        suptitle(['Baseline Firing Rate, Unit ' num2str(unitNum)])
        ylabel('Firing Rate (Hz)')
        
        export_fig(hF,'-pdf','-append',allNeuronSumFigName)
        delete(hF)
        
        
        %do exactly the same as above but for teh CV2 values
        hF = figure;
       % baseCV2s = {baselineCV2s{recordingsToAnalyse,unitNum}};
        baseCV2s = {baselineCV2s{:,unitNum}};
        baseCV2s(cellfun('isempty',baseCV2s))={nan};
        
        

muBaseCV = cellfun(@(x) nanmean(x,2),baseCV2s);
        semBaseCV = cellfun(@(x) nanstd(x,1,2)/sqrt(length(x)),baseCV2s);
        hA = axes('parent',hF);
        errorbar(vermNeurons,muBaseCV(vermNeurons),semBaseCV(vermNeurons),'color','b')
        hold on
        errorbar(latNeurons,muBaseCV(latNeurons),semBaseCV(latNeurons),'color','r')
        
        %now also plot line at overall mean
        line([0 200],[mean(muBaseCV) mean(muBaseCV)],'color','k','xliminclude','off')
        line([0 200],[mean(muBaseCV(vermNeurons)) mean(muBaseCV(vermNeurons))],'color','b','xliminclude','off')
        line([0 200],[mean(muBaseCV(latNeurons)) mean(muBaseCV(latNeurons))],'color','r','xliminclude','off')
        
        
        set(hA,'xtick',[1:length(recordingsToAnalyse)],...
            'xticklabel',{neuronList(recordingsToAnalyse).neuronName});
        xticklabel_rotate
        hL = createLegend([0 0 1;1 0 0],{'Vermis';'Lateral'},{'-','-'},hA,[0.7 0.7 0.2 0.2]);
        suptitle(['Baseline CV2s, Unit ' num2str(unitNum)])
        ylabel('CV2')
        
        export_fig(hF,'-pdf','-append',allNeuronSumFigName)
        delete(hF)
        
        
        %now do scatter graph comparing rate and cv2 with different colours
        %for lat and verm neurons
        scatCc =  [0 0 1;1 0 0];
        hF = figure;
        hA = axes('parent',hF);
        scatLabels = nan(max(recordingsToAnalyse),1);
        scatLabels(vermNeurons) = 1;
        scatLabels(latNeurons) = 2;
        scatLabels(isnan(scatLabels))=1; %this will make it a vermis neurons but as the dat is nan it won't acutally be displayed
        scatter(muBaseCV,muBaseRate,30,scatCc(scatLabels,:),'filled','parent',hA)
        xlabel('CV2')
        ylabel('Firing Rate (Hz)')
        export_fig(hF,'-pdf','-append','-opengl',allNeuronSumFigName)
        delete(hF)
    end
    
    save(allKuniResFile,'allNeuronResults','recList')
    
else
    
    load(allKuniResFile)
end

a=1;
%% needs adding to above cell

%     %now do a rOC rurve at each time point
%     tPs = {'saccRates','typeInstRates','dirInstRates'};
%     for tP = 1:3
%         thisTimePoint = tPs{tP};
%         thisProStats =  kuniResProStats.(thisTimePoint);
%         thisAntiStats =  kuniResAntiStats.(thisTimePoint);
%         
%         thisCombinedStats =  [[thisProStats ;ones(1,length(thisProStats))]  [thisAntiStats ;2*ones(1,length(thisAntiStats))]];
%         [X,Y,T,AUC,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(thisCombinedStats(2,:),thisCombinedStats(1,:),1);
%         figure; plot(X,Y)
%         line([0 1],[0 1],'color','k','linestyle','-.')
%         title([neuronList(neuronNum).neuronName ' ROC ' thisTimePoint])
%         export_fig(gcf,'-pdf','-append',thisNeuronSumFigName)
%     end
%% POPULATION BASED MI.

%here we will take the trials from all neurons and combine them together
%and then sort by saccade amplitude. Then we cut the trials into bins of
%0.5degree amplitude (mixed pro and anti) and for each bin calculate the
%mean of the pro and anti trials seperately. So we ahve a results matrix
%that is (numBins,2). Now we do the MI on that matrix to see if at the
%population level a DCN neuron using the mean rate output could
%discriminate between pro and anti saccades

if doPopMi
%TODO also should caompare between lateral and vermis and also need to
%think about what time window to use. For now we should stick with the
%saccade window
numNeurons = max(recordingsToAnalyse);
tempProTrials = cell(numNeurons,1);
tempAntiTrials = cell(numNeurons,1);
for neuronNum = recordingsToAnalyse
    %remove the lfpEoch filed as not present in all analysed
    %recordings
    if ~isempty(popRes(neuronNum,unitNum).alignedTrials)
        tempRes = popRes(neuronNum,unitNum).alignedTrials.Pro;
        if isfield(tempRes,'lfpEpoch')
            tempProTrials{neuronNum} = rmfield(tempRes,'lfpEpoch');
            
        else
            tempProTrials{neuronNum} =tempRes;
        end
        tempRes = popRes(neuronNum,unitNum).alignedTrials.Anti;
        if isfield(tempRes,'lfpEpoch')
            tempAntiTrials{neuronNum} = rmfield(tempRes,'lfpEpoch');
        else
            tempAntiTrials{neuronNum} = tempRes;
        end
        
    end
    
end


binEdges = 0:0.5:20;
   allProTrials = [tempProTrials{:}];
        allAntiTrials = [tempAntiTrials{:}];
        
        trialStructure = allProTrials;
        
        
        
        [binnedResults.Pro binnnedVals.Pro binnedCounts.Pro] = binByStatAndGetRates(allProTrials,sortByString,binEdges)
        [binnedResults.Anti binnnedVals.Anti binnedCounts.Anti] = binByStatAndGetRates(allAntiTrials,sortByString,binEdges)

windowInds = 1975:2200;
%TODO loop over different bin placements to check if particualr period is sensitive to difference
%TODO do for lat and verm seperately


        trialTypes = {'Pro','Anti'};
        numBins = length(binEdges)-1;
        forMi = nan(numBins*2,2); %first column i s value of statistic the second is th label (1 or 2)
        for trialTypeNum = 1:2
            thisTrialType = trialTypes{trialTypeNum}
            theseResults = binnedResults.(char(thisTrialType));
            startBinNum = (trialTypeNum-1)*numBins+1; %to place correctly in results matrix, 
           
            forMi(startBinNum:startBinNum+numBins-1,1) = cellfun(@(x) nanmean(x(windowInds(1):windowInds(2))),theseResults);
            forMi(startBinNum:startBinNum+numBins-1,2) = trialTypeNum*ones(1,numBins);
        end
        
        %now we need to remove any bins which have no trials in either
        %condition
        
        badBins = find([binnedCounts.Anti]==0 | [binnedCounts.Pro]==0);
        forMi(badBins+numBins,:) = []; %remove from anti trials first
         forMi(badBins,:) = []; %now pro bins
        
         
         [MI confLims sigLevel] = miAndConfLims(forMi(:,2),forMi(:,1),2000);
         
         
        %can put the slection of lat and eerm neurons here
        
                %now we need to decide a time period (trials here are aligned to
        %saccade=0 so lets use our noral saccade window of -0.025:0.2.
        %%TODO later we will use the sliding window approach
        %previous results extracted with a +-2s win so check new ones are
        %different and also the line above where 4001 is hardcoded needs to
        %change
        
%        rateMat = cell2mat({trialStructure.sdf}'); %TODO here could be swith for CV2
%         [val srtOrder] = sort([trialStructure.(sortByString)]); 
%         rateMat = rateMat(srtOrder,:);
%         
%       
%         
%         numBins = length(binEdges)-1;
%         
%         binnedResults = cell(numBins,1); %for keeping average whole traces
%         %binnedStats = nan(numBins,2);  %first colum mean, second is sd
%         binnnedVals = nan(numBins,1);
%         binnedCounts = zeros(numBins,1);
%         for binNum = 1:numBins
%             
%             theseTrials = find(val>=binEdges(binNum) & val<binEdges(binNum+1))
%             if ~isempty(theseTrials)
%                 
%                 if length(theseTrials)>1
%                     binnedResults{binNum} = mean(rateMat(theseTrials,:));
%                     binnnedVals(binNum) = mean(val(theseTrials));
%                     binnedCounts(binNum) = length(theseTrials);
%                 else
%                     binnedResults{binNum} = rateMat(theseTrials,:);
%                     binnnedVals(binNum) = val(theseTrials);
%                     binnedCounts(binNum) = 1;
%                 end
%                 
%                 
%             else
%                 binnedResults{binNum} = nan(1,4001);
%                 binnnedVals(binNum) = NaN;
%             end
%         end
%         

end      
    %% MI
%ide is to test for each neuron the mutual information about the stimulus
%(1v2 or 1:8 or 1:16?)
numMontCarlIter = 2000;
unitNum = 1;


if doMI
miVp = nan(max(recordingsToAnalyse),3,3); %[mi sacc, p sacc,mi p sacc],timePeriod
for neuronNum = recordingsToAnalyse;
    
    proSaccRates =  [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.saccRates];
    antiSaccRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.saccRates];
    proTypeInstRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.typeInstRates];
    antiTypeInstRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.typeInstRates];
    
    totalNumTrials = length(proSaccRates)+length(antiSaccRates);
    forMi = nan(totalNumTrials,3); %first column is trial type (1=pro,2=anti), second is saccade rate, 3rd is tpyeInstRate
    
    forMi(1:length(proSaccRates),1) = 1;
    forMi(1:length(proSaccRates),2) = proSaccRates;
    forMi(1:length(proSaccRates),3) = proTypeInstRates;
    forMi(length(proSaccRates)+1:end,1) = 2;
    forMi(length(proSaccRates)+1:end,2) = antiSaccRates;
    forMi(length(proSaccRates)+1:end,3) = antiTypeInstRates;
    
    
    hF = figure;
    hAx2 = subplot(1,3,2,'parent',hF,'colororder',[0 1 0;1 0 0]);
    hAx1 = subplot(1,3,1,'parent',hF,'colororder',[0 1 0;1 0 0],'nextplot','replacechildren');
    hAx3 = subplot(1,3,3,'parent',hF,'colororder',[0 1 0;1 0 0]);
    axes(hAx2)
    gscatter(forMi(:,2), forMi(:,3), forMi(:,1),'gr');
    xlabel('saccade period rate (hz)')
    ylabel('type inst period rate (hz)')
    
    matchingNum = min(length(proSaccRates),length(antiSaccRates)); %maximum number of trials in both
    hist([proSaccRates(1:matchingNum);antiSaccRates(1:matchingNum)]',100,'parent',hAx1);
    
    axes(hAx1)
    legend('1','2')
    %[hSacc pSacc] = ttest(proSaccRates(1:matchingNum),antiSaccRates(1:matchingNum));
    [pSacc hSacc] = signrank(proSaccRates(1:matchingNum),antiSaccRates(1:matchingNum));
    %saccMI = mutualinfo(forMi(:,1),forMi(:,2));
    testMI = MutualInformation(round(forMi(:,2)),forMi(:,1));
    
    set(get(hAx1,'xlabel'),'string','saccade period rate (hz)')
    
    %now using monte carlo reshuffling 2000 times recalculate the
    %MI each time and get confidence limits.
    [saccMI confLims saccSigLevel] = miAndConfLims(forMi(:,1),forMi(:,2),numMontCarlIter);
    
%     numLabels = length(forMi(:,1));
%     montCarlRes = nan(numMontCarlIter,1);
%     %tic
%     for MontCarlIterNum = 1:numMontCarlIter
%         shuffledLabels = forMi(randperm(numLabels),1);
%         montCarlRes(MontCarlIterNum) = mutualinfo(forMi(:,2),shuffledLabels);
%     end
%     %toc
%     %confLims = prctile(montCarlRes,[2.5 97.5]);
%     [val ord] = sort([montCarlRes; saccMI]);
%     posInList = find(ord==numMontCarlIter+1);
%    sigLevel = posInList/(numMontCarlIter+1);
    miVp(neuronNum,1,1) = saccMI;
    miVp(neuronNum,2,1) = pSacc;
    miVp(neuronNum,3,1) = saccSigLevel;
    
    set(get(hAx1,'title'),'string',['MI = ' num2str(saccMI) ', p = ' num2str(pSacc) ',MIp = ' num2str(miVp(neuronNum,3))])
    
    hist([proTypeInstRates(1:matchingNum);antiTypeInstRates(1:matchingNum)]',100,'parent',hAx3);
    
    %[hType pType] = ttest(proTypeInstRates(1:matchingNum),antiTypeInstRates(1:matchingNum));
    [pType hType] = signrank(proTypeInstRates(1:matchingNum),antiTypeInstRates(1:matchingNum));
    typeInstMI = mutualinfo(forMi(:,1),forMi(:,3));
    
    
    %now using monte carlo reshuffling 2000 times recalculate the
    %MI each time and get confidence limits.
    [typeInstMI confLims typeInstSigLevel] = miAndConfLims(forMi(:,1),forMi(:,3),numMontCarlIter);
    
%     numLabels = length(forMi(:,1));
%     montCarlRes = nan(numMontCarlIter,1);
%     %tic
%     for MontCarlIterNum = 1:numMontCarlIter
%         shuffledLabels = forMi(randperm(numLabels),1);
%         montCarlRes(MontCarlIterNum) = mutualinfo(forMi(:,3),shuffledLabels);
%     end
%     %toc
%     %confLims = prctile(montCarlRes,[2.5 97.5]);
%     [val ord] = sort([montCarlRes; saccMI]);
%     posInList = find(ord==numMontCarlIter+1);
%     sigLevel = posInList/(numMontCarlIter+1);
    
    miVp(neuronNum,1,3) = typeInstMI;
    miVp(neuronNum,2,3) = pType;
    miVp(neuronNum,3,3) = typeInstSigLevel;
    
    set(get(hAx3,'title'),'string',['MI = ' num2str(typeInstMI) ',p = ' num2str(pType) ',MIp = ' num2str(miVp(neuronNum,6))])
    set(get(hAx3,'xlabel'),'string','type inst period rate (hz)')
    suptitle(neuronList(neuronNum).neuronName)
    try
        
        export_fig(hF,'-pdf','-append',miSumFigName)
    catch
        export_fig(hF,'-pdf',miSumFigName)
    end
    delete(hF)
end
hF2 = figure;
hAx4 = subplot(1,2,1,'parent',hF2);
hAx5 = subplot(1,2,2,'parent',hF2);
scatter(miVp(:,3),miVp(:,2),'parent',hAx4)
set(get(hAx4,'xlabel'),'string','saccade MIp')
set(get(hAx4,'ylabel'),'string','saccade p')


scatter(miVp(:,6),miVp(:,5),'parent',hAx5)
set(get(hAx5,'xlabel'),'string','type inst MIp')
set(get(hAx5,'ylabel'),'string','type inst p')

%this time use the MI significant neurons as markers
testPeriods = {'Saccade Period','Direction Instruction','Type Instruction'}; %
%now loop over the 3 testing periods and porduce a figure for each
for testPeriodNum = [1 3]
    hFsacc = figure;
    hAsacc = axes('parent',hFsacc);
    hF2 = figure; %for lat and verm subplots
    
    
    muRates = allNeuronResults(max(recordingsToAnalyse),unitNum).muRates(:,:,testPeriodNum);
    emptyRecs = find(cellfun('isempty',muRates));
    muRates(emptyRecs) = {nan};
    
    proSaccadeMeanRates = [muRates{:,1}];
    antiSaccadeMeanRates = [muRates{:,2}];
    %proAntiSig = [muRates{:,3,testPeriodNum}];
  
        proAntiSig = miVp(:,3,testPeriodNum); %3rd colum for sacc MI significance, 6th for type inst
        
  
    proAntiSig = proAntiSig<0.05; %find those significant at given level
    
    set(hAsacc,'xlim',axLims(unitNum,:),'ylim',axLims(unitNum,:),'Nextplot','add')
    line(axLims(unitNum,:),axLims(unitNum,:),'color','k','parent',hAsacc); %LOI
    scatter(proSaccadeMeanRates,antiSaccadeMeanRates,20,'k','parent',hAsacc);
    %now scatter signifcant points with filled circles
    sigSacs = find(proAntiSig);
    scatter(proSaccadeMeanRates(sigSacs),antiSaccadeMeanRates(sigSacs),20,'k','filled','parent',hAsacc)
    
    set(get(hAsacc,'title'),'string',['MI All Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
    set(get(hAsacc,'xlabel'),'string','Pro Rate (Hz)')
    set(get(hAsacc,'ylabel'),'string','Anti Rate (Hz)')
    export_fig(hFsacc,'-pdf','-append',miSumFigName)
    
    
    %now split into vermis and lateral cerebellar recordings
    sigVermSaccs = intersect(sigSacs,vermNeurons);
    sigLatSaccs = intersect(sigSacs,latNeurons);
    
    hAx1 = subplot(1,2,1,'parent',hF2);
    set(hAx1,'xlim',axLims(unitNum,:),'ylim',axLims(unitNum,:),'Nextplot','add')
    line(axLims(unitNum,:),axLims(unitNum,:),'color','k','parent',hAx1); %LOI
    scatter(proSaccadeMeanRates(vermNeurons),antiSaccadeMeanRates(vermNeurons),20,'k','parent',hAx1);
    %now scatter signifcant points with filled circles
    %sigSacs = find(proAntiSig);
    scatter(proSaccadeMeanRates(sigVermSaccs),antiSaccadeMeanRates(sigVermSaccs),20,'k','filled','parent',hAx1);
    
    set(get(hAx1,'title'),'string',['MI Vermis Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
    set(get(hAx1,'xlabel'),'string','Pro Rate (Hz)')
    set(get(hAx1,'ylabel'),'string','Anti Rate (Hz)')
    %export_fig(hF2,'-pdf','-append',miSumFigName)
    
    
    hAx2 = subplot(1,2,2,'parent',hF2);
    set(hAx2,'xlim',axLims(unitNum,:),'ylim',axLims(unitNum,:),'Nextplot','add')
    line(axLims(unitNum,:),axLims(unitNum,:),'color','k','parent',hAx2); %LOI
    scatter(proSaccadeMeanRates(latNeurons),antiSaccadeMeanRates(latNeurons),20,'k','parent',hAx2);
    %now scatter signifcant points with filled circles
    %sigSacs = find(proAntiSig);
    scatter(proSaccadeMeanRates(sigLatSaccs),antiSaccadeMeanRates(sigLatSaccs),20,'k','filled','parent',hAx2)
    
    set(get(hAx2,'title'),'string',['MI Lateral Neurons ' testPeriods{testPeriodNum} ' unit ' num2str(unitNum)])
    set(get(hAx2,'xlabel'),'string','Pro Rate (Hz)')
    set(get(hAx2,'ylabel'),'string','Anti Rate (Hz)')
    export_fig(hF2,'-pdf','-append',miSumFigName)
    
    delete(hF2)
    %TODO for each time period tested also produce a tuned spider
end

end
%% SLIDING WINDOW MI INSTRUCTION


measureName = 'rate';
unitNum = 1;
numSlidingWindowSteps = 2000; %number of times to advance the window and test again.
%slidingWindowSize = 0.1; %both in ms
slidingWindowStepSize = 0.001;
slidingWindowStartPoint = -1; %in s [before after] saccade time
if doSlidWinMIinst
    %loop over different sized windows and check which give clearest results.
    windowSizesToTest = 0.1; %[0.1 0.2];
    windowSizeCounter = 1;
    %stInit = cell(length(windowSizesToTest),1);
    %allWindowResults = struct('results',stInit,'windowSize',stInit);
    
    
    for slidingWindowSize = windowSizesToTest;
        miSlidSumFigName = [resultsDir filesep '2016InstructslidingMI-' num2str(slidingWindowSize*1000) 'msWin' measureName '3.pdf'];
        
        
        thisWindowSizeResults = cell(max(recordingsToAnalyse),3);%first colum is results second is mis, 3rd is window size 
        allNeuronSignificantMI = nan(max(recordingsToAnalyse),numSlidingWindowSteps); %to keep track of which bits are above confidence
        for neuronNum = recordingsToAnalyse
            disp(['Sliding Window MI neuron '  num2str(neuronNum)])
            thisNeuronUnits = recList{2,neuronNum}(1,3);
            thisNeuronUnits = [thisNeuronUnits{:}];
            thisNeuronUnitNum = thisNeuronUnits(unitNum)
            %     a=1;
            %       proSaccRates =  [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.saccRates];
            %     antiSaccRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.saccRates];
            %     proTypeInstRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.typeInstRates];
            %     antiTypeInstRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.typeInstRates];
            %
            trialsToAnalyse = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
            numTrials  = length(trialsToAnalyse);
            theseTrials = wholeNeuronResults(neuronNum).allStableTrials;
            %TODO is currently aligned to fixation
            %intiliases results matrix
            thisNeuronResults = nan(numSlidingWindowSteps,numTrials);
            thisNeuronMIs = nan(numSlidingWindowSteps,4); %first is MI, second lowed conf, 3rd upper conf, 4th p
            
            
            
            thisNeuronLabels = nan(2,numTrials); %second row is full condtion code first is pro=1 anti=2
            thisNeuronLabels(2,:) = [theseTrials(trialsToAnalyse).conditionCode];
            thisNeuronLabels(1,:) = ismember(thisNeuronLabels(2,:),trialConditionCodes(:,1));
            thisNeuronLabels(1,thisNeuronLabels(1,:)==0) = 2;
            
            %here every direction should be labelled by its direction code
            
            
            
              % slMi = cell(numTrials,1);
                
                
              %go through each trail and get all spikes that fall within given window
              for trialNum = 1:numTrials %TODO this loop order is slower as recalcs CV2 for each trial each time
                   thisTrial = theseTrials(trialsToAnalyse(trialNum));
                  thisTrialCV2 = [calculateCV2(thisTrial.alignedSpikes{thisNeuronUnitNum}); nan]; %calcualte CV2 and add nan for last spike
                  
                  if any(thisTrialCV2)
                  for slidingStepNum = 1:numSlidingWindowSteps
                      
                      thisMiWindowStart = slidingWindowStartPoint+slidingWindowStepSize*(slidingStepNum-1);
                      thisMiWindow = [thisMiWindowStart thisMiWindowStart+slidingWindowSize]; %TODO update each iter
                      
                      cv2Med= nan;
                      
                     
                      miPeriodSpikesInds = cellfun(@(x) find(x>(thisTrial.goCueTime+thisMiWindow(1)) & x<(thisTrial.goCueTime+thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
                      
                      if strcmpi(measureName,'CV2')
                          
                          if length([miPeriodSpikesInds{:}])>=2
                              
                              cv2Med = nanmedian(thisTrialCV2([miPeriodSpikesInds{:}]));
                          end
                          
                          
                          thisNeuronResults(slidingStepNum,trialNum) = cv2Med;
                          
                      elseif strcmpi(measureName,'rate')
                          
                          if ~isempty(miPeriodSpikesInds)
                              thisNeuronResults(slidingStepNum,trialNum) = length([miPeriodSpikesInds{:}]);
                              
                          end
                      end
                  end
                  end
              end      
 
                
                for slidingStepNum = 1:numSlidingWindowSteps
                    %TODO need to remove any nan from both labels and results
                    %after building the matrix of trials then do MI caclulation
                    %for all timeWindows
                    goodTrials = find(~isnan(thisNeuronResults(slidingStepNum,:)));
                    
                    if ~isempty(goodTrials)
                        [thisMI thisConfLims thisSigLevel] = miAndConfLims(thisNeuronResults(slidingStepNum,goodTrials),thisNeuronLabels(1,goodTrials),numMontCarlIter);
                    else
                        thisMI = nan;
                        thisConfLims = [nan nan];
                        thisSigLevel = nan;
                    end
                    thisNeuronMIs(slidingStepNum,1) = thisMI;
                    thisNeuronMIs(slidingStepNum,2) = thisConfLims(1);
                    thisNeuronMIs(slidingStepNum,3) = thisConfLims(2);
                    thisNeuronMIs(slidingStepNum,4) = thisSigLevel;
                    
                end
            thisWindowSizeResults{neuronNum,1} = thisNeuronResults;
            thisWindowSizeResults{neuronNum,2} = thisNeuronMIs;
           
                        %need to get the times that MI is significant
            %TODO smooth first
            
            thisNeuronSigMI = thisNeuronMIs(:,1);
            thisNeuronSigMI(thisNeuronMIs(:,1)<thisNeuronMIs(:,3)) = nan;
            %now keep the significant bits of the MI trace and NaN the rest.
            allNeuronSignificantMI(neuronNum,:) = thisNeuronSigMI;
           
            %generate MI time vec and compensate for half window size so
            %that the results of each window is plotted in its midpoint
            %timewise.
            endPoint = slidingWindowStartPoint+slidingWindowStepSize*numSlidingWindowSteps;
            miTimeVec = linspace(slidingWindowStartPoint,endPoint,numSlidingWindowSteps);
            
            miTimeVec = miTimeVec+slidingWindowSize/2;
             thisWindowSizeResults{neuronNum,3} = miTimeVec;
             %now plot rate results
            hF = figure;
            hAx1 = subplot(2,1,1,'parent',hF);
            plot(miTimeVec,thisNeuronResults,'parent',hAx1)
            hold on;
            plot(miTimeVec,nanmean(thisNeuronResults,2),'k','parent',hAx1) %TODO should plot mean for pro and anti seperately
            plot(miTimeVec,nanmean(thisNeuronResults(:,thisNeuronLabels(1,:)==1),2),'g','parent',hAx1)
            plot(miTimeVec,nanmean(thisNeuronResults(:,thisNeuronLabels(1,:)==2),2),'r','parent',hAx1)
            
            hAx2 = subplot(2,1,2,'parent',hF)
            plot(miTimeVec,thisNeuronMIs(:,1),'k')
            hold on
            plot(miTimeVec,thisNeuronMIs(:,2),'k-.')
            plot(miTimeVec,thisNeuronMIs(:,3),'k-.')
            
            suptitle([neuronList(neuronNum).neuronName 'MI of ' measureName ' over time, aligned to go cue, ' num2str(1000*slidingWindowSize) 'ms window'])
            export_fig(hF,'-pdf','-append',miSlidSumFigName)
            delete(hF)
            
%             %now do exactly the same but when aligned to trialStart which is at
%             %fixation
%             for slidingStepNum = 1:numSlidingWindowSteps
%                 
%                 thisMiWindowStart = slidingWindowStartPoint+slidingWindowStepSize*(slidingStepNum-1);
%                 thisMiWindow = [thisMiWindowStart thisMiWindowStart+slidingWindowSize]; %TODO update each iter
%                 % slMi = cell(numTrials,1);
%                 
%                 
%                 %go through each trail and get all spikes that fall within given window
%                 for trialNum = 1:numTrials
%                     
%                     
%                     
%                     thisTrial = theseTrials(trialsToAnalyse(trialNum));
%                     cv2Med= nan;
%                     
%                     thisTrial = theseTrials(trialsToAnalyse(trialNum));
%                     miPeriodSpikesInds = cellfun(@(x) find(x>(thisTrial.saccadeTime+thisMiWindow(1)) & x<(thisTrial.saccadeTime+thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
%                     
%                     if strcmpi(measureName,'CV2')
%                         thisTrialCV2 = calculateCV2(thisTrial.alignedSpikes{thisNeuronUnitNum});
%                         if size(miPeriodSpikesInds,1)>=3
%                         
%                             cv2Med = nanmedian(thisTrialCV2([miPeriodSpikesInds{:}]));
%                         end
%                         
%                         
%                         thisNeuronResults(slidingStepNum,trialNum) = cv2Med;
%                         
%                     elseif strcmpi(measureName,'rate')
%                         
%                         if ~isempty(miPeriodSpikesInds)
%                          thisNeuronResults(slidingStepNum,trialNum) = length(miPeriodSpikesInds);
%                          
%                         end
%                     end
%                     
% %                     miPeriodSpikes = cellfun(@(x) x(x>(thisMiWindow(1)) & x<(thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
% %                     
% %                     miSpikeRate = length([miPeriodSpikes{:}])/slidingWindowSize;
% %                    % thisNeuronResults(slidingStepNum,trialNum) = miSpikeRate;
% %                      miPeriodSpikes = [miPeriodSpikes{:}];
% %                     
% %                     if size(miPeriodSpikes,1)>=3
% %                         theseISIs = diff(miPeriodSpikes);
% %                         
% %                         %fCV2 = @(x,y) 2*abs(x-y)/(x+y);
% %                         
% %                         isiMat = [circshift(theseISIs,1) theseISIs];
% %                         
% %                         
% %                         cv2Vec = arrayfun(@(x,y) fCV2(x,y),isiMat(:,1),isiMat(:,2));
% %                         cv2Vec(1) = nan; %first spike cv2 is not real
% %                         cv2Med = nanmedian(cv2Vec);
% %                         
% %                         thisNeuronResults(slidingStepNum,trialNum) = cv2Med;
% %                     end
%                     
%                     
%                     %after building the matrix of trials then do MI caclulation for this timeWindow
%                     
%                 end
%                 if ~isempty(goodTrials)
%                   goodTrials = find(~isnan(thisNeuronResults(slidingStepNum,:)));
%                    [thisMI thisConfLims thisSigLevel] = miAndConfLims(thisNeuronResults(slidingStepNum,goodTrials),thisNeuronLabels(1,goodTrials),numMontCarlIter);
%                
%                else
%                 thisMI = nan;
%                 thisConfLims = [nan nan];
%                 thisSigLevel = nan;
%                 end
%                 
%                 thisNeuronMIs(slidingStepNum,1) = thisMI;
%                 thisNeuronMIs(slidingStepNum,2) = thisConfLims(1);
%                 thisNeuronMIs(slidingStepNum,3) = thisConfLims(2);
%                 thisNeuronMIs(slidingStepNum,4) = thisSigLevel;
%             end
%             
%             
%             %now plot rate results
%             hF = figure;
%             hAx1 = subplot(2,1,1,'parent',hF);
%             plot(thisNeuronResults,'parent',hAx1)
%             hold on;
%             plot(nanmean(thisNeuronResults,2),'k','parent',hAx1) %TODO should plot mean for pro and anti seperately
%             plot(nanmean(thisNeuronResults(:,thisNeuronLabels(1,:)==1),2),'g','parent',hAx1)
%             plot(nanmean(thisNeuronResults(:,thisNeuronLabels(1,:)==2),2),'r','parent',hAx1)
%             
%             hAx2 = subplot(2,1,2,'parent',hF)
%             plot(thisNeuronMIs(:,1),'k')
%             hold on
%             plot(thisNeuronMIs(:,2),'k-.')
%             plot(thisNeuronMIs(:,3),'k-.')
%             
%             suptitle([neuronList(neuronNum).neuronName 'MI of ' measureName ' over time, aligned to fixation, ' num2str(1000*slidingWindowSize) 'ms window'])
%             export_fig(hF,'-pdf','-append',miSlidSumFigName)
%             delete(hF)

            
        end
        save([resultsDir filesep 'InstructslidingWidnowMI_' num2str(slidingWindowSize) 'msWin_' measureName '.mat'],'thisWindowSizeResults','slidingWindowSize')
        %save to big results struct
        %allWindowResults(windowSizeCounter).results = thisWindowSizeResults;
        %allWindowResults(windowSizeCounter).windowSize = slidingWindowSize;
        %windowSizeCounter = windowSizeCounter+1;
        %now make plot showing which bits were significant for both types
        %of neuron
        hF2 = figure;
        hAx1 = subplot(2,1,1,'parent',hF2)
        plot(miTimeVec,allNeuronSignificantMI(latNeurons,:),'color','r','parent',hAx1)
        hold on
        plot(miTimeVec,allNeuronSignificantMI(vermNeurons,:),'color','b','parent',hAx1)
        hAx2 = subplot(2,1,2,'parent',hF2)
        plot(miTimeVec,sum(~isnan(allNeuronSignificantMI(vermNeurons,:)))./length(vermNeurons),'color','b','linewidth',2,'parent',hAx2)
        hold on
        plot(miTimeVec,sum(~isnan(allNeuronSignificantMI(latNeurons,:)))./length(latNeurons),'color','r','linewidth',2,'parent',hAx2)
        legend('verm','lateral')
        suptitle(['All significant time periods pro v anti ' measureName num2str(1000*slidingWindowSize) 'ms window, go cue aligned'])
        export_fig(hF2,'-pdf','-append',miSlidSumFigName)
        delete(hF2)
        
        
    end
end
%% SLIDING WINDOW MI


measureName = 'rate';
unitNum = 1;
numSlidingWindowSteps = 2000; %number of times to advance the window and test again.
%slidingWindowSize = 0.1; %both in ms
slidingWindowStepSize = 0.001;
slidingWindowStartPoint = -1; %in s [before after] saccade time
if doSlidWinMI
    %loop over different sized windows and check which give clearest results.
    windowSizesToTest = 0.1; %[0.1 0.2];
    windowSizeCounter = 1;
    %stInit = cell(length(windowSizesToTest),1);
    %allWindowResults = struct('results',stInit,'windowSize',stInit);
    
    
    for slidingWindowSize = windowSizesToTest;
        miSlidSumFigName = [resultsDir filesep '2016slidingMI-' num2str(slidingWindowSize*1000) 'msWin' measureName '.pdf'];
        
        
        thisWindowSizeResults = cell(max(recordingsToAnalyse),3);%first colum is results second is mis, 3rd is window size 
        allNeuronSignificantMI = nan(max(recordingsToAnalyse),numSlidingWindowSteps); %to keep track of which bits are above confidence
        for neuronNum = recordingsToAnalyse
            disp(['Sliding Window MI neuron '  num2str(neuronNum)])
            thisNeuronUnits = recList{2,neuronNum}(1,3);
            thisNeuronUnits = [thisNeuronUnits{:}];
            thisNeuronUnitNum = thisNeuronUnits(unitNum)
            %     a=1;
            %       proSaccRates =  [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.saccRates];
            %     antiSaccRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.saccRates];
            %     proTypeInstRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.typeInstRates];
            %     antiTypeInstRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.typeInstRates];
            %
            trialsToAnalyse = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
            numTrials  = length(trialsToAnalyse);
            theseTrials = wholeNeuronResults(neuronNum).allStableTrials;
            
            %intiliases results matrix
            thisNeuronResults = nan(numSlidingWindowSteps,numTrials);
            thisNeuronMIs = nan(numSlidingWindowSteps,4); %first is MI, second lowed conf, 3rd upper conf, 4th p
            
            
            
            thisNeuronLabels = nan(2,numTrials); %second row is full condtion code first is pro=1 anti=2
            thisNeuronLabels(2,:) = [theseTrials(trialsToAnalyse).conditionCode];
            thisNeuronLabels(1,:) = ismember(thisNeuronLabels(2,:),trialConditionCodes(:,1));
            thisNeuronLabels(1,thisNeuronLabels(1,:)==0) = 2;
            
            %here every direction should be labelled by its direction code
            
            
            
              % slMi = cell(numTrials,1);
                
                
              %go through each trail and get all spikes that fall within given window
              for trialNum = 1:numTrials %TODO this loop order is slower as recalcs CV2 for each trial each time
                   thisTrial = theseTrials(trialsToAnalyse(trialNum));
                  thisTrialCV2 = [calculateCV2(thisTrial.alignedSpikes{thisNeuronUnitNum}); nan]; %calcualte CV2 and add nan for last spike
                  
                  if any(thisTrialCV2)
                  for slidingStepNum = 1:numSlidingWindowSteps
                      
                      thisMiWindowStart = slidingWindowStartPoint+slidingWindowStepSize*(slidingStepNum-1);
                      thisMiWindow = [thisMiWindowStart thisMiWindowStart+slidingWindowSize]; %TODO update each iter
                      
                      cv2Med= nan;
                      
                     
                      miPeriodSpikesInds = cellfun(@(x) find(x>(thisTrial.saccadeTime+thisMiWindow(1)) & x<(thisTrial.saccadeTime+thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
                      
                      if strcmpi(measureName,'CV2')
                          
                          if length([miPeriodSpikesInds{:}])>=2
                              
                              cv2Med = nanmedian(thisTrialCV2([miPeriodSpikesInds{:}]));
                          end
                          
                          
                          thisNeuronResults(slidingStepNum,trialNum) = cv2Med;
                          
                      elseif strcmpi(measureName,'rate')
                          
                          if ~isempty(miPeriodSpikesInds)
                              thisNeuronResults(slidingStepNum,trialNum) = length([miPeriodSpikesInds{:}]);
                              
                          end
                      end
                  end
                  end
              end      
 
                
                for slidingStepNum = 1:numSlidingWindowSteps
                    %TODO need to remove any nan from both labels and results
                    %after building the matrix of trials then do MI caclulation
                    %for all timeWindows
                    goodTrials = find(~isnan(thisNeuronResults(slidingStepNum,:)));
                    
                    if ~isempty(goodTrials)
                        [thisMI thisConfLims thisSigLevel] = miAndConfLims(thisNeuronResults(slidingStepNum,goodTrials),thisNeuronLabels(1,goodTrials),numMontCarlIter);
                    else
                        thisMI = nan;
                        thisConfLims = [nan nan];
                        thisSigLevel = nan;
                    end
                    thisNeuronMIs(slidingStepNum,1) = thisMI;
                    thisNeuronMIs(slidingStepNum,2) = thisConfLims(1);
                    thisNeuronMIs(slidingStepNum,3) = thisConfLims(2);
                    thisNeuronMIs(slidingStepNum,4) = thisSigLevel;
                    
                end
            thisWindowSizeResults{neuronNum,1} = thisNeuronResults;
            thisWindowSizeResults{neuronNum,2} = thisNeuronMIs;
           
                        %need to get the times that MI is significant
            %TODO smooth first
            
            thisNeuronSigMI = thisNeuronMIs(:,1);
            thisNeuronSigMI(thisNeuronMIs(:,1)<thisNeuronMIs(:,3)) = nan;
            %now keep the significant bits of the MI trace and NaN the rest.
            allNeuronSignificantMI(neuronNum,:) = thisNeuronSigMI;
           
            %generate MI time vec and compensate for half window size so
            %that the results of each window is plotted in its midpoint
            %timewise.
            endPoint = slidingWindowStartPoint+slidingWindowStepSize*numSlidingWindowSteps;
            miTimeVec = linspace(slidingWindowStartPoint,endPoint,numSlidingWindowSteps);
            
            miTimeVec = miTimeVec+slidingWindowSize/2;
             thisWindowSizeResults{neuronNum,3} = miTimeVec;
             %now plot rate results
            hF = figure;
            hAx1 = subplot(2,1,1,'parent',hF);
            plot(miTimeVec,thisNeuronResults,'parent',hAx1)
            hold on;
            plot(miTimeVec,nanmean(thisNeuronResults,2),'k','parent',hAx1) %TODO should plot mean for pro and anti seperately
            plot(miTimeVec,nanmean(thisNeuronResults(:,thisNeuronLabels(1,:)==1),2),'g','parent',hAx1)
            plot(miTimeVec,nanmean(thisNeuronResults(:,thisNeuronLabels(1,:)==2),2),'r','parent',hAx1)
            
            hAx2 = subplot(2,1,2,'parent',hF)
            plot(miTimeVec,thisNeuronMIs(:,1),'k')
            hold on
            plot(miTimeVec,thisNeuronMIs(:,2),'k-.')
            plot(miTimeVec,thisNeuronMIs(:,3),'k-.')
            
            suptitle([neuronList(neuronNum).neuronName 'MI of ' measureName ' over time, aligned to saccade, ' num2str(1000*slidingWindowSize) 'ms window'])
            export_fig(hF,'-pdf','-append',miSlidSumFigName)
            delete(hF)
            
%             %now do exactly the same but when aligned to trialStart which is at
%             %fixation
%             for slidingStepNum = 1:numSlidingWindowSteps
%                 
%                 thisMiWindowStart = slidingWindowStartPoint+slidingWindowStepSize*(slidingStepNum-1);
%                 thisMiWindow = [thisMiWindowStart thisMiWindowStart+slidingWindowSize]; %TODO update each iter
%                 % slMi = cell(numTrials,1);
%                 
%                 
%                 %go through each trail and get all spikes that fall within given window
%                 for trialNum = 1:numTrials
%                     
%                     
%                     
%                     thisTrial = theseTrials(trialsToAnalyse(trialNum));
%                     cv2Med= nan;
%                     
%                     thisTrial = theseTrials(trialsToAnalyse(trialNum));
%                     miPeriodSpikesInds = cellfun(@(x) find(x>(thisTrial.saccadeTime+thisMiWindow(1)) & x<(thisTrial.saccadeTime+thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
%                     
%                     if strcmpi(measureName,'CV2')
%                         thisTrialCV2 = calculateCV2(thisTrial.alignedSpikes{thisNeuronUnitNum});
%                         if size(miPeriodSpikesInds,1)>=3
%                         
%                             cv2Med = nanmedian(thisTrialCV2([miPeriodSpikesInds{:}]));
%                         end
%                         
%                         
%                         thisNeuronResults(slidingStepNum,trialNum) = cv2Med;
%                         
%                     elseif strcmpi(measureName,'rate')
%                         
%                         if ~isempty(miPeriodSpikesInds)
%                          thisNeuronResults(slidingStepNum,trialNum) = length(miPeriodSpikesInds);
%                          
%                         end
%                     end
%                     
% %                     miPeriodSpikes = cellfun(@(x) x(x>(thisMiWindow(1)) & x<(thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
% %                     
% %                     miSpikeRate = length([miPeriodSpikes{:}])/slidingWindowSize;
% %                    % thisNeuronResults(slidingStepNum,trialNum) = miSpikeRate;
% %                      miPeriodSpikes = [miPeriodSpikes{:}];
% %                     
% %                     if size(miPeriodSpikes,1)>=3
% %                         theseISIs = diff(miPeriodSpikes);
% %                         
% %                         %fCV2 = @(x,y) 2*abs(x-y)/(x+y);
% %                         
% %                         isiMat = [circshift(theseISIs,1) theseISIs];
% %                         
% %                         
% %                         cv2Vec = arrayfun(@(x,y) fCV2(x,y),isiMat(:,1),isiMat(:,2));
% %                         cv2Vec(1) = nan; %first spike cv2 is not real
% %                         cv2Med = nanmedian(cv2Vec);
% %                         
% %                         thisNeuronResults(slidingStepNum,trialNum) = cv2Med;
% %                     end
%                     
%                     
%                     %after building the matrix of trials then do MI caclulation for this timeWindow
%                     
%                 end
%                 if ~isempty(goodTrials)
%                   goodTrials = find(~isnan(thisNeuronResults(slidingStepNum,:)));
%                    [thisMI thisConfLims thisSigLevel] = miAndConfLims(thisNeuronResults(slidingStepNum,goodTrials),thisNeuronLabels(1,goodTrials),numMontCarlIter);
%                
%                else
%                 thisMI = nan;
%                 thisConfLims = [nan nan];
%                 thisSigLevel = nan;
%                 end
%                 
%                 thisNeuronMIs(slidingStepNum,1) = thisMI;
%                 thisNeuronMIs(slidingStepNum,2) = thisConfLims(1);
%                 thisNeuronMIs(slidingStepNum,3) = thisConfLims(2);
%                 thisNeuronMIs(slidingStepNum,4) = thisSigLevel;
%             end
%             
%             
%             %now plot rate results
%             hF = figure;
%             hAx1 = subplot(2,1,1,'parent',hF);
%             plot(thisNeuronResults,'parent',hAx1)
%             hold on;
%             plot(nanmean(thisNeuronResults,2),'k','parent',hAx1) %TODO should plot mean for pro and anti seperately
%             plot(nanmean(thisNeuronResults(:,thisNeuronLabels(1,:)==1),2),'g','parent',hAx1)
%             plot(nanmean(thisNeuronResults(:,thisNeuronLabels(1,:)==2),2),'r','parent',hAx1)
%             
%             hAx2 = subplot(2,1,2,'parent',hF)
%             plot(thisNeuronMIs(:,1),'k')
%             hold on
%             plot(thisNeuronMIs(:,2),'k-.')
%             plot(thisNeuronMIs(:,3),'k-.')
%             
%             suptitle([neuronList(neuronNum).neuronName 'MI of ' measureName ' over time, aligned to fixation, ' num2str(1000*slidingWindowSize) 'ms window'])
%             export_fig(hF,'-pdf','-append',miSlidSumFigName)
%             delete(hF)

            
        end
        save([resultsDir filesep 'slidingWidnowMI_' num2str(slidingWindowSize) 'msWin_' measureName '.mat'],'thisWindowSizeResults','slidingWindowSize')
        %save to big results struct
        %allWindowResults(windowSizeCounter).results = thisWindowSizeResults;
        %allWindowResults(windowSizeCounter).windowSize = slidingWindowSize;
        %windowSizeCounter = windowSizeCounter+1;
        %now make plot showing which bits were significant for both types
        %of neuron
        hF2 = figure;
        hAx1 = subplot(2,1,1,'parent',hF2)
        plot(miTimeVec,allNeuronSignificantMI(latNeurons,:),'color','r','parent',hAx1)
        hold on
        plot(miTimeVec,allNeuronSignificantMI(vermNeurons,:),'color','b','parent',hAx1)
        hAx2 = subplot(2,1,2,'parent',hF2)
        plot(miTimeVec,sum(~isnan(allNeuronSignificantMI(vermNeurons,:)))./length(vermNeurons),'color','b','linewidth',2,'parent',hAx2)
        hold on
        plot(miTimeVec,sum(~isnan(allNeuronSignificantMI(latNeurons,:)))./length(latNeurons),'color','r','linewidth',2,'parent',hAx2)
        legend('verm','lateral')
        suptitle(['All significant time periods pro v anti ' measureName num2str(1000*slidingWindowSize) 'ms window'])
        export_fig(hF2,'-pdf','-append',miSlidSumFigName)
        delete(hF2)
        
        
    end
end
%% SLIDING WINDOW MI BUT BY DIRECTION WITHIN PRO TRIALS

if doSlidWinMI

miSlidSumFigName = [resultsDir filesep '2016slidingMI-' num2str(slidingWindowSize*1000) 'msWinByDir ' measureName '.pdf'];
    
%allWinSizes = cell(5,1);
%winSizeCount = 1;
%loop over different sized windows and check which give clearest results.
for slidingWindowSize = 0.1 %[0.025 0.05 0.075 0.1 0.2];
    
    
    
      thisWindowSizeResults = cell(max(recordingsToAnalyse),3);%first colum is results second is mis, 3rd is window size 

    allNeuronSignificantMI = nan(max(recordingsToAnalyse),numSlidingWindowSteps); %to keep track of which bits are above confidence
    for neuronNum = recordingsToAnalyse
        thisNeuronUnits = recList{2,neuronNum}(1,3);
        thisNeuronUnits = [thisNeuronUnits{:}];
        thisNeuronUnitNum = thisNeuronUnits(unitNum);
        %     a=1;
        %       proSaccRates =  [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.saccRates];
        %     antiSaccRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.saccRates];
        %     proTypeInstRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corProTrials.typeInstRates];
        %     antiTypeInstRates = [allNeuronResults(neuronNum,unitNum).kuniResultsStats.corAntiTrials.typeInstRates];
        %
        trialsToAnalyse = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        numTrials  = length(trialsToAnalyse);
        theseTrials = wholeNeuronResults(neuronNum).allStableTrials;
        
        %intiliases results matrix
        thisNeuronResults = nan(numSlidingWindowSteps,numTrials);
        thisNeuronMIs = nan(numSlidingWindowSteps,4); %first is MI, second lowed conf, 3rd upper conf, 4th p
        
        
        
        thisNeuronLabels = nan(2,numTrials); %second row is full condtion code first is pro=1 anti=2
        thisNeuronLabels(2,:) = [theseTrials(trialsToAnalyse).conditionCode];
        thisNeuronLabels(1,:) = ismember(thisNeuronLabels(2,:),trialConditionCodes(:,1));
        thisNeuronLabels(1,thisNeuronLabels(1,:)==0) = 2;
        
        %now remove all of the anti trials from both the labels and the
        for trialNum = 1:numTrials %TODO this loop order is slower as recalcs CV2 for each trial each time
            thisTrial = theseTrials(trialsToAnalyse(trialNum));
            thisTrialCV2 = [calculateCV2(thisTrial.alignedSpikes{thisNeuronUnitNum}); nan]; %calcualte CV2 and add nan for last spike
            
            if any(thisTrialCV2)
                for slidingStepNum = 1:numSlidingWindowSteps
                    
                    thisMiWindowStart = slidingWindowStartPoint+slidingWindowStepSize*(slidingStepNum-1);
                    thisMiWindow = [thisMiWindowStart thisMiWindowStart+slidingWindowSize]; %TODO update each iter
                    
                    cv2Med= nan;
                    
                    
                    miPeriodSpikesInds = cellfun(@(x) find(x>(thisTrial.saccadeTime+thisMiWindow(1)) & x<(thisTrial.saccadeTime+thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
                    
                    if strcmpi(measureName,'CV2')
                        
                        if length([miPeriodSpikesInds{:}])>=2
                            
                            cv2Med = nanmedian(thisTrialCV2([miPeriodSpikesInds{:}]));
                        end
                        
                        
                        thisNeuronResults(slidingStepNum,trialNum) = cv2Med;
                        
                    elseif strcmpi(measureName,'rate')
                        
                        if ~isempty(miPeriodSpikesInds)
                            thisNeuronResults(slidingStepNum,trialNum) = length([miPeriodSpikesInds{:}]);
                            
                        end
                    end
                end
            end
        end
        
        
        for slidingStepNum = 1:numSlidingWindowSteps
            %TODO need to remove any nan from both labels and results
            %after building the matrix of trials then do MI caclulation
            %for all timeWindows
            goodTrials = find(~isnan(thisNeuronResults(slidingStepNum,:)));
            
            if ~isempty(goodTrials)
                [thisMI thisConfLims thisSigLevel] = miAndConfLims(thisNeuronResults(slidingStepNum,goodTrials),thisNeuronLabels(2,goodTrials),numMontCarlIter);
            else
                thisMI = nan;
                thisConfLims = [nan nan];
                thisSigLevel = nan;
            end
            thisNeuronMIs(slidingStepNum,1) = thisMI;
            thisNeuronMIs(slidingStepNum,2) = thisConfLims(1);
            thisNeuronMIs(slidingStepNum,3) = thisConfLims(2);
            thisNeuronMIs(slidingStepNum,4) = thisSigLevel;
            
        end
        thisWindowSizeResults{neuronNum,1} = thisNeuronResults;
        thisWindowSizeResults{neuronNum,2} = thisNeuronMIs;
        
        
%         for slidingStepNum = 1:numSlidingWindowSteps
%             
%             thisMiWindowStart = slidingWindowStartPoint+slidingWindowStepSize*(slidingStepNum-1);
%             thisMiWindow = [thisMiWindowStart thisMiWindowStart+slidingWindowSize]; %TODO update each iter
%             % slMi = cell(numTrials,1);
%             
%             
%             %go through each trail and get all spikes that fall within given window
%             for trialNum = 1:numTrials
%                 
%                 
%                 
%                 thisTrial = theseTrials(trialsToAnalyse(trialNum));
%                 
%                 
%                 miPeriodSpikes = cellfun(@(x) x(x>(thisTrial.saccadeTime+thisMiWindow(1)) & x<(thisTrial.saccadeTime+thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
%                 
%                 miSpikeRate = length([miPeriodSpikes{:}])/slidingWindowSize;
%                 
%                 if strcmpi(measureName,'CV2')
%                     
%                     thisNeuronResults(slidingStepNum,trialNum) = nanmedian(calculateCV2([miPeriodSpikes{:}]));
%                     
%                 else
%                 thisNeuronResults(slidingStepNum,trialNum) = miSpikeRate;
%                 
%                 end
%                 
%                 %after building the matrix of trials then do MI caclulation for this timeWindow
%                 
%             end
%             goodTrials = find(~isnan(thisNeuronResults(slidingStepNum,:)));
%             [thisMI thisConfLims thisSigLevel] = miAndConfLims(thisNeuronResults(slidingStepNum,goodTrials),thisNeuronLabels(2,goodTrials),numMontCarlIter);
%             
%             thisNeuronMIs(slidingStepNum,1) = thisMI;
%             thisNeuronMIs(slidingStepNum,2) = thisConfLims(1);
%             thisNeuronMIs(slidingStepNum,3) = thisConfLims(2);
%             thisNeuronMIs(slidingStepNum,4) = thisSigLevel;
%         end
%generate MI time vec and compensate for half window size so
            %that the results of each window is plotted in its midpoint
            %timewise.
            endPoint = slidingWindowStartPoint+slidingWindowStepSize*numSlidingWindowSteps;
            miTimeVec = linspace(slidingWindowStartPoint,endPoint,numSlidingWindowSteps);
            
            miTimeVec = miTimeVec+slidingWindowSize/2;
             thisWindowSizeResults{neuronNum,3} = miTimeVec;        


        %now plot rate results
        hF = figure;
        hAx1 = subplot(2,1,1,'parent',hF);
        plot(miTimeVec,thisNeuronResults,'parent',hAx1)
        hold on;
        plot(miTimeVec,mean(thisNeuronResults,2),'k','parent',hAx1) %TODO should plot mean for pro and anti seperately
        condCount = 1;
        for condCode = trialConditionCodes(:,1)'
            thisMu = mean(thisNeuronResults(:,thisNeuronLabels(2,:)==condCode),2);
            plot(miTimeVec,thisMu','color',cc(condCount,:),'parent',hAx1,'linewidth',2)
            condCount = condCount+1;
        end
        %plot(mean(thisNeuronResults(:,thisNeuronLabels(1,:)==2),2),'r','parent',hAx1)
        
        hAx2 = subplot(2,1,2,'parent',hF)
        plot(miTimeVec,thisNeuronMIs(:,1),'k')
        hold on
        plot(miTimeVec,thisNeuronMIs(:,2),'k-.')
        plot(miTimeVec,thisNeuronMIs(:,3),'k-.')
        
        suptitle([neuronList(neuronNum).neuronName 'MI of ' measureName ' over time for direction, aligned to saccade, ' num2str(1000*slidingWindowSize) 'ms window'])
        export_fig(hF,'-pdf','-append',miSlidSumFigName)
        delete(hF)
        
         thisNeuronSigMI = thisNeuronMIs(:,1);
         thisNeuronSigMI(thisNeuronMIs(:,1)<thisNeuronMIs(:,3)) = nan;
        %now keep the significant bits of the MI trace and NaN the rest.
        allNeuronSignificantMI(neuronNum,:) = thisNeuronSigMI;
        %         %now do exactly the same but when aligned to trialStart which is at
        %         %fixation
        %         for slidingStepNum = 1:numSlidingWindowSteps
        %
        %             thisMiWindowStart = slidingWindowStartPoint+slidingWindowStepSize*(slidingStepNum-1);
        %             thisMiWindow = [thisMiWindowStart thisMiWindowStart+slidingWindowSize]; %TODO update each iter
        %             % slMi = cell(numTrials,1);
        %
        %
        %             %go through each trail and get all spikes that fall within given window
        %             for trialNum = 1:numTrials
        %
        %
        %
        %                 thisTrial = theseTrials(trialsToAnalyse(trialNum));
        %
        %
        %                 miPeriodSpikes = cellfun(@(x) x(x>(thisMiWindow(1)) & x<(thisMiWindow(2))),thisTrial.alignedSpikes(thisNeuronUnitNum),'uniformOutput',false);
        %
        %                 miSpikeRate = length([miPeriodSpikes{:}])/slidingWindowSize;
        %                 thisNeuronResults(slidingStepNum,trialNum) = miSpikeRate;
        %
        %
        %
        %                 %after building the matrix of trials then do MI caclulation for this timeWindow
        %
        %             end
        %             [thisMI thisConfLims thisSigLevel] = miAndConfLims(thisNeuronResults(slidingStepNum,:),thisNeuronLabels(1,:),numMontCarlIter);
        %
        %             thisNeuronMIs(slidingStepNum,1) = thisMI;
        %             thisNeuronMIs(slidingStepNum,2) = thisConfLims(1);
        %             thisNeuronMIs(slidingStepNum,3) = thisConfLims(2);
        %             thisNeuronMIs(slidingStepNum,4) = thisSigLevel;
        %         end
        %
        %
        %         %now plot rate results
        %         hF = figure;
        %         hAx1 = subplot(2,1,1,'parent',hF);
        %         plot(thisNeuronResults,'parent',hAx1)
        %         hold on;
        %         plot(mean(thisNeuronResults,2),'k','parent',hAx1) %TODO should plot mean for pro and anti seperately
        %         plot(mean(thisNeuronResults(:,thisNeuronLabels(1,:)==1),2),'g','parent',hAx1)
        %         plot(mean(thisNeuronResults(:,thisNeuronLabels(1,:)==2),2),'r','parent',hAx1)
        %
        %         hAx2 = subplot(2,1,2,'parent',hF)
        %         plot(thisNeuronMIs(:,1),'k')
        %         hold on
        %         plot(thisNeuronMIs(:,2),'k-.')
        %         plot(thisNeuronMIs(:,3),'k-.')
        %
        %         suptitle([neuronList(neuronNum).neuronName 'MI over time, aligned to fixation, ' num2str(1000*slidingWindowSize) 'ms window'])
        %         export_fig(hF,'-pdf','-append',miSlidSumFigName)
        %         delete(hF)
        %TODO need to get the times that MI is significant
    end
    %TODO now plot all significant parts , one colour for verm 1 for lat
   
     save([resultsDir filesep 'slidingWidnowMIbyDir_' num2str(slidingWindowSize) 'msWin_' measureName '.mat'],'thisWindowSizeResults','slidingWindowSize')
       
    
    
    hF2 = figure;
    hAx1 = subplot(2,1,1,'parent',hF2)
    plot(miTimeVec,allNeuronSignificantMI(latNeurons,:),'color','r','parent',hAx1)
    hold on
    plot(miTimeVec,allNeuronSignificantMI(vermNeurons,:),'color','b','parent',hAx1)
    hAx2 = subplot(2,1,2,'parent',hF2)
    plot(miTimeVec,sum(~isnan(allNeuronSignificantMI(vermNeurons,:)))./length(vermNeurons),'color','b','linewidth',2,'parent',hAx2)
    hold on
    plot(miTimeVec,sum(~isnan(allNeuronSignificantMI(latNeurons,:)))./length(latNeurons),'color','r','linewidth',2,'parent',hAx2)
    legend('verm','lateral')
     suptitle(['All significant time periods by direction '  num2str(1000*slidingWindowSize) 'ms window'])
       
     export_fig(hF2,'-pdf','-append',miSlidSumFigName)
        delete(hF2)
        
end
end


%% QUICK LOAD SLID WIn MI
if doSlidWinMI

miTimeVec = thisWindowSizeResults{recordingsToAnalyse(1),3};
allNeuronSignificantMI=nan(max(recordingsToAnalyse),length(miTimeVec));
for neuronNum = recordingsToAnalyse
    
    thisNeuronMIs = thisWindowSizeResults{neuronNum,2};
    thisNeuronSigMI = thisNeuronMIs(:,1);
    thisNeuronSigMI = smooth(thisNeuronSigMI,15);
         thisNeuronSigMI(thisNeuronMIs(:,1)<thisNeuronMIs(:,3)) = nan;
        %now keep the significant bits of the MI trace and NaN the rest.
        allNeuronSignificantMI(neuronNum,:) = thisNeuronSigMI;
    
end
hF2 = figure;
    hAx1 = subplot(2,1,1,'parent',hF2)
    plot(miTimeVec,allNeuronSignificantMI(latNeurons,:),'color','r','parent',hAx1)
    hold on
    plot(miTimeVec,allNeuronSignificantMI(vermNeurons,:),'color','b','parent',hAx1)
    hAx2 = subplot(2,1,2,'parent',hF2)
    plot(miTimeVec,sum(~isnan(allNeuronSignificantMI(vermNeurons,:))),'color','b','linewidth',2,'parent',hAx2)
    hold on
    plot(miTimeVec,sum(~isnan(allNeuronSignificantMI(latNeurons,:))),'color','r','linewidth',2,'parent',hAx2)
    legend('verm','lateral')
     suptitle(['All significant time periods by dir by rate'  num2str(1000*slidingWindowSize) 'ms window'])
       set([hAx1 hAx2],'xlim',[-1 1])
     export_fig(hF2,'-pdf','-append',[resultsDir filesep 'sigPointsMiTest2015.pdf'])
end
%% MAX/MIN extraction


maxMinResFileName = [resultsDir filesep '-maxMinResults.mat'];
%[alignedTrials muRate] = alignAndExtractRates(trialStructure,selectedTrials,alignmentPoint,varargin)
if doMaxMin
 
%this uses the same code from the population response builder but instead
%will align all trials to saccade and then take 
 %trialTestTimes.saccPeriod = [-0.025 0.25]; %this is relative to saccade onset
%around saccade, and for each trial get the max,min and mean of the spike
%density function between those two times.

%we want to do 2 types of analysis:
% 1 for each trial get the time and value of the maximum/minumum of the sdf during
% the saccade period TODO do we need to limit to only peaks that are a few
% sds away?
% 2 get the time of the minimum in the mean sdf and then take the value of
% all trials at this time point
timeWindow = [4 4]; %for sdf traces which are sampled at 1KHz
saccadePerdiodInds = 1000*([-0.025 0.25]+timeWindow(1)); %inds in sdf trace to look for min and max

basePeriodInds = 1000*([-0.4 0]+timeWindow(1));
%sortByString = 'saccadeAmplitude';
%if exist(popResFileName,'file')~=2
stInit = cell(max(recordingsToAnalyse),2);
rangeRes = struct('alignedTrials',stInit,'maxMin',stInit);
 
for neuronNum = recordingsToAnalyse;
    thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
    unitCount =1;
    for unitNum  = thisNeuronUnits;
        %get the trials and the numbers of stable ones
        theseTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
        theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        
        theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
        theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
        muSdf = struct('Anti',[],'Pro',[]);
        %size(vertcat(theseStableTrials.alignedSpikes),2)<2;
        [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,theseAntiTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
        [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseProTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
        
        
        %also need trials aligned to fixation to get the baseline rate and
        %sd
        [baseTrials.Anti, baseMuSdf.Anti] = alignAndExtractRates(theseStableTrials,theseAntiTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow);
        [baseTrials.Pro, baseMuSdf.Pro] = alignAndExtractRates(theseStableTrials,theseProTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow);
        
        thisNeuronMaxMin = struct('Pro',[],'Anti',[]);
        %second method is find the minimum and maximum of the averaged
        %trace and get their time points
        for dirLabelNum = 1:2;
            thisLabel = byDirLabel{dirLabelNum,2};
            
            muTrace = muSdf.(thisLabel);
            
            muTraceSection = muTrace(saccadePerdiodInds(1):saccadePerdiodInds(2));
            
            [nadirVal nadirInd] = min(muTraceSection);
            [zenithVal zenithInd] = max(muTraceSection);
            
            nadirInd = nadirInd+saccadePerdiodInds(1);
            zenithInd = zenithInd+saccadePerdiodInds(1);
            %now stick all trials together
            trialTraces = vertcat(alignedTrials.(thisLabel).sdf); %all aligned trials, sdfs
            
            theseTrialValues.nadir = trialTraces(:,nadirInd);
            theseTrialValues.zenith = trialTraces(:,zenithInd);
             theseTrialValues.nadirInd = nadirInd;
            theseTrialValues.zenithInd = zenithInd;
            %now loop over these trials to get the trial by trial values
            theseBaseTrials = vertcat(baseTrials.(thisLabel).sdf);
            numTrials = size(trialTraces,1);
            stInit = cell(numTrials,1);
            trialByTrial = struct('baseMu',stInit,'baseSd',stInit,'lowVals',stInit,...
                'lowTimes',stInit,'highVals',stInit,'highTimes',stInit);
            for trialNum = 1:numTrials
               %get the mean and sd of the baseline period for this trial
                baseTraceSection = theseBaseTrials(trialNum,basePeriodInds(1):basePeriodInds(2));
                trialByTrial(trialNum).baseMu = mean(baseTraceSection);
                trialByTrial(trialNum).baseSd = std(baseTraceSection);
                
                
                %now get saccade period section and get the max and min
                %times and values
                thisTrialSaccPeriodTrace = trialTraces(trialNum,saccadePerdiodInds(1):saccadePerdiodInds(2));
                
                [highVal, highTime] = max(thisTrialSaccPeriodTrace);
                [lowVal, lowTime] = min(thisTrialSaccPeriodTrace);
                
                highTime = highTime+saccadePerdiodInds(1)-timeWindow(1)*1000; %make time relative to saccade onset
                lowTime = lowTime+saccadePerdiodInds(1)-timeWindow(1)*1000;
                
                
                trialByTrial(trialNum).lowVals = lowVal;
                trialByTrial(trialNum).lowTimes = lowTime;
                trialByTrial(trialNum).highVals = highVal;
                trialByTrial(trialNum).highTimes = highTime;
                
                %TODO we need to only select those peaks and troughs that
                %are above 2zscores from baseline rate and then count the
                %number of these
            end
           
            
            theseTrialValues.trialByTrial = trialByTrial;
            
            thisNeuronMaxMin.(thisLabel) = theseTrialValues;
        end
        
        
        thisPdfName = [resultsDir filesep maxMinFigNameEnd];
        
        %now plot the type 2 analysis results
        %proZenith v antiZenith, proNadir v antiNadir
        proZenithMu = mean(thisNeuronMaxMin.Pro.zenith);
        proZenithSem = std(thisNeuronMaxMin.Pro.zenith)/sqrt(length(thisNeuronMaxMin.Pro.zenith));
        proNadirMu = mean(thisNeuronMaxMin.Pro.nadir);
        proNadirSem = std(thisNeuronMaxMin.Pro.nadir)/sqrt(length(thisNeuronMaxMin.Pro.nadir));
         antiZenithMu = mean(thisNeuronMaxMin.Anti.zenith);
        antiZenithSem = std(thisNeuronMaxMin.Anti.zenith)/sqrt(length(thisNeuronMaxMin.Anti.zenith));
        antiNadirMu = mean(thisNeuronMaxMin.Anti.nadir);
        antiNadirSem = std(thisNeuronMaxMin.Anti.nadir)/sqrt(length(thisNeuronMaxMin.Anti.nadir));
        
        
        erIn = [proZenithSem antiZenithSem; proNadirSem antiNadirSem];
        cc = [0 1 0; 1 0 0; 0 1 0; 1 0 0];
        hF = figure;
        hAx = axes('parent',hF);
        [hPatches] = createGroupedBarGraph([proZenithMu antiZenithMu; proNadirMu antiNadirMu],erIn,hAx,cc,0.2);
        set(hAx,'xtick',[1 2],'xticklabel',{'Peak','Trough'})
        set(get(hAx,'ylabel'),'string','Firing Frequency (Hz)')
        set(get(hAx,'title'),'string',[neuronList(neuronNum).neuronName ' From mean trace,' num2str(unitNum)])
        
         
        try
            export_fig(hF,'-pdf','-append',thisPdfName)
        catch
            export_fig(hF,'-pdf',thisPdfName)
        end
        delete(hF)
        %now plot type 1 analysis which is the distribution of peak times,
        %trough times, peak values and trough values
        proHighVals = [thisNeuronMaxMin.Pro.trialByTrial.highVals];
        proLowVals = [thisNeuronMaxMin.Pro.trialByTrial.lowVals];
        proHighMu = mean(proHighVals);
        proLowMu = mean(proLowVals);
        proHighSem = std(proHighVals)/sqrt(length(proHighVals));
        proLowSem = std(proLowVals)/sqrt(length(proLowVals));
        %same for anti
        antiHighVals = [thisNeuronMaxMin.Anti.trialByTrial.highVals];
        antiLowVals = [thisNeuronMaxMin.Anti.trialByTrial.lowVals];
        antiHighMu = mean(antiHighVals);
        antiLowMu = mean(antiLowVals);
        antiHighSem = std(antiHighVals)/sqrt(length(antiHighVals));
        antiLowSem = std(antiLowVals)/sqrt(length(antiLowVals));
        
        erIn = [proHighSem antiHighSem; proLowSem antiLowSem];
        hF = figure;
        hAx = axes('parent',hF);
        [hPatches] = createGroupedBarGraph([proHighMu antiHighMu; proLowMu antiLowMu],erIn,hAx,cc,0.2);
        set(hAx,'xtick',[1 2],'xticklabel',{'Peak','Trough'})
        set(get(hAx,'ylabel'),'string','Firing Frequency (Hz)')
        set(get(hAx,'title'),'string',[neuronList(neuronNum).neuronName ' From individual traces, ' num2str(unitNum)])
        
        export_fig(hF,'-pdf','-append',thisPdfName)
         delete(hF)
        %now do exactly the same thing for teh times
        proHighTimes = [thisNeuronMaxMin.Pro.trialByTrial.highTimes];
        proLowTimes = [thisNeuronMaxMin.Pro.trialByTrial.lowTimes];
        proHighTimeMu = mean(proHighTimes);
        proLowTimeMu = mean(proLowTimes);
        proHighTimeSem = std(proHighTimes)/sqrt(length(proHighTimes));
        proLowTimeSem = std(proLowTimes)/sqrt(length(proLowTimes));
        %same for anti
        antiHighTimes = [thisNeuronMaxMin.Anti.trialByTrial.highTimes];
        antiLowTimes = [thisNeuronMaxMin.Anti.trialByTrial.lowTimes];
        antiHighTimeMu = mean(antiHighTimes);
        antiLowTimeMu = mean(antiLowTimes);
        antiHighTimeSem = std(antiHighTimes)/sqrt(length(antiHighTimes));
        antiLowTimeSem = std(antiLowTimes)/sqrt(length(antiLowTimes));
        
        erIn = [proHighTimeSem antiHighTimeSem; proLowTimeSem antiLowTimeSem];
        hF = figure;
        hAx = axes('parent',hF);
        [hPatches] = createGroupedBarGraph([proHighTimeMu antiHighTimeMu; proLowTimeMu antiLowTimeMu],erIn,hAx,cc,0.2);
        set(hAx,'xtick',[1 2],'xticklabel',{'Peak Time','Trough Time'})
        set(get(hAx,'ylabel'),'string','Event Time (relative to sacc onset)')
        set(get(hAx,'title'),'string',[neuronList(neuronNum).neuronName ' From individual traces,' num2str(unitNum)])
         export_fig(hF,'-pdf','-append',thisPdfName)
         delete(hF)
%         for dirLabelNum = 1:2;
%             thisLabel = byDirLabel{dirLabelNum,2};
%             [hA] = createSdfHeatMap(alignedTrials.(thisLabel),sortByString,[],'frRange',frDispRange(unitNum,:));
%             set(get(hA,'title'),'string',[neuronList(neuronNum).neuronName ' unit ' num2str(unitNum) ' ' thisLabel ' trials by ' sortByString])
%             %set(hA,'xlim',[3500 4500],'xtick',[3500 4000 4500],'xticklabel',{'-0.5','0','0.5'})
%             set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
%             
%             try
%                 export_fig(gcf,'-pdf','-append',expRend,popResPdfName)
%             catch
%                 export_fig(gcf,'-pdf',expRend,popResPdfName)
%             end
%         end
         close all
%thisNeuronMaxMin
        rangeRes(neuronNum,unitCount).alignedTrials = alignedTrials;
        rangeRes(neuronNum,unitCount).maxMin = thisNeuronMaxMin;
      
        unitCount= unitCount+1;
    end
end

save(maxMinResFileName,'rangeRes')
end
%% MAX/MIN BY DIRECRTION - same as above but now should do it for each target direction individually
%[alignedTrials muRate] = alignAndExtractRates(trialStructure,selectedTrials,alignmentPoint,varargin)

%this uses the same code from the population response builder but instead
%will align all trials to saccade and then take
%trialTestTimes.saccPeriod = [-0.025 0.25]; %this is relative to saccade onset
%around saccade, and for each trial get the max,min and mean of the spike
%density function between those two times.
if doMaxMinByDir
numDirections = 8;
validDirectionRows = find(~isnan(trialConditionCodes(:,1)))'; %which rows in trialConditonCodes are the 8 valid directions

%we want to do 2 types of analysis:
% 1 for each trial get the time and value of the maximum/minumum of the sdf during
% the saccade period TODO do we need to limit to only peaks that are a few
% sds away?
% 2 get the time of the minimum in the mean sdf and then take the value of
% all trials at this time point
timeWindow = [4 4]; %for sdf traces which are sampled at 1KHz
saccadePerdiodInds = 1000*([-0.025 0.25]+timeWindow(1)); %inds in sdf trace to look for min and max

basePeriodInds = 1000*([-0.4 0]+timeWindow(1));
%sortByString = 'saccadeAmplitude';
%if exist(popResFileName,'file')~=2
stInit = cell(max(recordingsToAnalyse),2);
rangeRes = struct('alignedTrials',stInit,'maxMin',stInit);

for neuronNum = recordingsToAnalyse;
    thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
    unitCount =1;
    for unitNum  = thisNeuronUnits(1); %only simple spikes for now
        %
        %thisUnitAllDirections = cell(1,8);
        stInit = cell(1,8);
        thisDirectionPeakTroughStats = struct('proZenithMu',stInit,'proZenithSem',stInit,'proNadirMu',stInit,'proNadirSem',[],...
                'antiZenithMu',stInit,'antiZenithSem',stInit,'antiNadirMu',[],'antiNadirSem',stInit,...
                'proHighMu',stInit,'proLowMu',stInit,'proHighSem',stInit,'proLowSem',stInit,...
                'antiHighMu',stInit,'antiLowMu',stInit,'antiHighSem',stInit,'antiLowSem',stInit,...
                'proHighTimeMu',stInit,'proLowTimeMu',stInit,'proHighTimeSem',stInit,'proLowTimeSem',stInit,...
                'antiHighTimeMu',stInit,'antiLowTimeMu',stInit,'antiHighTimeSem',stInit,'antiLowTimeSem',stInit,...
                'percValidProHigh',stInit,'percValidProLow',stInit,'percValidAntiHigh',stInit,'percValidAntiLow',stInit);
        for directionNum = 1:numDirections
            thisProCode = trialConditionCodes(validDirectionRows(directionNum),1);
            thisAntiCode = trialConditionCodes(validDirectionRows(directionNum),2);
            
            %4 lines form raster section to base trial finding on
            %         thisConditionCode = trialConditionCodes(plNum,dirLabelNum)
            %         tempLogCell = cellfun(@(x) x==thisConditionCode,{wholeNeuronResults(neuronNum).allStableTrials.conditionCode},'uniformoutput',false);
            %         theseDirTrialNums = find([tempLogCell{:}]);
            %         corDirTrialNums = intersect(wholeNeuronResults(neuronNum).selectedTrials.corSacTrials,theseDirTrialNums);
            %
            
            
            
            %get the trials and the numbers of stable ones
            theseCorTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
            theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
            
            tempLogCell = cellfun(@(x) x==thisProCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirProTrialNums = find([tempLogCell{:}]);
            thisDirCorProTrialNums = intersect(theseCorTrialNums,thisDirProTrialNums);
            
            tempLogCell = cellfun(@(x) x==thisAntiCode,{theseStableTrials.conditionCode},'uniformoutput',false);
            thisDirAntiTrialNums = find([tempLogCell{:}]);
            thisDirCorAntiTrialNums = intersect(theseCorTrialNums,thisDirAntiTrialNums);
            
            
            %         theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
            %         theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
            muSdf = struct('Anti',[],'Pro',[]);
            %size(vertcat(theseStableTrials.alignedSpikes),2)<2;
            [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,thisDirCorAntiTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
            [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,thisDirCorProTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
            
            
            %also need trials aligned to fixation to get the baseline rate and
            %sd
            [baseTrials.Anti, baseMuSdf.Anti] = alignAndExtractRates(theseStableTrials,thisDirCorAntiTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow);
            [baseTrials.Pro, baseMuSdf.Pro] = alignAndExtractRates(theseStableTrials,thisDirCorProTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow);
            
            thisDirMaxMin = struct('Pro',[],'Anti',[]);
            
            
             
            
            
            %second method is find the minimum and maximum of the averaged
            %trace and get their time points
            for dirLabelNum = 1:2;
                thisLabel = byDirLabel{dirLabelNum,2};
                
                muTrace = muSdf.(thisLabel);
                
                muTraceSection = muTrace(saccadePerdiodInds(1):saccadePerdiodInds(2));
                
                [nadirVal nadirInd] = min(muTraceSection);
                [zenithVal zenithInd] = max(muTraceSection);
                
                nadirInd = nadirInd+saccadePerdiodInds(1);
                zenithInd = zenithInd+saccadePerdiodInds(1);
                %now stick all trials together
                trialTraces = vertcat(alignedTrials.(thisLabel).sdf); %all aligned trials, sdfs
                
                theseTrialValues.nadir = trialTraces(:,nadirInd);
                theseTrialValues.zenith = trialTraces(:,zenithInd);
                theseTrialValues.nadirInd = nadirInd;
                theseTrialValues.zenithInd = zenithInd;
                %now loop over these trials to get the trial by trial values
                theseBaseTrials = vertcat(baseTrials.(thisLabel).sdf);
                numTrials = size(trialTraces,1);
                stInit = cell(numTrials,1);
                trialByTrial = struct('baseMu',stInit,'baseSd',stInit,'lowVals',stInit,...
                    'lowTimes',stInit,'highVals',stInit,'highTimes',stInit);
                for trialNum = 1:numTrials
                    %get the mean and sd of the baseline period for this trial
                    baseTraceSection = theseBaseTrials(trialNum,basePeriodInds(1):basePeriodInds(2));
                    trialByTrial(trialNum).baseMu = mean(baseTraceSection);
                    trialByTrial(trialNum).baseSd = std(baseTraceSection);
                    
                    
                    %limit for significant max/min is set at mu+-2*sd of
                    %baseline rate
                    lowThresh = trialByTrial(trialNum).baseMu-2*trialByTrial(trialNum).baseSd;
                    highThresh = trialByTrial(trialNum).baseMu+2*trialByTrial(trialNum).baseSd;
                    
                    %now get saccade period section and get the max and min
                    %times and values
                    thisTrialSaccPeriodTrace = trialTraces(trialNum,saccadePerdiodInds(1):saccadePerdiodInds(2));
                    
                    [highVal, highTime] = max(thisTrialSaccPeriodTrace);
                    [lowVal, lowTime] = min(thisTrialSaccPeriodTrace);
                    
                    if highVal>highThresh
                       
                        highTime = highTime+saccadePerdiodInds(1)-timeWindow(1)*1000; %make time relative to saccade onset
                    else
                         highVal = nan;
                         highTime = nan;
                    end
                    
                    if lowVal<lowThresh
                     lowTime = lowTime+saccadePerdiodInds(1)-timeWindow(1)*1000;
                    else
                        lowVal = nan;
                        lowTime = nan;
                    end
                    trialByTrial(trialNum).lowVals = lowVal;
                    trialByTrial(trialNum).lowTimes = lowTime;
                    trialByTrial(trialNum).highVals = highVal;
                    trialByTrial(trialNum).highTimes = highTime;
                    
                    %TODO we need to only select those peaks and troughs that
                    %are above 2zscores from baseline rate and then count the
                    %number of these
                end
                
                
                theseTrialValues.trialByTrial = trialByTrial;
                
                thisDirMaxMin.(thisLabel) = theseTrialValues;
            end
            
            
           % thisPdfName = [resultsDir filesep maxMinFigNameEnd];
            
            %now plot the type 2 analysis results
            %proZenith v antiZenith, proNadir v antiNadir
            proZenithMu = mean(thisDirMaxMin.Pro.zenith);
            proZenithSem = std(thisDirMaxMin.Pro.zenith)/sqrt(length(thisDirMaxMin.Pro.zenith));
            proNadirMu = mean(thisDirMaxMin.Pro.nadir);
            proNadirSem = std(thisDirMaxMin.Pro.nadir)/sqrt(length(thisDirMaxMin.Pro.nadir));
            antiZenithMu = mean(thisDirMaxMin.Anti.zenith);
            antiZenithSem = std(thisDirMaxMin.Anti.zenith)/sqrt(length(thisDirMaxMin.Anti.zenith));
            antiNadirMu = mean(thisDirMaxMin.Anti.nadir);
            antiNadirSem = std(thisDirMaxMin.Anti.nadir)/sqrt(length(thisDirMaxMin.Anti.nadir));
            
            
            %now copy all to struct
            thisDirectionPeakTroughStats(directionNum).proZenithMu = proZenithMu;
            thisDirectionPeakTroughStats(directionNum).proZenithSem = proZenithSem;
            thisDirectionPeakTroughStats(directionNum).proNadirMu = proNadirMu;
            thisDirectionPeakTroughStats(directionNum).proNadirSem = proNadirSem;
            thisDirectionPeakTroughStats(directionNum).antiZenithMu = antiZenithMu;
            thisDirectionPeakTroughStats(directionNum).antiZenithSem = antiZenithSem;
            thisDirectionPeakTroughStats(directionNum).antiNadirMu = antiNadirMu;
            thisDirectionPeakTroughStats(directionNum).antiNadirSem = antiNadirSem;
            
            
%             erIn = [proZenithSem antiZenithSem; proNadirSem antiNadirSem];
%             
%             hF = figure;
%             hAx = axes('parent',hF);
%             [hPatches] = createGroupedBarGraph([proZenithMu antiZenithMu; proNadirMu antiNadirMu],erIn,hAx,cc,0.2);
%             set(hAx,'xtick',[1 2],'xticklabel',{'Peak','Trough'})
%             set(get(hAx,'ylabel'),'string','Firing Frequency (Hz)')
%             set(get(hAx,'title'),'string',[neuronList(neuronNum).neuronName 'dir ' num2str(directionNum) ' From mean trace,' num2str(unitNum)])
%             
%             
%             try
%                 export_fig(hF,'-pdf','-append',thisPdfName)
%             catch
%                 export_fig(hF,'-pdf',thisPdfName)
%             end
%             delete(hF)
            %now plot type 1 analysis which is the distribution of peak times,
            %trough times, peak values and trough values
            proHighVals = [thisDirMaxMin.Pro.trialByTrial.highVals];
            proLowVals = [thisDirMaxMin.Pro.trialByTrial.lowVals];
            %find how many peaks or troughs were significant (nan if not)
            numValidProHighs = sum(~isnan(proHighVals));
            numValidProLows = sum(~isnan(proLowVals));
            
            %what if no valid trials?
            if numValidProLows>0               
                proLowMu = nanmean(proLowVals);
                proLowSem = nanstd(proLowVals)/sqrt(numValidProLows);
            else  
                proLowMu = nan;
                proLowSem = nan;
            end
            if numValidProHighs>0
                proHighMu = nanmean(proHighVals);
                proHighSem = nanstd(proHighVals)/sqrt(numValidProHighs);
            else
                proHighMu = nan;
                proHighSem = nan;
            end
            %same for anti
            antiHighVals = [thisDirMaxMin.Anti.trialByTrial.highVals];
            antiLowVals = [thisDirMaxMin.Anti.trialByTrial.lowVals];
            %find how many peaks or troughs were significant (nan if not)
            numValidAntiHighs = sum(~isnan(antiHighVals));
            numValidAntiLows = sum(~isnan(antiLowVals));
            
            
            
            if  numValidAntiLows>0
                antiLowMu = nanmean(antiLowVals);
                antiLowSem = nanstd(antiLowVals)/sqrt(numValidAntiLows);
            else
                antiLowMu = nan;
                antiLowSem = nan;
            end
            
            if numValidAntiHighs>0
                antiHighMu = nanmean(antiHighVals);
                antiHighSem = nanstd(antiHighVals)/sqrt(numValidAntiHighs);
            else
                antiHighMu = nan;
                antiHighSem = nan;
            end
            
            
            %now copy all to struct
            thisDirectionPeakTroughStats(directionNum).proHighMu = proHighMu;
            thisDirectionPeakTroughStats(directionNum).proLowMu = proLowMu;
            thisDirectionPeakTroughStats(directionNum).proHighSem = proHighSem;
            thisDirectionPeakTroughStats(directionNum).proLowSem = proLowSem;
            thisDirectionPeakTroughStats(directionNum).antiHighMu = antiHighMu;
            thisDirectionPeakTroughStats(directionNum).antiLowMu = antiLowMu;
            thisDirectionPeakTroughStats(directionNum).antiHighSem = antiHighSem;
            thisDirectionPeakTroughStats(directionNum).antiLowSem = antiLowSem;
            
            
            
%             erIn = [proHighSem antiHighSem; proLowSem antiLowSem];
%             hF = figure;
%             hAx = axes('parent',hF);
%             [hPatches] = createGroupedBarGraph([proHighMu antiHighMu; proLowMu antiLowMu],erIn,hAx,cc,0.2);
%             set(hAx,'xtick',[1 2],'xticklabel',{'Peak','Trough'})
%             set(get(hAx,'ylabel'),'string','Firing Frequency (Hz)')
%             set(get(hAx,'title'),'string',[neuronList(neuronNum).neuronName 'dir ' num2str(directionNum) ' From individual traces, ' num2str(unitNum)])
%             
%             export_fig(hF,'-pdf','-append',thisPdfName)
%             delete(hF)
            %now do exactly the same thing for teh times
            proHighTimes = [thisDirMaxMin.Pro.trialByTrial.highTimes];
            proLowTimes = [thisDirMaxMin.Pro.trialByTrial.lowTimes];
            
             if  numValidProLows>0
                proLowTimeMu = nanmean(proLowTimes);
                proLowTimeSem = nanstd(proLowTimes)/sqrt(numValidProLows);
            else
                proLowTimeMu = nan;
                proLowTimeSem = nan;
            end
            
            if numValidProHighs>0
                 proHighTimeMu = nanmean(proHighTimes);
                 proHighTimeSem = nanstd(proHighTimes)/sqrt(numValidProHighs);
            else
                proHighTimeMu = nan;
                proHighTimeSem = nan;
            end
            
           
           
           
           
            %same for anti
            antiHighTimes = [thisDirMaxMin.Anti.trialByTrial.highTimes];
            antiLowTimes = [thisDirMaxMin.Anti.trialByTrial.lowTimes];
             if  numValidAntiLows>0
                antiLowTimeMu = nanmean(antiLowTimes);
               antiLowTimeSem = nanstd(antiLowTimes)/sqrt(numValidAntiLows);
            else
               antiLowTimeMu = nan;
                antiLowTimeSem = nan;
            end
            
            if numValidAntiHighs>0
                 antiHighTimeMu = nanmean(antiHighTimes);
                 antiHighTimeSem = nanstd(antiHighTimes)/sqrt(numValidAntiHighs);
            else
                antiHighTimeMu = nan;
                antiHighTimeSem = nan;
            end
            
            
            
            
            
            
            %copy all to struct
            thisDirectionPeakTroughStats(directionNum).proHighTimeMu = proHighTimeMu;
            thisDirectionPeakTroughStats(directionNum).proLowTimeMu = proLowTimeMu;
            thisDirectionPeakTroughStats(directionNum).proHighTimeSem = proHighTimeSem;
            thisDirectionPeakTroughStats(directionNum).proLowTimeSem = proLowTimeSem;
            thisDirectionPeakTroughStats(directionNum).antiHighTimeMu = antiHighTimeMu;
            thisDirectionPeakTroughStats(directionNum).antiLowTimeMu = antiLowTimeMu;
            thisDirectionPeakTroughStats(directionNum).antiHighTimeSem = antiHighTimeSem;
            thisDirectionPeakTroughStats(directionNum).antiLowTimeSem = antiLowTimeSem;
            
            
            %now find the percentage of each type of vaild trial
            thisDirectionPeakTroughStats(directionNum).percValidProHigh = 100*numValidProHighs/length(proHighTimes);
            
            thisDirectionPeakTroughStats(directionNum).percValidProLow = 100*numValidProLows/length(proHighTimes);
            thisDirectionPeakTroughStats(directionNum).percValidAntiHigh = 100*numValidAntiHighs/length(antiLowTimes);
            thisDirectionPeakTroughStats(directionNum).percValidAntiLow = 100*numValidAntiLows/length(antiLowTimes);
            
            
%             erIn = [proHighTimeSem antiHighTimeSem; proLowTimeSem antiLowTimeSem];
%             hF = figure;
%             hAx = axes('parent',hF);
%             [hPatches] = createGroupedBarGraph([proHighTimeMu antiHighTimeMu; proLowTimeMu antiLowTimeMu],erIn,hAx,cc,0.2);
%             set(hAx,'xtick',[1 2],'xticklabel',{'Peak Time','Trough Time'})
%             set(get(hAx,'ylabel'),'string','Event Time (relative to sacc onset)')
%             set(get(hAx,'title'),'string',[neuronList(neuronNum).neuronName 'dir ' num2str(directionNum) ' From individual traces,' num2str(unitNum)])
%             export_fig(hF,'-pdf','-append',thisPdfName)
%             delete(hF)
            %         for dirLabelNum = 1:2;
            %             thisLabel = byDirLabel{dirLabelNum,2};
            %             [hA] = createSdfHeatMap(alignedTrials.(thisLabel),sortByString,[],'frRange',frDispRange(unitNum,:));
            %             set(get(hA,'title'),'string',[neuronList(neuronNum).neuronName ' unit ' num2str(unitNum) ' ' thisLabel ' trials by ' sortByString])
            %             %set(hA,'xlim',[3500 4500],'xtick',[3500 4000 4500],'xticklabel',{'-0.5','0','0.5'})
            %             set(hA,'xlim',[axisTicks(1) axisTicks(3)],'xtick',axisTicks,'xticklabel',axTickLabs)
            %
            %             try
            %                 export_fig(gcf,'-pdf','-append',expRend,popResPdfName)
            %             catch
            %                 export_fig(gcf,'-pdf',expRend,popResPdfName)
            %             end
            %         end
            close all
            %thisNeuronMaxMin
            rangeRes(neuronNum,unitCount).alignedTrials = alignedTrials;
            rangeRes(neuronNum,unitCount).maxMin(directionNum) = thisDirMaxMin;
            %thisUnitAllDirections{directionNum} = thisDirectionPeakTroughStats;
        end
        cc = [0 1 0; 1 0 0; 0 1 0; 1 0 0];
        anglesToTest = linspace(0,2*pi,9);
        anglesToTest = anglesToTest(1:8); %8 evenly spaced going counterclockwise
        statsToTest = {'proHigh','antiHigh';...
            'proLow','antiLow';...
            'proHighTime','antiHighTime';...
            'proLowTime','antiLowTime';...
            'proZenith','antiZenith';...
            'proNadir','antiNadir'}; %each row is comparison the columns are what to compare to each other
        numStatsToTest = size(statsToTest,1);
        for statToTest = 1:numStatsToTest
            
            %get data from struct for this unit, this neuron, all
            %directions
        theseProMus = [thisDirectionPeakTroughStats.([char(statsToTest{statToTest,1}) 'Mu'])];
        theseAntiMus = [thisDirectionPeakTroughStats.([char(statsToTest{statToTest,2}) 'Mu'])];
        theseProSems = [thisDirectionPeakTroughStats.([char(statsToTest{statToTest,1}) 'Sem'])];
        theseAntiSems = [thisDirectionPeakTroughStats.([char(statsToTest{statToTest,2}) 'Sem'])];
        hF = figure;
        hA = axes('parent',hF);
         %spider starts at 3 o'clock and goes anti-clockwise
            %plNum starts at 10 oclock and goes clockwise, vecotr below
            %converts
        polPlOrder = [5 3 2 1 4 6 7 8]; %this the order to take the directions in the struct to put them into the order for polar plotting
        proMusToPlot = theseProMus(polPlOrder);
        proSemsToPlot = theseProSems(polPlOrder);
        antiMusToPlot = theseAntiMus(polPlOrder);
        antiSemsToPlot = theseAntiSems(polPlOrder);
        
        allDataToPlot = {proMusToPlot antiMusToPlot};
        allErrsToPlot = {proSemsToPlot antiSemsToPlot};
        [hA] = createConfLimSpider(allDataToPlot,allErrsToPlot,anglesToTest,hA,[]);
        set(get(hA,'title'),'string',[neuronList(neuronNum).neuronName ',' num2str(unitNum) ', ' char(statsToTest{statToTest,1}) ' vs ' char(statsToTest{statToTest,2})])
            %spider starts at 3 o'clock and goes anti-clockwise
            %plNum starts at 10 oclock and goes clockwise
            %so flipud
%          ratesToPlot{dirLabelNum} = flipud(circshift(theseRatesToPlot,3));
%             errToPlot{dirLabelNum} = flipud(circshift(theseErrsToPlot,3));
%          [hA] = createConfLimSpider(ratesToPlot,errToPlot,anglesToTest,hA);
%         maxVal = max(max([ratesToPlot{:}]));
%         set(hA,'color','none','xlim',[-maxVal maxVal],'ylim',[-maxVal maxVal])
%         
%         
%         export_fig(hF2,'-pdf','-append',burstSumFigName)
%         

            export_fig(hF,'-pdf','-append',[resultsDir filesep maxMinByDirPdfName])
            
            %aslo export as a png for use in the chris style summary
            %figures, only for SS
            if unitNum==thisNeuronUnits(1)
            thisPngname = [resultsDir filesep neuronList(neuronNum).neuronName '-' char(statsToTest{statToTest,1}) ssDirTuningPngNameEnd];
            export_fig(hF,'-png',thisPngname)
            end
            delete(hF)
        end
        %now draw the 6 spider plots
        %mean trace peak
        
        %mean trace trough
        %individual peak
        %individual trough
        %individual peak time
        %individual trought time
        
        %now do a final one which plots
        
        unitCount= unitCount+1;
    end
end
end

%% INSTRUCTION PERIOD RASTERS

if doInstRasters
thisUnitNum  =1; %CS are marked on rasters anyway
for neuronNum = recordingsToAnalyse; %do 16 though end

    %build figure filename
    thisNeuronSumFigName = [resultsDir filesep neuronList(neuronNum).neuronName 'instructionPeriodRasters2.pdf'];
    
   
    %loop over all chosen rasters and create and export their figures
    for rastSelected = rastToPlot
        %also build png name
        thisNeuronThisRastPngName =  [resultsDir filesep neuronList(neuronNum).neuronName '-' (char(rastSelected)) instRastPngNameEnd]
        
        [hAx1 hAx2 hA2] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
            wholeNeuronResults(neuronNum).selectedTrials.(char(rastSelected)),[1 1],timeWindow,...
            'sortBy','bit6time','alignTo','goCueTime','markEvent','bit6time');
        title([neuronList(neuronNum).neuronName ' ' rastSelected ',sort by fixCue, align to goCue'])
        export_fig(gcf,'-pdf','-append',thisNeuronSumFigName)
        
      
         export_fig(gcf,'-png',thisNeuronThisRastPngName)
        
        delete(gcf)
%         
%         if ismember(neuronNum,exampleNeurons)
%             %export to eps
%            thisFigSaveName = [resultsDir filesep neuronList(neuronNum).neuronName '-' rastSelected{:} '.eps'];
%             export_fig(gcf,'-eps',thisFigSaveName)
%         end
%         
    end 
    
end
end
%% RASTER PLOTS
%create a raster pot (inclduing eye trace) for all conditions listed in
%cell below
%loop over all neurons and plot rasters for all, pro and anti trials. Then
%alos by direction
smallTimeWindow = [1 1];
if doRasters
thisUnitNum  =1; %CS are marked on rasters anyway
for neuronNum = recordingsToAnalyse; %do 16 though end
    
    
    %build figure filename
    thisNeuronSumFigName = [resultsDir filesep neuronList(neuronNum).neuronName indNeuronSumFigEnd];
    %loop over all chosen rasters and create and export their figures
    for rastSelected = rastToPlot
        
        %this raster is sorted yb go cue time as default
        [hAx1 hAx2 hA2] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
            wholeNeuronResults(neuronNum).selectedTrials.(char(rastSelected)),smallTimeWindow,timeWindow)
        title([neuronList(neuronNum).neuronName ' ' rastSelected])
        export_fig(gcf,'-pdf','-append',thisNeuronSumFigName)
        
        
        if ismember(neuronNum,exampleNeurons)
            %export to eps
           thisFigSaveName = [resultsDir filesep neuronList(neuronNum).neuronName '-' rastSelected{:} '.eps'];
            export_fig(gcf,'-eps',thisFigSaveName)
        end
        delete(gcf)
        
        thisNeuronThisRastPngName =  [resultsDir filesep neuronList(neuronNum).neuronName '-' (char(rastSelected)) unsortedRasterPngName]
        
        %now make a raster with no sroting and just keeping the original
        %trail order to see any cahgnes over time
          [hAx1 hAx2 hA2] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,...
            wholeNeuronResults(neuronNum).selectedTrials.(char(rastSelected)),smallTimeWindow,timeWindow,'sortBy','');
        
        
        title([neuronList(neuronNum).neuronName ' ' rastSelected])
       export_fig(gcf,'-png',thisNeuronThisRastPngName)
       delete(gcf)
    end
    
    
    %now alos produce rasters split into each condition and direction
    %seperately for pro saccades then anti
    
    for dirLabelNum = 1:2
       
        hF2 = figure('units','normalized','position',[-0.0278    1.2240    0.5146    0.9688]);
        limMat = nan(9,3); %to keep track of firing rate axis lims, last colum is handle
        %murRateByDir = nan(9,1); %to keep the saccade period firing rae for each direction
        for plNum = find(~isnan(trialConditionCodes(:,dirLabelNum)))'
            
            %set current figure and placce axis
            set(0, 'currentfigure', hF2);
            %hA = subaxis(3,3,plNum,'SpacingVert',0.05,'MR',0.05);
            hA = subaxis(3,3,plNum,'SpacingVert',0.2,'MR',0.05);
            set(hA,'position',get(hA,'position')-[0 0.075 0 0])
            hA2 = axes('XAxisLocation','bottom','YAxisLocation','right','color','none',...
                'Xcolor',[1 1 1],'Ycolor',[1 0 0],'XTick',[],'parent',hF2,'position',get(hA,'position'));
            hA3 = copyobj(hA,hF2);
            set(hA3,'position',get(hA3,'position')+[0 0.15 0 0])
            
            
            %now get all relevant trials
            thisConditionCode = trialConditionCodes(plNum,dirLabelNum)
            tempLogCell = cellfun(@(x) x==thisConditionCode,{wholeNeuronResults(neuronNum).allStableTrials.conditionCode},'uniformoutput',false);
            theseDirTrialNums = find([tempLogCell{:}]);
            corDirTrialNums = intersect(wholeNeuronResults(neuronNum).selectedTrials.corSacTrials,theseDirTrialNums);
            if ~isempty(corDirTrialNums)
            %create raster
            [hAx1 hAx2 hAx22] = rasterFromTrials(wholeNeuronResults(neuronNum).allStableTrials,corDirTrialNums,smallTimeWindow,timeWindow)
            %copy to new figure and make pretty
            newHand = copyobj(allchild(hAx2),hA);
            set(hA,'Xticklabel','','Yticklabel','')
            set(hA,'ylim',[0 length(corDirTrialNums)+1])
            newHand2 = copyobj(allchild(hAx22),hA2);
            if plNum ~=9
                set(hA2,'Xticklabel','','Yticklabel','')
                set(get(hA2,'ylabel'),'string','')
            end
            
            limMat(plNum,1:2) = get(hAx22,'ylim');
            limMat(plNum,3) = hA2;
            %now copy eye traces
            newHand = copyobj(allchild(hAx1),hA3);
            
            set(hA3,'Xticklabel','','Yticklabel','')
            set(hA3,'ylim',[-15 15])
            delete(get(hAx22,'parent'))
            end
            
            %         %now get the firing rate in the saccade period using the kunimatsu
            %         %anyslsis function
            %         [kuniResDir kuniResDirStats] = kunimatsuAnalysis(allStableTrials,corDirTrialNums,timeWindow,'',thisUnitNum);
            %         murRateByDir(plNum) = kuniResDirStats.muRates(4);
            %         delete(gcf)
        end
        set(limMat(~isnan(limMat(:,3)),3),'ylim',[min(limMat(:,1)) max(limMat(:,2))])
        set(findobj(hF2,'type','axes'),'xlim',[-smallTimeWindow(1) smallTimeWindow(2)])
        set(hA2,'ytickmode','auto','xticklabelmode','auto')
        suptitle([neuronList(neuronNum).neuronName ' ' byDirLabel{dirLabelNum,2}])
        %TODO automate for other window size choices
        set(get(hA2,'xlabel'),'string','Time (s)','color','k')
        set(hA,'xtick',[-0.2:0.2:0.4])
        set(hA,'xticklabel',{'-0.2','0','0.2', '0.4'})
        set(get(hA2,'ylabel'),'string','Firing Rate (Hz)')
        set(get(hA3,'ylabel'),'string','Eye Position')
        set(hA3,'ytick',[-15:15:15])
        set(hA3,'yticklabel',{'-15','15','15'})
        set(hA3,'yaxislocation','right')
        
        
        %now add spider plot to middle showing saccade period firing rate for
        %each direction
        theseMuRates = allNeuronResults(neuronNum,thisUnitNum).kuniResults.muRateByDir(:,dirLabelNum);
        hA = subaxis(3,3,5,'SpacingVert',0.2,'MR',0.05);
        ratesToPlot = theseMuRates(~isnan(theseMuRates));
        %spider starts at 3 o'clock and goes anti-clockwise
        %plNum starts at 10 oclock and goes clockwise
        %so flipud
        ratesToPlot = flipud(circshift(ratesToPlot,3));
        [f, ca, o] = createSpider(ratesToPlot ,'',[0 max(ratesToPlot)],repmat({''},8,1),'');
        hT = copyobj(ca,hF2);
        set(hT,'units','normalized','position',[get(hA,'position') + [0 -0.1 0 0.1]])
        
        
        delete(hA)
        if ismember(neuronNum,exampleNeurons)
            thisFigSaveName = [resultsDir filesep neuronList(neuronNum).neuronName '-' byDirLabel{dirLabelNum,2} '-rastByDir.eps'];
            export_fig(hF2,'-eps',thisFigSaveName)
        end
        export_fig(hF2,'-pdf','-append',[resultsDir filesep rasterByDirPdfName])
        
    end
    
    
    close all
    
end
end
%% CS COUNTS BY DIRECTION

% same as above but now should do it for each target direction individually
% but just counting CS no max/min stuff
%[alignedTrials muRate] = alignAndExtractRates(trialStructure,selectedTrials,alignmentPoint,varargin)

%this uses the same code from the population response builder but instead
%will align all trials to saccade and then take
%trialTestTimes.saccPeriod = [-0.025 0.25]; %this is relative to saccade
%onset,%now switched to 50:200 post saccadonset to match the herzfield
%paper, however that should actually be by saccade offset.

%around saccade, and for each trial get the max,min and mean of the spike
%density function between those two times.
if doCsByDir
csCountByDirPdfName = 'csCountsByDirection.pdf';
numDirections = 8;
validDirectionRows = find(~isnan(trialConditionCodes(:,1)))'; %which rows in trialConditonCodes are the 8 valid directions

%we want to do 2 types of analysis:
% 1 for each trial get the time and value of the maximum/minumum of the sdf during
% the saccade period TODO do we need to limit to only peaks that are a few
% sds away?
% 2 get the time of the minimum in the mean sdf and then take the value of
% all trials at this time point
timeWindow = [4 4]; %for sdf traces which are sampled at 1KHz
saccadePerdiodInds = 1000*([0.05 0.2]+timeWindow(1)); %inds in sdf trace to look for min and max

basePeriodInds = 1000*([-0.4 0]+timeWindow(1));
%sortByString = 'saccadeAmplitude';
%if exist(popResFileName,'file')~=2
stInit = cell(max(recordingsToAnalyse),2);
countRes = struct('alignedTrials',stInit,'maxMin',stInit);

for neuronNum = recordingsToAnalyse;
    thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
    unitCount =1;
    if length(thisNeuronUnits)>1
        for unitNum  = thisNeuronUnits(2); %only simple spikes for now
            %
            %thisUnitAllDirections = cell(1,8);
            thisDirectionCountStats = cell(1,8);
            %thisDirectionCountStats = struct('proCount',stInit,'antiCount',stInit);
            %                 'proHighMu',stInit,'proLowMu',stInit,'proHighSem',stInit,'proLowSem',stInit,...
            %                 'antiHighMu',stInit,'antiLowMu',stInit,'antiHighSem',stInit,'antiLowSem',stInit,...
            %                 'proHighTimeMu',stInit,'proLowTimeMu',stInit,'proHighTimeSem',stInit,'proLowTimeSem',stInit,...
            %                 'antiHighTimeMu',stInit,'antiLowTimeMu',stInit,'antiHighTimeSem',stInit,'antiLowTimeSem',stInit,...
            %                 'percValidProHigh',stInit,'percValidProLow',stInit,'percValidAntiHigh',stInit,'percValidAntiLow',stInit);
            for directionNum = 1:numDirections
                thisProCode = trialConditionCodes(validDirectionRows(directionNum),1);
                thisAntiCode = trialConditionCodes(validDirectionRows(directionNum),2);
                
                %4 lines form raster section to base trial finding on
                %         thisConditionCode = trialConditionCodes(plNum,dirLabelNum)
                %         tempLogCell = cellfun(@(x) x==thisConditionCode,{wholeNeuronResults(neuronNum).allStableTrials.conditionCode},'uniformoutput',false);
                %         theseDirTrialNums = find([tempLogCell{:}]);
                %         corDirTrialNums = intersect(wholeNeuronResults(neuronNum).selectedTrials.corSacTrials,theseDirTrialNums);
                %
                
                
                
                %get the trials and the numbers of stable ones
                theseCorTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
                theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
                
                tempLogCell = cellfun(@(x) x==thisProCode,{theseStableTrials.conditionCode},'uniformoutput',false);
                thisDirProTrialNums = find([tempLogCell{:}]);
                thisDirCorProTrialNums = intersect(theseCorTrialNums,thisDirProTrialNums);
                
                tempLogCell = cellfun(@(x) x==thisAntiCode,{theseStableTrials.conditionCode},'uniformoutput',false);
                thisDirAntiTrialNums = find([tempLogCell{:}]);
                thisDirCorAntiTrialNums = intersect(theseCorTrialNums,thisDirAntiTrialNums);
                
                
                %         theseAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
                %         theseProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
                % muSdf = struct('Anti',[],'Pro',[]);
                %size(vertcat(theseStableTrials.alignedSpikes),2)<2;
                [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,thisDirCorAntiTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
                [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,thisDirCorProTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
                
                countPeriod = [-0.025 0.25];
                %now count the number of complex spikes in the given bin
                [csCount.Anti] = countCsInBin(theseStableTrials,thisDirCorAntiTrialNums,unitNum,countPeriod);
                [csCount.Pro] = countCsInBin(theseStableTrials,thisDirCorProTrialNums,unitNum,countPeriod);
                %             %also need trials aligned to fixation to get the baseline rate and
                %             %sd
                %             [baseTrials.Anti, baseMuSdf.Anti] = alignAndExtractRates(theseStableTrials,thisDirCorAntiTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow);
                %             [baseTrials.Pro, baseMuSdf.Pro] = alignAndExtractRates(theseStableTrials,thisDirCorProTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow);
                %
                %  thisDirMaxMin = struct('Pro',[],'Anti',[]);
                
                thisDirectionCountStats{directionNum} = csCount;
            end
            close all
            countRes(neuronNum,unitCount).directionCsCounts = thisDirectionCountStats;
            %thisNeuronMaxMin
            % countRes(neuronNum,unitCount).alignedTrials = alignedTrials;
            % countRes(neuronNum,unitCount).maxMin(directionNum) = thisDirMaxMin;
            %thisUnitAllDirections{directionNum} = thisDirectionPeakTroughStats;
        end
        cc = [0 1 0; 1 0 0; 0 1 0; 1 0 0];
        anglesToTest = linspace(0,2*pi,9);
        anglesToTest = anglesToTest(1:8); %8 evenly spaced going counterclockwise
        %         statsToTest = {'proHigh','antiHigh';...
        %             'proLow','antiLow';...
        %             'proHighTime','antiHighTime';...
        %             'proLowTime','antiLowTime';...
        %             'proZenith','antiZenith';...
        %             'proNadir','antiNadir'}; %each row is comparison the columns are what to compare to each other
        %         numStatsToTest = size(statsToTest,1);
        %         for statToTest = 1:numStatsToTest
        
        %get data from struct for this unit, this neuron, all
        %directions
        theseCsCounts = [thisDirectionCountStats{:}];
        theseProCounts = {theseCsCounts.Pro};
        theseAntiCounts = {theseCsCounts.Anti};
        %each of these cells now is 1xnumDirection, each cell contains a 1,numTrials
        %vector of counts of cs in the bin specified
        csProbByDirection.Pro = cellfun(@(x) sum(x)/length(x),theseProCounts);
        csProbByDirection.Anti = cellfun(@(x) sum(x)/length(x),theseAntiCounts);
        %         theseProSems = [thisDirectionPeakTroughStats.([char(statsToTest{statToTest,1}) 'Sem'])];
        %         theseAntiSems = [thisDirectionPeakTroughStats.([char(statsToTest{statToTest,2}) 'Sem'])];
        hF = figure;
        hA = axes('parent',hF);
        %spider starts at 3 o'clock and goes anti-clockwise
        %plNum starts at 10 oclock and goes clockwise, vecotr below
        %converts
        polPlOrder = [5 3 2 1 4 6 7 8]; %this the order to take the directions in the struct to put them into the order for polar plotting
        proMusToPlot = csProbByDirection.Pro(polPlOrder);
        proSemsToPlot = zeros(1,8); %no error just count
        antiMusToPlot = csProbByDirection.Anti(polPlOrder);
        antiSemsToPlot = zeros(1,8);
        
        allDataToPlot = {proMusToPlot antiMusToPlot};
        allErrsToPlot = {proSemsToPlot antiSemsToPlot};
        [hA] = createConfLimSpider(allDataToPlot,allErrsToPlot,anglesToTest,hA,[]);
        set(get(hA,'title'),'string',[neuronList(neuronNum).neuronName ', CS probability in saccade period by direction'])
        %spider starts at 3 o'clock and goes anti-clockwise
        %plNum starts at 10 oclock and goes clockwise
        %so flipud
        %          ratesToPlot{dirLabelNum} = flipud(circshift(theseRatesToPlot,3));
        %             errToPlot{dirLabelNum} = flipud(circshift(theseErrsToPlot,3));
        %          [hA] = createConfLimSpider(ratesToPlot,errToPlot,anglesToTest,hA);
        %         maxVal = max(max([ratesToPlot{:}]));
        %         set(hA,'color','none','xlim',[-maxVal maxVal],'ylim',[-maxVal maxVal])
        %
        %
        %         export_fig(hF2,'-pdf','-append',burstSumFigName)
        %
        
        export_fig(hF,'-pdf','-append',[resultsDir filesep csCountByDirPdfName])
        
        %alos export as .png for Chris summary figures
        export_fig(hF,'-png',[resultsDir filesep neuronList(neuronNum).neuronName csDirTuningPngNameEnd])
        delete(hF)
        %         end
        %now draw the 6 spider plots
        %mean trace peak
        
        %mean trace trough
        %individual peak
        %individual trough
        %individual peak time
        %individual trought time
        
        unitCount= unitCount+1;
    end
    
end

end
%%  new spike densiyty plots as currently they come from teh population section which takes ages
%saccadeSdfPngNameEnds = {'Saccade-muSdfSSResort.png','Saccade-muSdfCSResort.png'}; %for each unit
%instructSdfPngNameEnds = {'Instruct-muSdfSSResort.png','Instruct-muSdfCSResort.png'}; %for each unit
%TODO for both the 
%1 confLimSdf for aligned to saccade onset for both units, pro and anti
%separate colours. expand the axis to go -1s
%2 aligned to fix Cue
allNeuronsSdf.Pro = cell(max(recordingsToAnalyse),1);
allNeuronsSdf.Anti = cell(max(recordingsToAnalyse),1);
for neuronNum = recordingsToAnalyse
    disp(['Num ' num2str(neuronNum)])
     thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
    unitCount =1;
    for unitNum  = thisNeuronUnits
  
        %get the trials and the numbers of stable ones
                theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
                theseCorAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
                theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        
                if unitCount==1
                [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseCorProTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
                [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,theseCorAntiTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow);
                else
                     [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseCorProTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow,'binSize',0.02);
                [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,theseCorAntiTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow,'binSize',0.02);
              
                end
                
                if unitCount==1
                allNeuronsSdf.Pro{neuronNum} =  muSdf.Pro;
                allNeuronsSdf.Anti{neuronNum} =  muSdf.Pro;
              end
                [hF] = createConfLimSdf(alignedTrials,[1 1]);
              line([0 0],[0 200],'color','k','parent',gca,'yliminclude','off')
              export_fig(hF,'-png',[resultsDir filesep neuronList(neuronNum).neuronName char(saccadeSdfPngNameEnds{unitCount})])
        delete(hF)
        if unitCount==1
                 [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseCorProTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow);
                [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,theseCorAntiTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow);
             
        else
                   [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseCorProTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow,'binSize',0.02);
                [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,theseCorAntiTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow,'binSize',0.02);
            
            
        end
                [hF] = createConfLimSdf(alignedTrials,[1 1]);
                line([0 0],[0 200],'color','k','parent',gca,'yliminclude','off')
              export_fig(hF,'-png',[resultsDir filesep neuronList(neuronNum).neuronName char(instructSdfPngNameEnds{unitCount})])
        delete(hF)
              unitCount =  unitCount+1;
    end
end

%now get all of the burst neurons
burstNeurons = find(strcmpi(allNeuronClassification,'burst'));
pauseNeurons = find(strcmpi(allNeuronClassification,'pause'));
nsNeurons = find(strcmpi(allNeuronClassification,'not sig'));
allSdfs.burst = {allNeuronsSdf.Pro{burstNeurons}};
allSdfs.pause = {allNeuronsSdf.Pro{pauseNeurons}};
allSdfs.burst = [allSdfs.burst{:}]
allSdfs.pause = [allSdfs.pause{:}]
%saccade starts at 4000, so baseline is 3800:3900


burstVerm = intersect(burstNeurons,vermNeurons);
pauseVerm = intersect(pauseNeurons,vermNeurons);
nsVerm = intersect(nsNeurons,vermNeurons);

timeVec = linspace(-4,4,8001);
hF = figure;
allSdfs.burst = {allNeuronsSdf.Pro{burstVerm}};
allSdfs.pause = {allNeuronsSdf.Pro{pauseVerm}};
allSdfs.ns= {allNeuronsSdf.Pro{nsVerm}};
allSdfs.burst = [allSdfs.burst{:}];
allSdfs.pause = [allSdfs.pause{:}];
allSdfs.ns= [allSdfs.ns{:}];

baselineRates.burst = mean(allSdfs.burst(3800:3900,:));
baselineRates.pause = mean(allSdfs.pause(3800:3900,:));
baselineRates.ns = mean(allSdfs.ns(3800:3900,:));
if ~isempty(allSdfs.burst)
line(timeVec,bsxfun(@minus,allSdfs.burst,baselineRates.burst),'color','b')
end
if ~isempty(allSdfs.pause)
line(timeVec,bsxfun(@minus,allSdfs.pause,baselineRates.pause),'color','r')
end
if ~isempty(allSdfs.ns)
line(timeVec,bsxfun(@minus,allSdfs.ns,baselineRates.ns),'color','k')
end
title('Vermis Neurons category')
xlim([-1 1])
xlabel('time (s)')
ylabel('rate difference from baseline (Hz)')
export_fig(hF,'-pdf','-append',[resultsDir filesep 'neuronClass.pdf'])





burstLat = intersect(burstNeurons,latNeurons);
pauseLat = intersect(pauseNeurons,latNeurons);
nsLat = intersect(nsNeurons,latNeurons);

hF = figure;
allSdfs.burst = {allNeuronsSdf.Pro{burstLat}};
allSdfs.pause = {allNeuronsSdf.Pro{pauseLat}};
allSdfs.ns= {allNeuronsSdf.Pro{nsLat}};
allSdfs.burst = [allSdfs.burst{:}];
allSdfs.pause = [allSdfs.pause{:}];
allSdfs.ns= [allSdfs.ns{:}];
baselineRates.burst = mean(allSdfs.burst(3800:3900,:));
baselineRates.pause = mean(allSdfs.pause(3800:3900,:));
baselineRates.ns = mean(allSdfs.ns(3800:3900,:));
if ~isempty(allSdfs.burst)
line(timeVec,bsxfun(@minus,allSdfs.burst,baselineRates.burst),'color','b')
end
if ~isempty(allSdfs.pause)
line(timeVec,bsxfun(@minus,allSdfs.pause,baselineRates.pause),'color','r')
end
if ~isempty(allSdfs.ns)
line(timeVec,bsxfun(@minus,allSdfs.ns,baselineRates.ns),'color','k')
end
title('Lateral Neurons category')
xlim([-1 1])
xlabel('time (s)')
ylabel('rate difference from baseline (Hz)')
 export_fig(hF,'-pdf','-append',[resultsDir filesep 'neuronClass.pdf'])
      
%% FINDING THE TIME OF MAXIMAL CS RESPONSE

stInit= cell(max(recordingsToAnalyse),1);
peakRelativeToSaccade = struct('Pro',stInit,'Anti',stInit,'All',stInit);
peakRelativeToFix = struct('Pro',stInit,'Anti',stInit,'All',stInit);

for neuronNum = recordingsToAnalyse
    
     thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
    unitCount =1;
    for unitNum  = thisNeuronUnits(2)
  
        %get the trials and the numbers of stable ones
                theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
                theseCorAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
                theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
        
                
                [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseCorProTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow,'binSize',0.02);
                [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,theseCorAntiTrialNums,'saccadeTime','unitNum',unitNum,'timeWindow',timeWindow,'binSize',0.02);
             %smooth to 30ms
             %where to look in mean trace
             timetoSearch = [3500:4500]; %Ssamples, %TODO code as proper option
               % smPro = smooth(muSdf.Pro,90); %TODO switch to gausiian
                smPro = muSdf.Pro;
                [val ind] = max(smPro(timetoSearch));
                peakRelativeToSaccade(neuronNum).Pro = (ind+timetoSearch(1)-timeWindow(1)*1000)/1000;
                
               % smAnti = smooth(muSdf.Anti,90); %TODO switch to gausiian
                smAnti = muSdf.Anti;
                [val ind] = max(smAnti(timetoSearch));
                peakRelativeToSaccade(neuronNum).Anti = (ind+timetoSearch(1)-timeWindow(1)*1000)/1000;
                   
                smAll = nanmean([muSdf.Pro muSdf.Anti],2);
                [val ind] = max(smAll(timetoSearch));
                peakRelativeToSaccade(neuronNum).All = (ind+timetoSearch(1)-timeWindow(1)*1000)/1000;
                
                [hF] = createConfLimSdf(alignedTrials,[1 1]);
              line([0 0],[0 200],'color','k','parent',gca,'yliminclude','off')
              line([peakRelativeToSaccade(neuronNum).Pro peakRelativeToSaccade(neuronNum).Pro],[0 200],'color','g','parent',gca,'yliminclude','off','linewidth',2);
              %now also mark with vertical lines the found peaks
              line([peakRelativeToFix(neuronNum).All peakRelativeToSaccade(neuronNum).All],[0 200],'color','m','parent',gca,'yliminclude','off','linewidth',2);
              line([peakRelativeToSaccade(neuronNum).Anti peakRelativeToSaccade(neuronNum).Anti],[0 200],'color','r','parent',gca,'yliminclude','off','linewidth',2);
             % export_fig(hF,'-png',[resultsDir filesep neuronList(neuronNum).neuronName '-csTimingTest.png'])
              
              title([neuronList(neuronNum).neuronName ' Saccade'])
             export_fig(hF,'-pdf','-append',[resultsDir filesep csTimingPdfName])
        delete(hF)
        
         %now do the same fo relative to isntruction
        
                 [alignedTrials.Pro, muSdf.Pro] = alignAndExtractRates(theseStableTrials,theseCorProTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow,'binSize',0.02);
                [alignedTrials.Anti, muSdf.Anti] = alignAndExtractRates(theseStableTrials,theseCorAntiTrialNums,'fixation','unitNum',unitNum,'timeWindow',timeWindow,'binSize',0.02);
              
                 %smooth to 30ms
                %smPro = smooth(muSdf.Pro,30); %TODO switch to gausiian
                smPro = muSdf.Pro;
                [val ind] = max(smPro(timetoSearch));
                peakRelativeToFix(neuronNum).Pro = (ind+timetoSearch(1)-timeWindow(1)*1000)/1000;
                
               % smAnti = smooth(muSdf.Anti,30); %TODO switch to gausiian
               smAnti =  muSdf.Anti;
               [val ind] = max(smAnti(timetoSearch));
                peakRelativeToFix(neuronNum).Anti = (ind+timetoSearch(1)-timeWindow(1)*1000)/1000;
                
               % smAll = smooth(sum([muSdf.Pro muSdf.Anti],2),30)
                smAll = nanmean([muSdf.Pro muSdf.Anti],2);
                [val ind] = max(smAll(timetoSearch));
                peakRelativeToFix(neuronNum).All = (ind+timetoSearch(1)-timeWindow(1)*1000)/1000;
                
                [hF] = createConfLimSdf(alignedTrials,[1 1]);
                line([0 0],[0 200],'color','k','parent',gca,'yliminclude','off')
              
                line([peakRelativeToFix(neuronNum).Pro peakRelativeToFix(neuronNum).Pro],[0 200],'color','g','parent',gca,'yliminclude','off','linewidth',2);
              %now also mark with vertical lines the found peaks
          line([peakRelativeToFix(neuronNum).Anti peakRelativeToFix(neuronNum).Anti],[0 200],'color','r','parent',gca,'yliminclude','off','linewidth',2);
          line([peakRelativeToFix(neuronNum).All peakRelativeToFix(neuronNum).All],[0 200],'color','m','parent',gca,'yliminclude','off','linewidth',2);
               title([neuronList(neuronNum).neuronName ' Instruction'])
          export_fig(hF,'-pdf','-append',[resultsDir filesep csTimingPdfName])
               % export_fig(hF,'-png',[resultsDir filesep neuronList(neuronNum).neuronName '-csTimingTest.png'])
        delete(hF)
              unitCount =  unitCount+1;
    end
end

hF = figure;
hA(1) = subplot(3,1,1,'parent',hF,'nextplot','add','xlim',[-0.5 0.5]); %vermis
hA(2) = subplot(3,1,2,'parent',hF,'nextplot','add','xlim',[-0.5 0.5]); %lateral
hA(3) = subplot(3,1,3,'parent',hF,'nextplot','add','xlim',[-0.5 0.5]); %all
%now do figure with 3 subplots verm, lat and all, with different colours
%for pro and anti. 
[counts,centers] = hist([peakRelativeToSaccade(vermNeurons).All],50);
bar(centers,counts,'parent',hA(1),'m')
[counts,centers] = hist([peakRelativeToSaccade(vermNeurons).Pro],50);
bar(centers,counts,'parent',hA(1),'g')

[counts,centers] = hist([peakRelativeToSaccade(vermNeurons).Anti],50);
bar(centers,counts,'parent',hA(1),'r')


set(get(hA(1),'title'),'string','Vermis Neurons')


[counts,centers] = hist([peakRelativeToSaccade(latNeurons).All],50);
bar(centers,counts,'parent',hA(2),'m')
[counts,centers] = hist([peakRelativeToSaccade(latNeurons).Pro],50);
bar(centers,counts,'parent',hA(2),'g')

[counts,centers] = hist([peakRelativeToSaccade(latNeurons).Anti],50);
bar(centers,counts,'parent',hA(2),'r')


set(get(hA(2),'title'),'string','Lateral Neurons')


[counts,centers] = hist([peakRelativeToSaccade.All],50);
bar(centers,counts,'parent',hA(3),'m')
[counts,centers] = hist([peakRelativeToSaccade.Pro],50);
bar(centers,counts,'parent',hA(3),'g')

[counts,centers] = hist([peakRelativeToSaccade.Anti],50);
bar(centers,counts,'parent',hA(3),'r')

set(get(hA(3),'title'),'string','All Neurons')
suptitle('Relative to Saccade')
 export_fig(hF,'-pdf','-append',[resultsDir filesep csTimingPdfName])
               
        delete(hF)



hF = figure;
hA(1) = subplot(3,1,1,'parent',hF,'nextplot','add','xlim',[-0.5 0.5]); %vermis
hA(2) = subplot(3,1,2,'parent',hF,'nextplot','add','xlim',[-0.5 0.5]); %lateral
hA(3) = subplot(3,1,3,'parent',hF,'nextplot','add','xlim',[-0.5 0.5]); %all
%now do figure with 3 subplots verm, lat and all, with different colours
%for pro and anti. 
[counts,centers] = hist([peakRelativeToFix(vermNeurons).All],50);
bar(centers,counts,'parent',hA(1),'m')
[counts,centers] = hist([peakRelativeToFix(vermNeurons).Pro],50);
bar(centers,counts,'parent',hA(1),'g')

[counts,centers] = hist([peakRelativeToFix(vermNeurons).Anti],50);
bar(centers,counts,'parent',hA(1),'r')


set(get(hA(1),'title'),'string','Vermis Neurons')


[counts,centers] = hist([peakRelativeToFix(latNeurons).All],50);
bar(centers,counts,'parent',hA(2),'m')
[counts,centers] = hist([peakRelativeToFix(latNeurons).Pro],50);
bar(centers,counts,'parent',hA(2),'g')

[counts,centers] = hist([peakRelativeToFix(latNeurons).Anti],50);
bar(centers,counts,'parent',hA(2),'r')


set(get(hA(2),'title'),'string','Lateral Neurons')

[counts,centers] = hist([peakRelativeToFix.All],50);
bar(centers,counts,'parent',hA(3),'m')
[counts,centers] = hist([peakRelativeToFix.Pro],50);
bar(centers,counts,'parent',hA(3),'g')

[counts,centers] = hist([peakRelativeToFix.Anti],50);
bar(centers,counts,'parent',hA(3),'r')


set(get(hA(3),'title'),'string','All Neurons')
suptitle('Relative to Fixation')

 export_fig(hF,'-pdf','-append',[resultsDir filesep csTimingPdfName])
              
        delete(hF)
        
        
        %TODO can also try a simple sum across all smoothed sfd's 
%% PCA of sdf responses and then label by loation to see if we can see clear clusters in pcatest.m


%% count cs in mutliple bins.
%this should be similar to the cs count cell above but should now work over
%4 windows.
%1 - instruction - 0:150ms post fixation
%2 - presaccade - -200:0ms saccade onset
%3 - saccade - during saccade, need to find end points
%4 - post-saccade - 0:200ms post end of saccade, need to find endpoints
%trialConditionCodes is in subplot order

if doSpikesInBins
inAngleOrder = [6 3 2 1 4 7 8 9]; %order of trialConditionCodes(:,dirLabelNum) in terms of starting at angle 0 (horizontal right) and counting coutner clockwise

allNeuronDirectionCountResults = cell(max(recordingsToAnalyse),1);
allNeuronDirectionCountResultsProCsAlign = allNeuronDirectionCountResults;
for neuronNum = recordingsToAnalyse
    thisNeuronFigFileName = [resultsDir  filesep recList{1,neuronNum} '-binCountsAndAlignToMax092.pdf'];
    
    %([pross,anti ss, pro cs, anti cs],direction)
    directionCountInit = nan(4,8);
    directionCountResults = struct('instructionBin',directionCountInit,'preSaccBin',directionCountInit,'periSaccBin',directionCountInit,'postSaccBin',directionCountInit);
    theseCorrectTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corSacTrials;
    %      theseCorProTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corProTrials;
    %                 theseCorAntiTrialNums = wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials;
    theseStableTrials = wholeNeuronResults(neuronNum).allStableTrials;
    
    thisNeuronUnits = neuronList(neuronNum).fileList{1,3};
    
    for  binName = {'instructionBin','preSaccBin','periSaccBin','postSaccBin'}
        unitCount =1;
        
        theseResults = nan(4,8);
        for unitNum  = thisNeuronUnits
            
            [countResults axHand] = countSpikesInTrialBins(theseStableTrials,theseCorrectTrialNums,timeWindow,'',unitNum)
            
            
            
            %now get the average for each of the 4 within condtionCodes (ie direction
            %and pro anti)
            
            for directionNum = 1:8
                for typeCount = 1:2%pro and anti
                    theseTrialNums = find([countResults.conditionCode] == trialConditionCodes(inAngleOrder(directionNum),typeCount));
                    if ~isempty(theseTrialNums)
                        theseCounts = countResults(theseTrialNums);
                        theseResults(typeCount+2*(unitCount-1),directionNum) = mean([theseCounts.(char(binName))]);
                    end
                    
                end
            end
            unitCount = unitCount+1;
        end
        directionCountResults.(char(binName)) = theseResults;
        
    end
    allNeuronDirectionCountResults{neuronNum} = directionCountResults;
    
    
    for  binName = {'instructionBin','preSaccBin','periSaccBin','postSaccBin'}
        %now show all as they are on two separate overlaying axis
        hF = figure;
        hA(1) = axes('parent',hF,'units','normalized','position',[0.1 0.05 0.8 0.35]); %ss
        hA(2) = axes('parent',hF,'units','normalized','position',[0.1 0.5 0.8 0.35]); %cs
        linkaxes(hA,'x')
        set(hA(1),'xlim',[0 9])
        angAxLabs = {'3','4.5','6','7.5','9','10.5','12','1.5'};
        angAxTicks = [1:8];
        theseResults = allNeuronDirectionCountResults{neuronNum}.(char(binName));
        line(angAxTicks,theseResults(1,:),'parent',hA(1),'color','g','linestyle','-','linewidth',2)
        line(angAxTicks,theseResults(2,:),'parent',hA(1),'color','r','linestyle','-','linewidth',2)
        line(angAxTicks,theseResults(3,:),'parent',hA(2),'color','g','linestyle','-','linewidth',2)
        line(angAxTicks,theseResults(4,:),'parent',hA(2),'color','r','linestyle','-','linewidth',2)
        set(get(hA(1),'title'),'string','SS')
        set(get(hA(2),'title'),'string','CS')
        set(hA(2),'xtick',[])
        set(hA(1),'xtick',angAxTicks,'xticklabel',angAxLabs)
        suptitle([neuronList(neuronNum).neuronName '-' char(binName) '-raw'])
        export_fig(thisNeuronFigFileName, '-pdf','-append', hF);
        delete(hF)
        
        
        %now do aplot which aligns all to the highest CS direction in pro
        %trials
        alignmentTarget = 4; %colun to put peak in
        alignmentRow = 3; %proCS %TODO could loop here
        
        [peakVal peakDir] = max(theseResults(alignmentRow,:));
        if numel(peakDir)~=1
            a=1;
        end
        %so we now want to do a circshift to move peakDir to position 4 and
        %all results as well
        alignedResults = circshift(theseResults,[0 alignmentTarget-peakDir]); %0 shift in fist dim as much as needed in second
        
        hF = figure;
        hA(1) = axes('parent',hF,'units','normalized','position',[0.1 0.05 0.8 0.35]); %ss
        hA(2) = axes('parent',hF,'units','normalized','position',[0.1 0.5 0.8 0.35]); %cs
        linkaxes(hA,'x')
        set(hA(1),'xlim',[0 9])
        angAxLabs = {'3','4.5','6','7.5','9','10.5','12','1.5'};
        angAxTicks = [1:8];
        theseResults = allNeuronDirectionCountResults{neuronNum}.(char(binName));
        line(angAxTicks,alignedResults(1,:),'parent',hA(1),'color','g','linestyle','-','linewidth',2)
        line(angAxTicks,alignedResults(2,:),'parent',hA(1),'color','r','linestyle','-','linewidth',2)
        line(angAxTicks, alignedResults(3,:),'parent',hA(2),'color','g','linestyle','-','linewidth',2)
        line(angAxTicks, alignedResults(4,:),'parent',hA(2),'color','r','linestyle','-','linewidth',2)
        set(get(hA(1),'title'),'string','SS')
        set(get(hA(2),'title'),'string','CS')
        set(hA(2),'xtick',[])
        set(hA(1),'xtick',angAxTicks,'xticklabel',angAxLabs)
        suptitle([neuronList(neuronNum).neuronName '-' char(binName) '-alignedToProCs'])
        export_fig(thisNeuronFigFileName, '-pdf','-append', hF);
        delete(hF)
        
        
        allNeuronDirectionCountResultsProCsAlign{neuronNum}.(char(binName)) = alignedResults;
    end
end


%now show all nonaligned on top of each other
for  binName = {'instructionBin','preSaccBin','periSaccBin','postSaccBin'}
    %first all neurons
    alignedResults = [allNeuronDirectionCountResultsProCsAlign{recordingsToAnalyse}];
    
    
    thisBinAlignedResults = cat(3,alignedResults.(char(binName)));
    
    proSS = squeeze(thisBinAlignedResults(1,:,:)); %now each column is  neuron and the rows are th 8 directions
    antiSS = squeeze(thisBinAlignedResults(2,:,:));
    proCS = squeeze(thisBinAlignedResults(3,:,:));
    antiCS = squeeze(thisBinAlignedResults(4,:,:));
    
    muProSS = mean(proSS,2);
    semProSS = std(proSS,1,2)./sqrt(length(recordingsToAnalyse));
    uppProSS = muProSS+semProSS;
    lowProSS = muProSS-semProSS;
    
    muAntiSS = mean(antiSS,2);
    semAntiSS = std(antiSS,1,2)./sqrt(length(recordingsToAnalyse));
    uppAntiSS = muAntiSS+semAntiSS;
    lowAntiSS = muAntiSS-semAntiSS;
    
    muProCS = mean(proCS,2);
    semProCS = std(proCS,1,2)./sqrt(length(recordingsToAnalyse));
    uppProCS = muProCS+semProCS;
    lowProCS = muProCS-semProCS;
    
    muAntiCS = mean(antiCS,2);
    semAntiCS = std(antiCS,1,2)./sqrt(length(recordingsToAnalyse));
    uppAntiCS = muAntiCS+semAntiCS;
    lowAntiCS = muAntiCS-semAntiCS;
    
    hF = figure;
    hA(1) = axes('parent',hF,'units','normalized','position',[0.1 0.05 0.8 0.35],'ylim',[50 100]); %ss
    hA(2) = axes('parent',hF,'units','normalized','position',[0.1 0.5 0.8 0.35],'ylim',[0 3]); %cs
    linkaxes(hA,'x')
    set(hA(1),'xlim',[0 9])
    angAxLabs = {'3','4.5','6','7.5','9','10.5','12','1.5'};
    angAxTicks = [1:8];
    
    line(angAxTicks,muProSS,'parent',hA(1),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiSS,'parent',hA(1),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muProCS,'parent',hA(2),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiCS,'parent',hA(2),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    %TODO now need to add patches for sems instead of dotted lines
    
    set(get(hA(1),'title'),'string','SS')
    set(get(hA(2),'title'),'string','CS')
    set(hA(2),'xtick',[])
    set(hA(1),'xtick',angAxTicks,'xticklabel',angAxLabs)
    suptitle(['All Neurons alignedToMaxProCs-' char(binName)])
    export_fig( allNeuronPdfName, '-pdf','-append', hF);
    delete(hF)
    
    %now do exactly the same to the raw results
    rawResults = [allNeuronDirectionCountResults{recordingsToAnalyse}];
    
    
    
    thisBinRawResults = cat(3,rawResults.(char(binName)));
    
    proSS = squeeze(thisBinRawResults(1,:,:)); %now each column is  neuron and the rows are th 8 directions
    antiSS = squeeze(thisBinRawResults(2,:,:));
    proCS = squeeze(thisBinRawResults(3,:,:));
    antiCS = squeeze(thisBinRawResults(4,:,:));
    
    muProSS = mean(proSS,2);
    semProSS = std(proSS,1,2)./sqrt(length(recordingsToAnalyse));
    uppProSS = muProSS+semProSS;
    lowProSS = muProSS-semProSS;
    
    muAntiSS = mean(antiSS,2);
    semAntiSS = std(antiSS,1,2)./sqrt(length(recordingsToAnalyse));
    uppAntiSS = muAntiSS+semAntiSS;
    lowAntiSS = muAntiSS-semAntiSS;
    
    muProCS = mean(proCS,2);
    semProCS = std(proCS,1,2)./sqrt(length(recordingsToAnalyse));
    uppProCS = muProCS+semProCS;
    lowProCS = muProCS-semProCS;
    
    muAntiCS = mean(antiCS,2);
    semAntiCS = std(antiCS,1,2)./sqrt(length(recordingsToAnalyse));
    uppAntiCS = muAntiCS+semAntiCS;
    lowAntiCS = muAntiCS-semAntiCS;
    
    hF = figure;
    hA(1) = axes('parent',hF,'units','normalized','position',[0.1 0.05 0.8 0.35],'ylim',[50 100]); %ss
    hA(2) = axes('parent',hF,'units','normalized','position',[0.1 0.5 0.8 0.35],'ylim',[0 3]); %cs
    linkaxes(hA,'x')
    set(hA(1),'xlim',[0 9])
    angAxLabs = {'3','4.5','6','7.5','9','10.5','12','1.5'};
    angAxTicks = [1:8];
    
    line(angAxTicks,muProSS,'parent',hA(1),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiSS,'parent',hA(1),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muProCS,'parent',hA(2),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiCS,'parent',hA(2),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    
    set(get(hA(1),'title'),'string','SS')
    set(get(hA(2),'title'),'string','CS')
    set(hA(2),'xtick',[])
    set(hA(1),'xtick',angAxTicks,'xticklabel',angAxLabs)
    suptitle(['All Neurons no alignment-' char(binName)])
    export_fig( allNeuronPdfName, '-pdf','-append', hF);
    delete(hF)
    
    
    %now vermis neurons
    
    alignedResults = [allNeuronDirectionCountResultsProCsAlign{vermNeurons}];
    
    
    thisBinAlignedResults = cat(3,alignedResults.(char(binName)));
    
    proSS = squeeze(thisBinAlignedResults(1,:,:)); %now each column is  neuron and the rows are th 8 directions
    antiSS = squeeze(thisBinAlignedResults(2,:,:));
    proCS = squeeze(thisBinAlignedResults(3,:,:));
    antiCS = squeeze(thisBinAlignedResults(4,:,:));
    
    muProSS = mean(proSS,2);
    semProSS = std(proSS,1,2)./sqrt(length(vermNeurons));
    uppProSS = muProSS+semProSS;
    lowProSS = muProSS-semProSS;
    
    muAntiSS = mean(antiSS,2);
    semAntiSS = std(antiSS,1,2)./sqrt(length(vermNeurons));
    uppAntiSS = muAntiSS+semAntiSS;
    lowAntiSS = muAntiSS-semAntiSS;
    
    muProCS = mean(proCS,2);
    semProCS = std(proCS,1,2)./sqrt(length(vermNeurons));
    uppProCS = muProCS+semProCS;
    lowProCS = muProCS-semProCS;
    
    muAntiCS = mean(antiCS,2);
    semAntiCS = std(antiCS,1,2)./sqrt(length(vermNeurons));
    uppAntiCS = muAntiCS+semAntiCS;
    lowAntiCS = muAntiCS-semAntiCS;
    
    hF = figure;
    hA(1) = axes('parent',hF,'units','normalized','position',[0.1 0.05 0.8 0.35],'ylim',[50 100]); %ss
    hA(2) = axes('parent',hF,'units','normalized','position',[0.1 0.5 0.8 0.35],'ylim',[0 3]); %cs
    linkaxes(hA,'x')
    set(hA(1),'xlim',[0 9])
    angAxLabs = {'3','4.5','6','7.5','9','10.5','12','1.5'};
    angAxTicks = [1:8];
    
    line(angAxTicks,muProSS,'parent',hA(1),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiSS,'parent',hA(1),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muProCS,'parent',hA(2),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiCS,'parent',hA(2),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    %TODO now need to add patches for sems instead of dotted lines
    
    set(get(hA(1),'title'),'string','SS')
    set(get(hA(2),'title'),'string','CS')
    set(hA(2),'xtick',[])
    set(hA(1),'xtick',angAxTicks,'xticklabel',angAxLabs)
    suptitle(['Vermis Neurons alignedToMaxProCs-' char(binName)])
    export_fig( allNeuronPdfName, '-pdf','-append', hF);
    delete(hF)
    
    %now do exactly the same to the raw results
    rawResults = [allNeuronDirectionCountResults{vermNeurons}];
    
    
    
    thisBinRawResults = cat(3,rawResults.(char(binName)));
    
    proSS = squeeze(thisBinRawResults(1,:,:)); %now each column is  neuron and the rows are th 8 directions
    antiSS = squeeze(thisBinRawResults(2,:,:));
    proCS = squeeze(thisBinRawResults(3,:,:));
    antiCS = squeeze(thisBinRawResults(4,:,:));
    
    muProSS = mean(proSS,2);
    semProSS = std(proSS,1,2)./sqrt(length(vermNeurons));
    uppProSS = muProSS+semProSS;
    lowProSS = muProSS-semProSS;
    
    muAntiSS = mean(antiSS,2);
    semAntiSS = std(antiSS,1,2)./sqrt(length(vermNeurons));
    uppAntiSS = muAntiSS+semAntiSS;
    lowAntiSS = muAntiSS-semAntiSS;
    
    muProCS = mean(proCS,2);
    semProCS = std(proCS,1,2)./sqrt(length(vermNeurons));
    uppProCS = muProCS+semProCS;
    lowProCS = muProCS-semProCS;
    
    muAntiCS = mean(antiCS,2);
    semAntiCS = std(antiCS,1,2)./sqrt(length(vermNeurons));
    uppAntiCS = muAntiCS+semAntiCS;
    lowAntiCS = muAntiCS-semAntiCS;
    
    hF = figure;
    hA(1) = axes('parent',hF,'units','normalized','position',[0.1 0.05 0.8 0.35],'ylim',[50 100]); %ss
    hA(2) = axes('parent',hF,'units','normalized','position',[0.1 0.5 0.8 0.35],'ylim',[0 3]); %cs
    linkaxes(hA,'x')
    set(hA(1),'xlim',[0 9])
    angAxLabs = {'3','4.5','6','7.5','9','10.5','12','1.5'};
    angAxTicks = [1:8];
    
    line(angAxTicks,muProSS,'parent',hA(1),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiSS,'parent',hA(1),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muProCS,'parent',hA(2),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiCS,'parent',hA(2),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    
    set(get(hA(1),'title'),'string','SS')
    set(get(hA(2),'title'),'string','CS')
    set(hA(2),'xtick',[])
    set(hA(1),'xtick',angAxTicks,'xticklabel',angAxLabs)
    suptitle(['Vermis Neurons no alignment-' char(binName)])
    export_fig( allNeuronPdfName, '-pdf','-append', hF);
    delete(hF)
    
    
    %now lateral neurons
    alignedResults = [allNeuronDirectionCountResultsProCsAlign{latNeurons}];
    
    
    thisBinAlignedResults = cat(3,alignedResults.(char(binName)));
    
    proSS = squeeze(thisBinAlignedResults(1,:,:)); %now each column is  neuron and the rows are th 8 directions
    antiSS = squeeze(thisBinAlignedResults(2,:,:));
    proCS = squeeze(thisBinAlignedResults(3,:,:));
    antiCS = squeeze(thisBinAlignedResults(4,:,:));
    
    muProSS = mean(proSS,2);
    semProSS = std(proSS,1,2)./sqrt(length(latNeurons));
    uppProSS = muProSS+semProSS;
    lowProSS = muProSS-semProSS;
    
    muAntiSS = mean(antiSS,2);
    semAntiSS = std(antiSS,1,2)./sqrt(length(latNeurons));
    uppAntiSS = muAntiSS+semAntiSS;
    lowAntiSS = muAntiSS-semAntiSS;
    
    muProCS = mean(proCS,2);
    semProCS = std(proCS,1,2)./sqrt(length(latNeurons));
    uppProCS = muProCS+semProCS;
    lowProCS = muProCS-semProCS;
    
    muAntiCS = mean(antiCS,2);
    semAntiCS = std(antiCS,1,2)./sqrt(length(latNeurons));
    uppAntiCS = muAntiCS+semAntiCS;
    lowAntiCS = muAntiCS-semAntiCS;
    
    hF = figure;
    hA(1) = axes('parent',hF,'units','normalized','position',[0.1 0.05 0.8 0.35],'ylim',[50 100]); %ss
    hA(2) = axes('parent',hF,'units','normalized','position',[0.1 0.5 0.8 0.35],'ylim',[0 3]); %cs
    linkaxes(hA,'x')
    set(hA(1),'xlim',[0 9])
    angAxLabs = {'3','4.5','6','7.5','9','10.5','12','1.5'};
    angAxTicks = [1:8];
    
    line(angAxTicks,muProSS,'parent',hA(1),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiSS,'parent',hA(1),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muProCS,'parent',hA(2),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiCS,'parent',hA(2),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    %TODO now need to add patches for sems instead of dotted lines
    
    set(get(hA(1),'title'),'string','SS')
    set(get(hA(2),'title'),'string','CS')
    set(hA(2),'xtick',[])
    set(hA(1),'xtick',angAxTicks,'xticklabel',angAxLabs)
    suptitle(['Lateral Neurons alignedToMaxProCs-' char(binName)])
    export_fig( allNeuronPdfName, '-pdf','-append', hF);
    delete(hF)
    
    %now do exactly the same to the raw results
    rawResults = [allNeuronDirectionCountResults{latNeurons}];
    
    
    
    thisBinRawResults = cat(3,rawResults.(char(binName)));
    
    proSS = squeeze(thisBinRawResults(1,:,:)); %now each column is  neuron and the rows are th 8 directions
    antiSS = squeeze(thisBinRawResults(2,:,:));
    proCS = squeeze(thisBinRawResults(3,:,:));
    antiCS = squeeze(thisBinRawResults(4,:,:));
    
    muProSS = mean(proSS,2);
    semProSS = std(proSS,1,2)./sqrt(length(latNeurons));
    uppProSS = muProSS+semProSS;
    lowProSS = muProSS-semProSS;
    
    muAntiSS = mean(antiSS,2);
    semAntiSS = std(antiSS,1,2)./sqrt(length(latNeurons));
    uppAntiSS = muAntiSS+semAntiSS;
    lowAntiSS = muAntiSS-semAntiSS;
    
    muProCS = mean(proCS,2);
    semProCS = std(proCS,1,2)./sqrt(length(latNeurons));
    uppProCS = muProCS+semProCS;
    lowProCS = muProCS-semProCS;
    
    muAntiCS = mean(antiCS,2);
    semAntiCS = std(antiCS,1,2)./sqrt(length(latNeurons));
    uppAntiCS = muAntiCS+semAntiCS;
    lowAntiCS = muAntiCS-semAntiCS;
    
    hF = figure;
    hA(1) = axes('parent',hF,'units','normalized','position',[0.1 0.05 0.8 0.35],'ylim',[50 100]); %ss
    hA(2) = axes('parent',hF,'units','normalized','position',[0.1 0.5 0.8 0.35],'ylim',[0 3]); %cs
    linkaxes(hA,'x')
    set(hA(1),'xlim',[0 9])
    angAxLabs = {'3','4.5','6','7.5','9','10.5','12','1.5'};
    angAxTicks = [1:8];
    
    line(angAxTicks,muProSS,'parent',hA(1),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProSS,'parent',hA(1),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiSS,'parent',hA(1),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiSS,'parent',hA(1),'color','r','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muProCS,'parent',hA(2),'color','g','linestyle','-','linewidth',2)
    line(angAxTicks,uppProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    line(angAxTicks,lowProCS,'parent',hA(2),'color','g','linestyle','-.','linewidth',2)
    
    line(angAxTicks,muAntiCS,'parent',hA(2),'color','r','linestyle','-','linewidth',2)
    line(angAxTicks,uppAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    line(angAxTicks,lowAntiCS,'parent',hA(2),'color','r','linestyle','-.','linewidth',2)
    
    set(get(hA(1),'title'),'string','SS')
    set(get(hA(2),'title'),'string','CS')
    set(hA(2),'xtick',[])
    set(hA(1),'xtick',angAxTicks,'xticklabel',angAxLabs)
    suptitle(['Lateral Neurons no alignment-' char(binName)])
    export_fig( allNeuronPdfName, '-pdf','-append', hF);
    delete(hF)
end
end
%%
%now get the percentage correct trials and correct trials with saccades for
%each neuron.
stInit =cell(max(recordingsToAnalyse),1);
trialCounts = struct('proNumTrials',stInit,'proNumCorTrials',stInit,...
    'antiNumTrials',stInit,'antiNumCorTrials',stInit);
   
    %total, correct, correctsacc, pro, pro-correct,pro
for neuronNum = recordingsToAnalyse
    trialCounts(neuronNum).proNumTrials = numel(wholeNeuronResults(neuronNum).selectedTrials.proTrials)   ;
    trialCounts(neuronNum).proNumCorTrials =numel(wholeNeuronResults(neuronNum).selectedTrials.corProTrials);
    trialCounts(neuronNum).antiNumCorTrials = numel(wholeNeuronResults(neuronNum).selectedTrials.corAntiTrials);
    trialCounts(neuronNum).antiNumTrials = numel(wholeNeuronResults(neuronNum).selectedTrials.antiTrials) ;

   

end

 correctPercentages = cell(length(recordingsToAnalyse),3); %neuron name, %corPro, %corAnti 

proCorPer = [trialCounts.proNumCorTrials]./[trialCounts.proNumTrials]*100;
antiCorPer = [trialCounts.antiNumCorTrials]./[trialCounts.antiNumTrials]*100;

correctPercentages(:,1) = {neuronList(recordingsToAnalyse).neuronName};
correctPercentages(:,2) = num2cell(proCorPer);
correctPercentages(:,3) = num2cell(antiCorPer);

xlswrite([resultsDir filesep 'neuronCorrectTrialPercents.xlsx'],correctPercentages)