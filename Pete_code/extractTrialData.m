function [trialData] = extractTrialData(thisFileInfo,dataStartTime,bitNum,secondBitNum,processDir,neuronList)


fileNameBase = thisFileInfo.fileName(1:end-4);

%load stimulus list file %TODO should warn if it doesnt exist and
%launch the stimDir2List function somehow
temp = strfind(fileNameBase, '_');
if isempty(temp)
 temp =   strfind(fileNameBase, '-');
end

try
load([processDir filesep fileNameBase(1:temp(end)-1) '-stimList']);
catch em
    
     
    load(neuronList(thisFileInfo.neuronNum).stimListFile)
end


%load ports (event channels) from AO .mat file (output from their
%converter)
%the new version labels things CInPort
%  try
 w = whos('-file',thisFileInfo.fileName);
  if ~isempty(find(strcmpi({w.name},'CPort_005')))
     load(thisFileInfo.fileName,'CPort_005','CPort_003_KHz','CPort_003');
     fs = CPort_003_KHz*1000;
     inFlag = 0;
     portNameString = 'CPort';
     bitPortNum = 5;
     wordPortNum = 3;
  elseif ~isempty(find(strcmpi({w.name},'CPort_004')))
     load(thisFileInfo.fileName,'CPort_005','CPort_003_KHz','CPort_003');
     fs = CPort_003_KHz*1000;
     inFlag = 0;
     portNameString = 'CPort';
     bitPortNum = 3;
     wordPortNum = 4;
 elseif ~isempty(find(strcmpi({w.name},'CInPort_005')))
     load(thisFileInfo.fileName,'CInPort_005','CInPort_003_KHz','CInPort_003');
     fs = CInPort_003_KHz*1000;
     inFlag = 1; %to indicate if this needs to be inserted later in code
     portNameString = 'CInPort';
     bitPortNum = 5;
     wordPortNum = 3;
  
 
 else
     load(thisFileInfo.fileName,'CInPort_001','CInPort_001_KHz','CInPort_002');
     fs = CInPort_001_KHz*1000;
     inFlag = 1; %to indicate if this needs to be inserted later in code
     portNameString = 'CInPort';
     bitPortNum = 2;
     wordPortNum = 1;
 end
%  catch em
%      load(thisFileInfo.fileName)
%      error_log
%      keyboard
%  end

    %!!!!!!!!!!!!FROM NOW ON ALL TIMES WITHIN FILE ARE RELATIVE TO
            %START OF FILE (spk start time) eg START OF FILE = 0s!!!!!!!!
            %TODO check if bitport starts earlier

%either read word port and bit port or reconstruct word port from
%bitport with MANUAL starting file (output of another function I
%can't remember so need switch below to launch it here if there is no
%starting file given when needed)

if ~exist([portNameString '_' sprintf('%03d',bitPortNum)],'var') %PH 171016, tried to make more general but not sure if should check for bit or word port
%if ~exist([portNameString '_005'],'var') 
    %TODO case with no port 5
   if inFlag
        bitPort = CInPort_003;
   else
      
        bitPort = CPort_003;
   end
    A = bitPort(2,:);
    timeA = bitPort(1,:);
    %get bits
    for i = 1:length(A)
        B(:,i) = bitget(A(i), 1:8);
    end
    
    %recreate wordport
    if ~isempty(thisFileInfo.stimStartFile)
        firsstimfile = thisFileInfo.stimStartFile;
        bit2times = timeA([0 diff(A)]==2);
        wordPort = [bit2times;[firsstimfile:firsstimfile+length(bit2times)-1]];
    else
        keyboard %TODO should launch finding functions
    end
else
    if inFlag
%         bitPort = CInPort_005;
%         wordPort = CInPort_003;
        %PH 171016 changed to eval to deal with Nico moving the ports to 1
        %and 2
         eval(['bitPort = double(CInPort_' sprintf('%03d',bitPortNum) ');']);
      eval(['wordPort  = double(CInPort_' sprintf('%03d',wordPortNum) ');']);
   else
      wordPort = CPort_003;
      if exist('CPort_005')
        bitPort = CPort_005;
      else 
          bitPort= CPort_004;
      end
      
   end
   
    A = bitPort(2,:);
    timeA = bitPort(1,:);
    %get bits
    for i = 1:length(A)
        B(:,i) = bitget(A(i), 1:8);
    end
    %trim wordPort and convert to stim file nums
    
    wordPort = wordPort (:,wordPort(2,:) > 32768);
    wordPort(2,:) = wordPort(2,:) - 32768;
    wordPort = wordPort(:,[1 diff(wordPort(2,:))]>0);
end

%bit port is all the actual codes sent to the bit port
%word port is the stimFileNumber and the sample (full fs) at which
%bit2 was sent to the SNR

numTrials = size(wordPort,2);
eventData = cell(9,numTrials);
%TODO: this should be a structure but needs a lot recoding
eventData(1,:) = num2cell((wordPort(1,:)-dataStartTime*fs)/fs); %first row-time in s relative to spk start
eventData(2,:) = num2cell(wordPort(2,:)); %second is stim file number
timeA = (bitPort(1,:)-dataStartTime*fs)/fs;


for trialNum = 1:numTrials
    %stimFileNum = eventData{2,1}+trialNum-1; %was this line, why?
    stimFileNum = eventData{2,trialNum};
    if stimFileNum<=numel(stimListFile) %catch misreads and ignore
    eventData{3,trialNum} = stimListFile(stimFileNum).Cor; %correct or incorrect
    eventData{4,trialNum} = stimListFile(stimFileNum).cond; %condition number (direction and anti or pro)
    
    %0 in bits signify end of trials so
    %in bit port find the next and previous occurence of 0 - watchf or ends!!!!
    [i j] = min(abs(eventData{1,trialNum}-timeA)); %find closest time in bit port
    trialStartInd = find(A(1:j)==0,1,'last')+1; %find start of trial searching backwards
    trialEndInd = find(A(j:end)==0,1,'first')+j-1; %find end of trial searching forwards
    
    eventData{9,trialNum} = timeA(trialStartInd); %start time s
    eventData{10,trialNum} = timeA(trialEndInd); %end time s
    
    eventData{5,trialNum} = [timeA(trialStartInd:trialEndInd); B(:,trialStartInd:trialEndInd)];
    
    
    
    %text(eventData{1,trialNum},40,num2str(eventData{4,trialNum})) %TODO
    %if there are bits in trial then find the occurence of the two
    %askde for bits
    if ~isempty(eventData{5,trialNum})
        %line([eventData{5,trialNum}(1,1)
        %eventData{5,trialNum}(1,end-1)],[50 50],'color','k') %TODO
        bitInd = find(eventData{5,trialNum}(bitNum+1,:)==1,1,'first');
        if ~isempty(bitInd)
            %TODO line([eventData{5,trialNum}(1,bitInd) eventData{5,trialNum}(1,bitInd)],[45 55],'color','k');
            eventData{6,trialNum} = eventData{5,trialNum}(1,bitInd);
        end
        secondBitInd = find(eventData{5,trialNum}(secondBitNum+1,:)==1,1,'first');
        if ~isempty(secondBitInd)
            %TODO line([eventData{5,trialNum}(1,secondBitInd) eventData{5,trialNum}(1,secondBitInd)],[60 70],'color','g');
            eventData{7,trialNum} = eventData{5,trialNum}(1,secondBitInd);
        end
        %eventData{8,trialNum} = eventData{5,trialNum}(1,1);
    end
    end
    
end
eventData(:,cellfun('isempty',eventData(3,:))) =[]; %trim trials missing a corect reponse as they are system errors

%this is shadow of eventData as a structure which will replace it

trialData = struct('bit2time',eventData(1,:),'stimFileNumber',eventData(2,:),...
    'correctResponse',eventData(3,:),'conditionCode',eventData(4,:),...
    'trialBits',eventData(5,:),'bit8time',eventData(6,:),...
    'bit3time',eventData(7,:) ,'trialStart',eventData(9,:),...
    'trialEnd',eventData(10,:));