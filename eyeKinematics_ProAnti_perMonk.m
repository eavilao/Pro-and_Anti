function eyeKin = eyeKinematics_ProAnti_perMonk(units)

recArea = 'lateral'; 
monkey = 'moshe'; %'moshe'; 'mickey'

%%
for cell = 1:length(units)   % PICK 
    %n_area(cell) = strcmp(units(cell).area, recArea);
    monk(cell) = strcmp(units(cell).monk, monkey);
end

%% PICK 
%indx_monk = find(n_area & monk); % recArea and monk
indx_monk = find(monk); %  monk
%indx_monk = find(n_area); % recArea and monk

%%
for cellNum = 1:length(indx_monk)
    %prosacc
    for trialNum = 1:length([units(indx_monk(cellNum)).pro.behav.trial.saccAmplitude])
        eyeKin(1,cellNum).proAmp(trialNum,1) = units(indx_monk(cellNum)).pro.behav.trial(trialNum).saccAmplitude;
        eyeKin(1,cellNum).proDur(trialNum,1) = units(indx_monk(cellNum)).pro.behav.trial(trialNum).saccDuration;
        eyeKin(1,cellNum).proPV(trialNum,1) = units(indx_monk(cellNum)).pro.behav.trial(trialNum).saccPeakVel;
        eyeKin(1,cellNum).proRT(trialNum,1) = units(indx_monk(cellNum)).pro.behav.trial(trialNum).reactionTime;
        
    end
    %antisacc
     for trialNum = 1:length([units(indx_monk(cellNum)).anti.behav.trial.saccAmplitude])
        eyeKin(1,cellNum).antiAmp(trialNum,1) = units(indx_monk(cellNum)).anti.behav.trial(trialNum).saccAmplitude;
        eyeKin(1,cellNum).antiDur(trialNum,1) = units(indx_monk(cellNum)).anti.behav.trial(trialNum).saccDuration;
        eyeKin(1,cellNum).antiPV(trialNum,1) = units(indx_monk(cellNum)).anti.behav.trial(trialNum).saccPeakVel;
        eyeKin(1,cellNum).antiRT(trialNum,1) = units(indx_monk(cellNum)).anti.behav.trial(trialNum).reactionTime;
     end
end

proAmp = vertcat(eyeKin(1,:).proAmp); 
antiAmp = vertcat(eyeKin(1,:).antiAmp); 

proDur = vertcat(eyeKin(1,:).proDur); 
antiDur = vertcat(eyeKin(1,:).antiDur); 

proPV = vertcat(eyeKin(1,:).proPV); 
antiPV = vertcat(eyeKin(1,:).antiPV); 

proRT = vertcat(eyeKin(1,:).proRT)*1000; 
antiRT = vertcat(eyeKin(1,:).antiRT)*1000; 

