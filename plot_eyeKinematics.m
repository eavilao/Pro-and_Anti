% Saccade kinematics

%IN = trialData.m
% This code takes in trialData and plots eye data and computes stats for
% all neurons. 

%gather
for cellNum = 1:length(trialData)
    %prosacc
    for trialNum = 1:length([trialData(cellNum).pro.trial.saccAmplitude])
        eyeKin(1,cellNum).proAmp(trialNum,1) = trialData(cellNum).pro.trial(trialNum).saccAmplitude;
        eyeKin(1,cellNum).proDur(trialNum,1) = trialData(cellNum).pro.trial(trialNum).saccDuration;
        eyeKin(1,cellNum).proPV(trialNum,1) = trialData(cellNum).pro.trial(trialNum).saccPeakVel;
    end
    %antisacc
     for trialNum = 1:length([trialData(cellNum).anti.trial.saccAmplitude])
        eyeKin(1,cellNum).antiAmp(trialNum,1) = trialData(cellNum).anti.trial(trialNum).saccAmplitude;
        eyeKin(1,cellNum).antiDur(trialNum,1) = trialData(cellNum).anti.trial(trialNum).saccDuration;
        eyeKin(1,cellNum).antiPV(trialNum,1) = trialData(cellNum).anti.trial(trialNum).saccPeakVel;
     end
end

proAmp = vertcat(eyeKin(1,:).proAmp); 
antiAmp = vertcat(eyeKin(1,:).antiAmp); 

proDur = vertcat(eyeKin(1,:).proDur); 
antiDur = vertcat(eyeKin(1,:).antiDur); 

proPV = vertcat(eyeKin(1,:).proPV); 
antiPV = vertcat(eyeKin(1,:).antiPV); 

%% plot

% plot amp
figure; hold on
h1 = histfit(proAmp,20,'kernel');
h2 = histfit(antiAmp,20,'kernel');
xlabel('Saccade amplitude (deg)')
ylabel('Number of trials')
set (gca, 'TickDir', 'out','FontSize', 18);
alpha(0.25)
set(h1(1),'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
set(h2(1),'FaceColor', [0 1 0],'EdgeColor', [0 1 0]);
set(h1(2),'Color',[1 0 0]);
set(h2(2),'Color',[0 1 0]);

%plot dur
figure; hold on
h1 = histfit(proDur,20,'kernel');
h2 = histfit(antiDur,20,'kernel');
xlabel('Saccade duration (ms)')
ylabel('Number of trials')
set (gca, 'TickDir', 'out','FontSize', 18);
alpha(0.25)
set(h1(1),'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
set(h2(1),'FaceColor', [0 1 0],'EdgeColor', [0 1 0]);
set(h1(2),'Color',[1 0 0]);
set(h2(2),'Color',[0 1 0]);

% plot PV
figure; hold on
h1 = histfit(proPV,20,'kernel');
h2 = histfit(antiPV,20,'kernel');
xlabel('Saccade peak velocity (deg/s)')
ylabel('Number of trials')
set (gca, 'TickDir', 'out','FontSize', 18);
alpha(0.25)
set(h1(1),'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
set(h2(1),'FaceColor', [0 1 0],'EdgeColor', [0 1 0]);
set(h1(2),'Color',[1 0 0]);
set(h2(2),'Color',[0 1 0]);