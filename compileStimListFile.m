% compile stimlist
% make one stimlistfile in old recordings that have all trials saved as
% seperate files


clear

listdir = uigetdir;

% collect all trials
trials = getAllMatFiles(listdir);

% index usefull and useless trails
for i = 1:size(trials,1)
    
  indd(i,:) =  ~any(strfind(trials{i},'._')); % zero's are bad trials
  
end

% only use good files
goodtrials =trials(indd==1);


% append all trials in the known struct-file
for i = 1:size(goodtrials,1)
    load(goodtrials{i}, 'trialData');
    
    if exist('trialData') 
    stimListFile(i) = trialData;
    end
    
    clearvars trialData
end

% extract the date of recording

datestring =  stimListFile(1).hour(1:8);
try
datestr(datestring,'yymmdd')
catch
    
    disp('improper date')
    keyboard
end

% makeup a senisble file name 
if ~isempty(regexpi(listdir,'mickey'))
    filename = ['Mickey_', datestring, '-stimList'];
    save(filename,'stimListFile')
      disp(['saved as ' filename])

elseif ~isempty(regexpi(listdir,'moshe'))
     filename = ['Moshe_', datestring, '-stimList'];
    save(filename,'stimListFile')
    disp(['saved as ' filename])
else
    keyboard
end

    clear
    
        

