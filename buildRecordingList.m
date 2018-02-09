% Build recording list 

% collect all files and folders 
clear
allfiles = getAllFiles(cd);


% get all folders with Neurons: taged on N1, N2 Nn..,. 
folders = regexp(allfiles(:,1),'N[123456789]+');
index_neurons = find(~cellfun(@isempty,folders));

 folders(cellfun(@isempty,folders))=[];
%  Nindex = cell2mat(folders)

% get marked stimlist folders
 stimListFolders = allfiles(~cellfun(@isempty, regexpi(allfiles,'stimListFiles\')));
 % prune strings to get real folders
for i = 1: size(stimListFolders,1)
    slash = strfind(stimListFolders{i},'stimListFiles\');
    stimListFolders{i,:} = stimListFolders{i}(1:slash+13);
end
stimListFolders=unique(stimListFolders);
  
% get all containing files
stimfolder=cell(4,1);

for ii = 1:4
    stimfolder{ii} = getAllFiles(stimListFolders{ii});
end
stimfiles = vertcat(stimfolder{:,:});
stimfiles(  ~cellfun(@isempty, regexpi(stimfiles,'\._')))=[];

clearvars stimfolder

if ~isempty(index_neurons)
allfiles = allfiles(index_neurons);
end



% cut foldernames to only keep folders that are unique to cells
cell_folders=cell(size(allfiles));
for i = 1: size(allfiles,1)
    cell_folders{i} = allfiles{i}(1:folders{i}+1);
end
cell_folders = unique(cell_folders);
clearvars allfiles index_neurons

% throw out cell in spikeSorted directory to prevent double counting these
% cells
cell_folders(~cellfun(@isempty, regexpi(cell_folders,'spikesorted')))=[];
cell_folders(~cellfun(@isempty, regexpi(cell_folders,'stimlist')))=[];
 
for i = 1: size(cell_folders,1)
    slash = strfind(cell_folders{i},'\');
    cell_names{i,:} = cell_folders{i}(slash(end)+1:end);
end

if size(cell_names) ~= size(unique(cell_names))
    disp('duplicate cells detected')
    keyboard
end


% first parts of neuron list: path to data, cellname,monkey name, area
% for i= 1:size(cell_names,1)
    for i= 1:10

    % path and name of cells
    neuronList(i).path = cell_folders{i};
    neuronList(i).neuronName = cell_names{i};
    
    % detect area from path string
    if regexpi(cell_folders{i},'Vermis')
        neuronList(i).area = 'vermis';
    elseif regexpi(cell_folders{i},'lateral')
        neuronList(i).area = 'lateral';
    else 
        disp('no area detected')
        keyboard
    end
    
    % detect monkey name from path string 
    if regexpi(cell_folders{i},'moshe')
        neuronList(i).monkey = 'moshe';
    elseif regexpi(cell_folders{i},'mickey')
        neuronList(i).monkey = 'mickey';
    else 
        disp('no monkey name detected')
        keyboard
    end

    
    % now get the date of the recording
    year = regexp(cell_names{i},'1+[234567]+'); % 2012 - 2017
    day = regexp(cell_names{i},'[123]+[0123456789]+') ;% 1 - 31
    daystring = cell_names{i}(year(1):day(end)+1);
    try
    recDay = datetime(daystring,'InputFormat','yy_MM_dd');
    catch
        keyboard
    end 
    neuronList(i).recordingDate = recDay;
    
    % compare recording date with newer (nico) recordings, these have
    % eyechannels on different input. Newer than 2016-01-01 
    if recDay > datetime('16_01_01','InputFormat','yy_MM_dd')
        neuronList(i).eyeChannels = [1,2];
    else 
        neuronList(i).eyeChannels = [];
    end
        
    if exist( cell_folders{i})==7
    recordings = getAllFiles(cell_folders{i});
    elseif exist([cell_folders{i},'ProAndAnti'])
    recordings =   getAllFiles([cell_folders{i},'ProAndAnti']);
    else
        disp('cant locate folder:')
        disp(cell_folders{i})
        keyboard
    end
    
    % find all rsd files
     rsdfiles = recordings;
     matfiles = recordings;
     matfiles(cellfun(@isempty, regexpi(matfiles,'mat')))=[];
     rsdfiles(cellfun(@isempty, regexpi(rsdfiles,'rsd')))=[];
     
     if isempty(rsdfiles)
         % no rsd files throw neuron from neuronList - add to error list 
         disp('cell not spike-sorted')
         error{i,1}= cell_folders{i};
         error{i,2}= 'not sorted';
         
         neuronList(i)=[];
         continue
     end
     
     % get recording channel number from .rsd string
     str_rsd = cell2mat(regexpi(rsdfiles(1),'.rsd'));
     recording_channel = str2num(rsdfiles{1}(str_rsd-2:str_rsd));
     
     for ii = 1:size(rsdfiles,1)
         rsd{ii,:}=rsdfiles{ii}(1:end-7);
     end
     for ii = 1:size(matfiles,1)
         mat{ii,:}=matfiles{ii}(1:end-4);
     end

     
     for ii = 1:size(rsd,1)
       
         % indexes of rsd files  in matfile-list e.g. ind == 1 matfile(1) == rsd 1 
         ind(ii)  = find(strcmp(mat,rsd{ii}));
%         
     end
      
     % slightly messy indexing 
     index = zeros(size(mat,1),1);
     index(ind)=1;
         spikeMatFiles = mat(ind);
          fileList = cell(size(spikeMatFiles,1),7);

 for ii = 1:size(spikeMatFiles,1)
     fileSeperator = strfind(spikeMatFiles,'\');
     fileList{ii,1} =  [spikeMatFiles{1}(fileSeperator{1}(end)+1:end),'.mat'];
     fileList{ii,2} = recording_channel;
     fileList{ii,3} = [1,2];
     fileList{ii,6} = nan;
 end
     neuronList(i).fileList = fileList; 

% neuronList(39).neuronName = '17_01_20_N2';    % 
% neuronList(39).fileList = {
% % File Name                 spk channel number  [unit numbers] startTime stopTime stimStartFile resortFlag
% 'F170120-0005.mat'          6              [1 2]     []          []      NaN    0
% };
% neuronList(39).group = 'Pause_Burst'; %not real 
% neuronList(39).area = 'vermis';
% neuronList(39).depth = 31.6; % real
% neuronList(39).gridLocation = [23 5]; %real
% neuronList(39).eyeChannels = [1 2];

% get all files that end in [_number].txt 
% folders3 = regexp(allfiles(:,1),'_[123456789]+.txt');
% index_txt = find(~cellfun(@isempty,folders3));
% allfiles = allfiles(index_txt);

     otherMatFiles = cell(size(mat(index==0),1),3);
     otherMatFiles(:,1) =  mat(index==0);
     
     % search non-spike matfiles for stimList
     if ~isempty(otherMatFiles)
         for ii= i:size(otherMatFiles,1)
             
             % get variables that are present in remaining mat files 
             thisFile = matfile(otherMatFiles{ii,1});
             thisFileInfo = whos(thisFile);
             variablenames =  [thisFileInfo(:).name];
                          
             out = regexpi(variablenames,{'ans','epd','stimlist'});
             
             if ~isempty(cell2mat(out)) && isempty( out{1}) % ans could potential be in a string due to coincidence
                 otherMatFiles{ii,2} = 'stimListFile';
                 stimListInfo = dir([otherMatFiles{ii,1},'.mat']); 
                 otherMatFiles{ii,3} = stimListInfo.bytes; % add filesize to later find largest stimList(with most trials)
             elseif ~isempty(out{1})  && strcmp(thisFileInfo(1).name,'ans' )% if ans has been found validate that its really a stimlist
                 
                 
                 otherMatFiles{ii,2} = 'stimListFile'; 
                 stimListInfo = dir([otherMatFiles{ii,1},'.mat']);
                 otherMatFiles{ii,3} = stimListInfo.bytes;
             end
         end
     end
      
     
     candidates_placeholder = stimfiles;
     % find files in stimListFolders that contain the right date
     candidates_placeholder(cellfun(@isempty, regexpi(stimfiles,daystring)))=[];
     % prepare empty cell array that can later be concatinated
     candidates= cell(size(candidates_placeholder,1),3);
     candidates(:,1) = candidates_placeholder;
     
     for ii= i:size(candidates,1)
         
         % get variables that are present in remaining mat files
         thisFile = matfile(candidates{ii,1});
         thisFileInfo = whos(thisFile);
         variablenames =  [thisFileInfo(:).name];
         out = regexpi(variablenames,{'ans','epd','stimlist'});
         if ~isempty(cell2mat(out)) && isempty( out{1}) % ans could potential be in a string due to coincidence
             candidates{ii,2} = 'stimListFile';
             stimListInfo = dir([candidates{ii,1} ]);
             candidates{ii,3} = stimListInfo.bytes; % add filesize to later find largest stimList(with most trials)
         elseif ~isempty(out{1})  && strcmp(thisFileInfo(1).name,'ans' )% if ans has been found validate that its really a stimlist
             
             
             candidates{ii,2} = 'stimListFile';
             stimListInfo = dir([candidates{ii,1} ]);
             candidates{ii,3} = stimListInfo.bytes;
         end
         
     end
     
     possibleStimListFiles = [candidates;otherMatFiles];
     
     realCandidates = [possibleStimListFiles{:,2}];
     if isempty(realCandidates)
         % no rsd files throw neuron from neuronList - add to error list 
         disp('no stimlistfile')
         error{i,1}= cell_folders{i};
         error{i,2}= 'no stimListFile';
         
         neuronList(i)=[];
         continue
     end
     [m, mi] =  max([possibleStimListFiles{:,3}]);
     neuronList(i).stimListFile = possibleStimListFiles{mi,1};
%     take rsd - find matching mat files - these are the spike files. 
%     
%      any leftover files are likely the stimlistfile -> check by seeing if it contains variable: epd or stimlistfile or ans. 
%     find biggest stimlist file in mb's 
    
end


% eyeChannels = [1 2] % if recording newer then 2016

% categories=  dir(cd);
% 
% for i = 3:size(categories,1)
%     monkey = regexp(
%     location 