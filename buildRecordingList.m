% Build recording list 

% collect all files and folders 

allfiles = getAllFiles(cd);


% get all folders with Neurons: taged on N1, N2 Nn..,. 
folders = regexp(allfiles(:,1),'N[123456789]+');
index_neurons = find(~cellfun(@isempty,folders));

 folders(cellfun(@isempty,folders))=[];
%  Nindex = cell2mat(folders)
 
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
for i= 1:size(cell_names,1)
    
    % path and name of cells
    neuronList(i).path = cell_folders{i};
    neuronList(i).name = cell_names{i};
    
    % detect area from path string
    if regexpi(cell_folder{i},'Vermis')
        neuronList(i).area = 'vermis';
    elseif regexpi(cell_folder{i},'lateral')
        neuronList(i).area = 'lateral';
    else 
        disp('no area detected')
        keyboard
    end
    
    % detect monkey name from path string 
    if regexpi(cell_folder{i},'moshe')
        neuronList(i).monkey = 'moshe';
    elseif regexpi(cell_folder{i},'mickey')
        neuronList(i).area = 'mickey';
    else 
        disp('no monkey name detected')
        keyboard
    end
        
        
       
end


% get all files that end in [_number].txt 
% folders3 = regexp(allfiles(:,1),'_[123456789]+.txt');
% index_txt = find(~cellfun(@isempty,folders3));
% allfiles = allfiles(index_txt);

% categories=  dir(cd);
% 
% for i = 3:size(categories,1)
%     monkey = regexp(
%     location 