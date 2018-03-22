
[filename,path,~]=uigetfile()

filename= [path,filename]
% num = xlsread(filename)
[num,txt,raw] = xlsread(filename)

neurons_depths = cell(size(raw))
for i = 1:size(neurons_depths,1)
    neurons_depths{i,1} =  datetime(txt{i,1}(2:end-1),'InputFormat','yy/MM/dd');
    neurons_depths{i,2} = txt{i,2}(2:end-1);
    if ~isnan(raw{i,3})
        neurons_depths{i,3} = raw{i,3};
    else
        neurons_depths{i,3} = [];
    end
    
    
    if ~isnan(raw{i,4})
        neurons_depths{i,4} = str2num(raw{i,4});
    else
        neurons_depths{i,4} =[];
    end
end

nans = (cellfun(@isnan, neurons_depths(:,[3,4]),'UniformOutput',false))

for i =1:size(neuronList,2)
    
    for ii=1:size(neurons_depths,1)
        if strcmp(neuronList(i).neuronName,neurons_depths{ii,2})
             
            neuronList(i).depth = neurons_depths{ii,3};
            neuronList(i).gridLocation = neurons_depths{ii,4};
        end 
    end
end

nans = cell2mat(cellfun(@isnan, raw(:,4),'UniformOutput',false))


