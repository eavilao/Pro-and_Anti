
neuronListverm= neuronListVerm.neuronList
neuronCell  = squeeze(struct2cell(neuronList));
neuronCellVerm = squeeze(struct2cell(neuronListverm));



% loop over all cell fields and match neuronName, then copy grid and depth
for i = 1:size(neuronCellVerm,2)
    for ii= size(neuronCell)
        if strcmp(neuronname , neuronname verm)
            grid loc 
            depth
        end 
    end
end

celly2(:,1)=neuronCell(5,:)
xx= find(~cellfun(@isempty, celly2(:,1)))
formatOut = 'yy/mm/dd';
n=1
for i = xx'
    cells{n,1} = datestr(celly2{i,1},formatOut);
    cells{n,2} = celly2{i,2};
    n=n+1;
end

sortrows(cells,1)