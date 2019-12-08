function [names] = extract_neuron_names(neuronList)

disp('          extracting neuron names...')
for i = 1:length(neuronList)
    if ~isempty(neuronList(i))
        neuron_id = neuronList(i).neuronName;
        
        fid=fopen('neuron_list', 'at');
        fprintf(fid,'\n');
        
        fprintf(fid, '\n%s ', neuron_id);
        fclose(fid);
    end
end

disp('  >>>>   Done    <<<<  ')

end 