function createTable_pro_anti(units)
% Create stats table
total_cells = length(units);

%% vermis 
for cellNum = 1:length(units)
    indx_area(cellNum) = strcmp(units(cellNum).area, 'vermis');
end
indx_vermis = find(indx_area);
total_vermis = length(indx_vermis); %%%%%%%%%%%%% to table

% how many significant diff for pro vs anti and base vs windows

for i=1:length(indx_vermis)
indx_sign_sacc(i)= units(indx_vermis(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk; 
indx_sign_instr(i)= units(indx_vermis(i)).stats.instr.flags.proVsAnti_instr_ks_nspk; 
indx_sign_sacc_base_pro(i)= units(indx_vermis(i)).stats.pro.flags.saccVSbase_nspk; 
indx_sign_instr_base_pro(i)= units(indx_vermis(i)).stats.pro.flags.instrVSbase_nspk;
indx_sign_sacc_base_anti(i)= units(indx_vermis(i)).stats.anti.flags.saccVSbase_nspk; 
indx_sign_instr_base_anti(i)= units(indx_vermis(i)).stats.anti.flags.instrVSbase_nspk;
end 
vermis_diff_pro_anti_sacc = sum(indx_sign_sacc); %%%%%%%%%%%%% to table
proportion_vermis_sacc = vermis_diff_pro_anti_sacc/length(indx_vermis); %%%%%%%%%%%%% to table
vermis_diff_pro_anti_instr = sum(indx_sign_instr); %%%%%%%%%%%%% to table
proportion_vermis_instr = vermis_diff_pro_anti_instr/length(indx_vermis); %%%%%%%%%%%%% to table
vermis_diff_sacc_base_pro = sum(indx_sign_sacc_base_pro); %%%%%%%%%%%%% to table
vermis_diff_instr_base_pro = sum(indx_sign_instr_base_pro); %%%%%%%%%%%%% to table
vermis_diff_sacc_base_anti = sum(indx_sign_sacc_base_anti); %%%%%%%%%%%%% to table
vermis_diff_instr_base_anti = sum(indx_sign_instr_base_anti); %%%%%%%%%%%%% to table

%% lateral
clear indx_area
for cellNum = 1:length(units)
    indx_area(cellNum) = strcmp(units(cellNum).area, 'lateral');
end
indx_lateral = find(indx_area);
total_lateral = length(indx_lateral); %%%%%%%%%%%%% to table

for i=1:length(indx_lateral)
indx_sign_sacc_lat(i)= units(indx_lateral(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk; 
indx_sign_instr_lat(i)= units(indx_lateral(i)).stats.instr.flags.proVsAnti_instr_ks_nspk; 
indx_sign_sacc_base_pro_lat(i)= units(indx_lateral(i)).stats.pro.flags.instrVSbase_nspk; 
indx_sign_instr_base_pro_lat(i)= units(indx_lateral(i)).stats.pro.flags.instrVSbase_nspk;
indx_sign_sacc_base_anti_lat(i)= units(indx_lateral(i)).stats.anti.flags.instrVSbase_nspk; 
indx_sign_instr_base_anti_lat(i)= units(indx_lateral(i)).stats.anti.flags.instrVSbase_nspk;
end 

lateral_diff_pro_anti_sacc = sum(indx_sign_sacc_lat); %%%%%%%%%%%%% to table
proportion_lateral_sacc = lateral_diff_pro_anti_sacc/length(indx_lateral); %%%%%%%%%%%%% to table
lateral_diff_pro_anti_instr = sum(indx_sign_instr_lat); %%%%%%%%%%%%% to table
proportion_lateral_instr = lateral_diff_pro_anti_instr/length(indx_lateral); %%%%%%%%%%%%% to table
lateral_diff_sacc_base_pro = sum(indx_sign_sacc_base_pro_lat); %%%%%%%%%%%%% to table
lateral_diff_instr_base_pro = sum(indx_sign_instr_base_pro_lat); %%%%%%%%%%%%% to table
lateral_diff_sacc_base_anti = sum(indx_sign_sacc_base_anti_lat); %%%%%%%%%%%%% to table
lateral_diff_instr_base_anti = sum(indx_sign_instr_base_anti_lat); %%%%%%%%%%%%% to table

% print table
disp('>>> VERMIS <<<')
disp(['Total num of cells = ' num2str(total_cells)]); 
disp(['Total num vermis = ' num2str(total_vermis)]); 
disp(['Sign pro vs anti instr = ' num2str(vermis_diff_pro_anti_instr)]);
disp(['Proportion sign pro vs anti instr = ' num2str(proportion_vermis_instr)]);
disp(['Sign pro vs anti sacc = ' num2str(vermis_diff_pro_anti_sacc)]);
disp(['Proportion sign pro vs anti sacc = ' num2str(proportion_vermis_sacc)]);
disp(['Sign diff base vs instr pro = ' num2str(vermis_diff_instr_base_pro)]);
disp(['Sign diff base vs sacc pro = ' num2str(vermis_diff_sacc_base_pro)]);
disp(['Sign diff base vs instr anti = ' num2str(vermis_diff_instr_base_anti)]);
disp(['Sign diff base vs sacc anti = ' num2str(vermis_diff_sacc_base_anti)]);
disp('                             ')
disp('>>> LATERAL <<<')
disp(['Total num of cells = ' num2str(total_cells)]); 
disp(['Total num lateral = ' num2str(total_lateral)]); 
disp(['Sign pro vs anti instr = ' num2str(lateral_diff_pro_anti_instr)]);
disp(['Proportion sign pro vs anti instr = ' num2str(proportion_lateral_instr)]);
disp(['Sign pro vs anti sacc = ' num2str(lateral_diff_pro_anti_sacc)]);
disp(['Proportion sign pro vs anti sacc = ' num2str(proportion_lateral_sacc)]);
disp(['Sign diff base vs instr pro = ' num2str(lateral_diff_instr_base_pro)]);
disp(['Sign diff base vs sacc pro = ' num2str(lateral_diff_sacc_base_pro)]);
disp(['Sign diff base vs instr anti = ' num2str(lateral_diff_instr_base_anti)]);
disp(['Sign diff base vs sacc anti = ' num2str(lateral_diff_sacc_base_anti)]);

