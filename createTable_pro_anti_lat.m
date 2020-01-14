function selected = createTable_pro_anti_lat(units)

clear indx_area
for cellNum = 1:length(units)
    indx_area(cellNum) = strcmp(units(cellNum).area, 'lateral');
end
indx_lateral = find(indx_area);
total_lateral = length(indx_lateral); %%%%%%%%%%%%% to table

for i=1:length(indx_lateral)
indx_sign_sacc(i)= units(indx_lateral(i)).stats.sacc.flags.proVsAnti_sacc_ks_nspk; 
indx_sign_instr(i)= units(indx_lateral(i)).stats.instr.flags.proVsAnti_instr_ks_nspk; 
indx_sign_sacc_base_pro(i)= units(indx_lateral(i)).stats.pro.flags.saccVSbase_nspk; 
indx_sign_instr_base_pro(i)= units(indx_lateral(i)).stats.pro.flags.instrVSbase_nspk;
indx_sign_sacc_base_anti(i)= units(indx_lateral(i)).stats.anti.flags.saccVSbase_nspk;
indx_sign_instr_base_anti(i)= units(indx_lateral(i)).stats.anti.flags.instrVSbase_nspk;
indx_sign_instrDir_instr_pro(i) = units(indx_lateral(i)).stats.pro.flags.instrDirVSinstr_nspk;
indx_sign_instrDir_instr_anti(i) = units(indx_lateral(i)).stats.anti.flags.instrDirVSinstr_nspk;
indx_sign_goCue_instrDir_pro(i) = units(indx_lateral(i)).stats.pro.flags.goCueVSinstrDir_nspk;
indx_sign_goCue_instrDir_anti(i) = units(indx_lateral(i)).stats.anti.flags.goCueVSinstrDir_nspk;
indx_sign_saccVSinstrDir_nspk_pro(i) = units(indx_lateral(i)).stats.pro.flags.saccVSinstrDir_nspk;
indx_sign_saccVSinstrDir_nspk_anti(i) = units(indx_lateral(i)).stats.anti.flags.saccVSinstrDir_nspk;
indx_exc_sacc_pro(i) = units(indx_lateral(i)).pro.neural.exc; 
indx_sup_sacc_pro(i) = units(indx_lateral(i)).pro.neural.sup; 
indx_exc_sacc_anti(i) = units(indx_lateral(i)).anti.neural.exc; 
indx_sup_sacc_anti(i) = units(indx_lateral(i)).anti.neural.sup;
indx_exc_instr_pro(i) = units(indx_lateral(i)).pro.neural.exc_instr;
indx_exc_instr_anti(i) = units(indx_lateral(i)).anti.neural.exc_instr;
indx_sup_instr_pro(i) = units(indx_lateral(i)).pro.neural.sup_instr;
indx_sup_instr_anti(i) = units(indx_lateral(i)).anti.neural.sup_instr;
indx_sacc_pro_exc_sup_sign(i) = units(indx_lateral(i)).pro.neural.categ.flag;
indx_sacc_anti_exc_sup_sign(i) = units(indx_lateral(i)).anti.neural.categ.flag;
pval_sacc_pro_exc_sup_sign(i) = units(indx_lateral(i)).pro.neural.categ.pval;
pval_sacc_anti_exc_sup_sign(i) = units(indx_lateral(i)).anti.neural.categ.pval;
indx_instr_pro(i) = units(indx_lateral(i)).stats.pro.flags.instrVSbase_nspk; 
indx_instr_anti(i) = units(indx_lateral(i)).stats.anti.flags.instrVSbase_nspk; 
indx_instr_back_pro(i) = units(indx_lateral(i)).stats.pro.flags.instr_backVSbase_nspk; 
indx_instr_back_anti(i) = units(indx_lateral(i)).stats.anti.flags.instr_backVSbase_nspk; 

num_corr_pro(i) = length(units(i).pro.indx_correctProTrials);  % calculate total number of trials for pro and anti
num_corr_anti(i) = length(units(i).anti.indx_correctAntiTrials);  
end 

pro_only_instr_lat = sum(indx_sign_instr_base_pro & ~indx_sign_instr_base_anti);
anti_only_instr_lat = sum(indx_sign_instr_base_anti & ~indx_sign_instr_base_pro);
pro_anti_lat = sum(indx_sign_instr_base_anti & indx_sign_instr_base_pro);
pro_only_sacc_lat = sum(indx_sign_sacc_base_pro & ~indx_sign_sacc_base_anti);
% pro_only_sacc_exc_omv = sum(indx_exc_sacc_pro & ~indx_sign_sacc_base_anti & indx_sacc_pro_exc_sup_sign);  indx_pro_only_sacc_exc_omv = indx_exc_sacc_pro & ~indx_sign_sacc_base_anti & indx_sacc_pro_exc_sup_sign;
% pro_only_sacc_sup_omv = sum(indx_sup_sacc_pro & ~indx_sign_sacc_base_anti  & indx_sacc_pro_exc_sup_sign);  indx_only_sacc_sup_omv = indx_sup_sacc_pro & ~indx_sign_sacc_base_anti  & indx_sacc_pro_exc_sup_sign;
anti_only_sacc_lat = sum(indx_sign_sacc_base_anti & ~indx_sign_sacc_base_pro);
pro_anti_sacc_lat = sum(indx_sign_sacc_base_anti & indx_sign_sacc_base_pro);
pro_only_instrDir_lat = sum(indx_sign_instrDir_instr_pro & ~indx_sign_instrDir_instr_anti);
anti_only_instrDir_lat = sum(indx_sign_instrDir_instr_anti & ~indx_sign_instrDir_instr_pro);
pro_anti_instrDir_lat = sum(indx_sign_instrDir_instr_pro & indx_sign_instrDir_instr_anti);
pro_only_goCue_lat = sum(indx_sign_goCue_instrDir_pro & ~indx_sign_goCue_instrDir_anti);
anti_only_goCue_lat = sum(indx_sign_goCue_instrDir_anti & ~indx_sign_goCue_instrDir_pro);
pro_anti_goCue_lat = sum(indx_sign_goCue_instrDir_pro & indx_sign_goCue_instrDir_anti);

p_val = 0.05; 
% pro_only_sacc_exc_omv = sum(indx_exc_sacc_pro & ~indx_sign_sacc_base_anti & indx_sacc_pro_exc_sup_sign);  indx_pro_only_sacc_exc_omv = indx_exc_sacc_pro & ~indx_sign_sacc_base_anti & indx_sacc_pro_exc_sup_sign;
% pro_only_sacc_exc_omv =  sum((pval_sacc_pro_exc_sup_sign<p_val) & ~(pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_pro);  indx_pro_only_sacc_exc_omv = (pval_sacc_pro_exc_sup_sign<p_val) & ~(pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_pro;
% pro_only_sacc_sup_omv = sum((pval_sacc_pro_exc_sup_sign<p_val) & ~(pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro);  indx_pro_only_sacc_sup_omv = (pval_sacc_pro_exc_sup_sign<p_val) & ~(pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro;
% anti_only_sacc_exc_omv = sum(~(pval_sacc_pro_exc_sup_sign<p_val) & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_anti);  indx_anti_only_sacc_exc_omv = ~(pval_sacc_pro_exc_sup_sign<p_val) & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_anti;
% anti_only_sacc_sup_omv = sum(~(pval_sacc_pro_exc_sup_sign<p_val) & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_anti);  indx_anti_only_sacc_sup_omv = ~(pval_sacc_pro_exc_sup_sign<p_val) & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_anti;
% pro_anti_sacc_exc_omv = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & ~indx_sup_sacc_anti);  indx_pro_anti_sacc_exc_omv = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & ~indx_sup_sacc_anti; 
% pro_anti_sacc_sup_omv = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro & indx_sup_sacc_anti);  indx_pro_anti_sacc_sup_omv = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro & indx_sup_sacc_anti;   

pro_only_sacc_exc_lat =  sum((pval_sacc_pro_exc_sup_sign<p_val) & ~indx_sacc_anti_exc_sup_sign & ~indx_sup_sacc_pro);  indx_pro_only_sacc_exc_lat = (pval_sacc_pro_exc_sup_sign<p_val) & ~indx_sacc_anti_exc_sup_sign & ~indx_sup_sacc_pro;
pro_only_sacc_sup_lat = sum((pval_sacc_pro_exc_sup_sign<p_val) & ~indx_sacc_anti_exc_sup_sign & ~indx_exc_sacc_pro);  indx_pro_only_sacc_sup_lat = (pval_sacc_pro_exc_sup_sign<p_val) & ~indx_sacc_anti_exc_sup_sign & ~indx_exc_sacc_pro;
anti_only_sacc_exc_lat = sum(~indx_sacc_pro_exc_sup_sign & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_anti);  indx_anti_only_sacc_exc_lat = ~indx_sacc_pro_exc_sup_sign & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_sup_sacc_anti;
anti_only_sacc_sup_lat = sum(~indx_sacc_pro_exc_sup_sign & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_anti);  indx_anti_only_sacc_sup_lat = ~indx_sacc_pro_exc_sup_sign & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_anti;
pro_anti_sacc_exc_lat = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & ~indx_sup_sacc_anti);  indx_pro_anti_sacc_exc_lat = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & ~indx_sup_sacc_anti; 
pro_anti_sacc_sup_lat = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro & indx_sup_sacc_anti);  indx_pro_anti_sacc_sup_lat = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & ~indx_exc_sacc_pro & indx_sup_sacc_anti;   
both_anti_exc_sup_lat = sum(pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & indx_sup_sacc_anti); indx_mixed = pval_sacc_pro_exc_sup_sign<p_val & (pval_sacc_anti_exc_sup_sign<p_val) & indx_exc_sacc_pro & indx_sup_sacc_anti; 

selected.sacc.pro.exc = indx_lateral(indx_pro_only_sacc_exc_lat); 
selected.sacc.pro.sup = indx_lateral(indx_pro_only_sacc_sup_lat); 
selected.sacc.anti.exc = indx_lateral(indx_anti_only_sacc_exc_lat);
selected.sacc.anti.sup = indx_lateral(indx_anti_only_sacc_sup_lat);
selected.sacc.both.exc = indx_lateral(indx_pro_anti_sacc_exc_lat);
selected.sacc.both.sup = indx_lateral(indx_pro_anti_sacc_sup_lat);
selected.sacc.both.mixed = indx_lateral(indx_mixed);
selected.sacc.all.pro.exc = indx_lateral(pval_sacc_pro_exc_sup_sign<p_val & ~indx_sup_sacc_pro); 
selected.sacc.all.anti.exc = indx_lateral(pval_sacc_anti_exc_sup_sign<p_val & ~indx_sup_sacc_anti); 
selected.sacc.all.pro.sup = indx_lateral(pval_sacc_pro_exc_sup_sign<p_val & ~indx_exc_sacc_pro); 
selected.sacc.all.anti.sup = indx_lateral(pval_sacc_anti_exc_sup_sign<p_val & ~indx_exc_sacc_anti);

%% instruction too
instr_pro_only_sacc_exc_lat = sum(indx_pro_only_sacc_exc_lat & indx_instr_back_pro); 
instr_pro_only_sacc_sup_lat = sum (indx_pro_only_sacc_sup_lat & indx_instr_back_pro); 
instr_anti_only_sacc_exc_lat = sum(indx_anti_only_sacc_exc_lat & indx_instr_back_anti); 
instr_anti_only_sacc_sup_lat = sum(indx_anti_only_sacc_sup_lat & indx_instr_back_anti); 
instr_pro_anti_sacc_exc_lat = sum(indx_pro_anti_sacc_exc_lat & indx_instr_pro & indx_instr_back_anti); 
instr_pro_anti_sacc_sup_lat = sum(indx_pro_anti_sacc_sup_lat & indx_instr_pro & indx_instr_back_anti); 
%instr_mixed = sum(indx_mixed & indx_instr_back_anti & indx_instr_back_pro); 

selected.instr_back.pro.exc = indx_lateral(indx_pro_only_sacc_exc_lat & indx_instr_back_pro); 
selected.instr_back.pro.sup = indx_lateral(indx_pro_only_sacc_sup_lat & indx_instr_back_pro); 
selected.instr_back.anti.exc = indx_lateral(indx_anti_only_sacc_exc_lat & indx_instr_back_anti);
selected.instr_back.anti.sup = indx_lateral(indx_anti_only_sacc_sup_lat & indx_instr_back_anti);
selected.instr_back.both.exc = indx_lateral(indx_pro_anti_sacc_exc_lat & indx_instr_back_pro & indx_instr_back_anti);
selected.instr_back.both.sup = indx_lateral(indx_pro_anti_sacc_sup_lat & indx_instr_back_pro & indx_instr_back_anti);
selected.instr_back.both.mixed = indx_lateral(indx_mixed);
selected.instr_back.all.pro.exc = indx_lateral(pval_sacc_pro_exc_sup_sign<p_val & ~indx_sup_sacc_pro & indx_instr_back_pro);
selected.instr_back.all.anti.exc = indx_lateral(pval_sacc_anti_exc_sup_sign<p_val & ~indx_sup_sacc_anti & indx_instr_back_anti);
selected.instr_back.all.pro.sup = indx_lateral(pval_sacc_pro_exc_sup_sign<p_val & ~indx_exc_sacc_pro & indx_instr_back_pro); 
selected.instr_back.all.anti.sup = indx_lateral(pval_sacc_anti_exc_sup_sign<p_val & ~indx_exc_sacc_anti & indx_instr_back_anti);

selected.instr.pro.exc = indx_lateral(indx_pro_only_sacc_exc_lat & indx_instr_pro); 
selected.instr.pro.sup = indx_lateral(indx_pro_only_sacc_sup_lat & indx_instr_pro); 
selected.instr.anti.exc = indx_lateral(indx_anti_only_sacc_exc_lat & indx_instr_anti);
selected.instr.anti.sup = indx_lateral(indx_anti_only_sacc_sup_lat & indx_instr_anti);
selected.instr.both.exc = indx_lateral(indx_pro_anti_sacc_exc_lat & indx_instr_pro & indx_instr_anti);
selected.instr.both.sup = indx_lateral(indx_pro_anti_sacc_sup_lat & indx_instr_pro & indx_instr_anti);
selected.instr.both.mixed = indx_lateral(indx_mixed);
selected.instr.all.pro.exc = indx_lateral(pval_sacc_pro_exc_sup_sign<p_val & ~indx_sup_sacc_pro & indx_instr_pro);
selected.instr.all.anti.exc = indx_lateral(pval_sacc_anti_exc_sup_sign<p_val & ~indx_sup_sacc_anti & indx_instr_anti);
selected.instr.all.pro.sup = indx_lateral(pval_sacc_pro_exc_sup_sign<p_val & ~indx_exc_sacc_pro & indx_instr_pro); 
selected.instr.all.anti.sup = indx_lateral(pval_sacc_anti_exc_sup_sign<p_val & ~indx_exc_sacc_anti & indx_instr_anti);


indx_lateral(~indx_pro_only_sacc_exc_lat & ~indx_pro_only_sacc_sup_lat & ~indx_anti_only_sacc_exc_lat & ~indx_anti_only_sacc_sup_lat & ~indx_pro_anti_sacc_exc_lat & ~indx_pro_anti_sacc_sup_lat & indx_instr_pro & indx_instr_anti); 
%% 
% pro_sacc_instrDir_omv = sum(indx_sign_saccVSinstrDir_nspk_pro & ~indx_sign_saccVSinstrDir_nspk_anti)
% anti_sacc_instrDir_omv = sum(indx_sign_saccVSinstrDir_nspk_anti & ~indx_sign_saccVSinstrDir_nspk_pro)
% pro_anti_sacc_instrDir_omv = sum(indx_sign_saccVSinstrDir_nspk_pro & indx_sign_saccVSinstrDir_nspk_anti)

lat_diff_pro_anti_sacc = sum(indx_sign_sacc); %%%%%%%%%%%%% to table
proportion_lat_sacc = lat_diff_pro_anti_sacc/length(indx_lateral); %%%%%%%%%%%%% to table
lat_diff_pro_anti_instr = sum(indx_sign_instr); %%%%%%%%%%%%% to table
proportion_lat_instr = lat_diff_pro_anti_instr/length(indx_lateral); %%%%%%%%%%%%% to table
lat_diff_sacc_base_pro = sum(indx_sign_sacc_base_pro); %%%%%%%%%%%%% to table
lat_diff_instr_base_pro = sum(indx_sign_instr_base_pro); %%%%%%%%%%%%% to table
lat_diff_sacc_base_anti = sum(indx_sign_sacc_base_anti); %%%%%%%%%%%%% to table
lat_diff_instr_base_anti = sum(indx_sign_instr_base_anti); %%%%%%%%%%%%% to table
lat_diff_instrDir_base_pro = sum(indx_sign_instrDir_instr_pro); %%%%%%%%%%%%% to table
lat_diff_instrDir_base_anti = sum(indx_sign_instrDir_instr_anti); %%%%%%%%%%%%% to table
lat_diff_goCue_instrDir_pro = sum(indx_sign_goCue_instrDir_pro); %%%%%%%%%%%%% to table
lat_diff_goCue_instrDir_anti = sum(indx_sign_goCue_instrDir_anti); %%%%%%%%%%%%% to table
lat_diff_saccVSinstrDir_nspk_pro = sum(indx_sign_saccVSinstrDir_nspk_pro); %%%%%%%%%%%%% to table
lat_diff_saccVSinstrDir_nspk_anti = sum(indx_sign_saccVSinstrDir_nspk_anti); %%%%%%%%%%%%% to table

%% instr
%indx_lat_sign = (indx_pro_only_sacc_exc_lat | indx_only_sacc_sup_lat | indx_anti_only_sacc_exc_lat | indx_anti_only_sacc_sup_lat | indx_pro_anti_sacc_exc_lat | indx_pro_anti_sacc_sup_lat); 
% instr_pro_pick = sum(indx_vermis_sign & indx_sign_instr_base_pro & ~indx_sign_instr_base_anti); 
% instr_pro_exc_pick = sum(indx_lat_sign & indx_sign_instr_base_pro & ~indx_sign_instr_base_anti & indx_exc_instr_pro);
% instr_pro_sup_pick = sum(indx_lat_sign & indx_sign_instr_base_pro & ~indx_sign_instr_base_anti & indx_sup_instr_pro);
% instr_anti_pick = sum(indx_lat_sign & ~indx_sign_instr_base_pro & indx_sign_instr_base_anti); 
% instr_anti_exc_pick = sum(indx_lat_sign & ~indx_sign_instr_base_pro & indx_sign_instr_base_anti & indx_exc_instr_anti); 
% instr_anti_supp_pick = sum(indx_lat_sign & ~indx_sign_instr_base_pro & indx_sign_instr_base_anti & indx_sup_instr_anti);
% instr_pro_anti_pick = sum(indx_lat_sign & indx_sign_instr_base_pro & indx_sign_instr_base_anti); 
% instr_pro_anti_exc_pick = sum(indx_lat_sign & indx_sign_instr_base_pro & indx_sign_instr_base_anti & indx_exc_instr_pro);
% instr_pro_anti_sup_pick = sum(indx_lat_sign & indx_sign_instr_base_pro & indx_sign_instr_base_anti & indx_sup_instr_pro);
% instr_sign_pick = sum(indx_lat_sign & indx_sign_instr );
%% 
  % pro only omv
  % anti only omv
  % pro only lat
  % anti only lat

%%
% print table
% disp('>>> VERMIS <<<')
% disp(['Total num of cells = ' num2str(total_cells)]); 
% disp(['Total num vermis = ' num2str(total_vermis)]); 
% disp(['Sign pro vs anti instr = ' num2str(vermis_diff_pro_anti_instr)]);
% disp(['Proportion sign pro vs anti instr = ' num2str(proportion_vermis_instr)]);
% disp(['Sign pro vs anti sacc = ' num2str(vermis_diff_pro_anti_sacc)]);
% disp(['Proportion sign pro vs anti sacc = ' num2str(proportion_vermis_sacc)]);
% disp(['Sign diff base vs instr pro = ' num2str(vermis_diff_instr_base_pro)]);
% disp(['Sign diff base vs sacc pro = ' num2str(vermis_diff_sacc_base_pro)]);
% disp(['Sign diff base vs instr anti = ' num2str(vermis_diff_instr_base_anti)]);
% disp(['Sign diff base vs sacc anti = ' num2str(vermis_diff_sacc_base_anti)]);
% disp(['Sign diff base vs instrDir pro = ' num2str(vermis_diff_instrDir_base_pro)]);
% disp(['Sign diff base vs instrDir anti = ' num2str(vermis_diff_instrDir_base_anti)]);
% disp(['Sign diff instrDir vs goCue pro = ' num2str(vermis_diff_goCue_instrDir_pro)]);
% disp(['Sign diff instrDir vs goCue anti = ' num2str(vermis_diff_goCue_instrDir_anti)]);
% disp('                             ')
% disp('>>> LATERAL <<<')
% disp(['Total num of cells = ' num2str(total_cells)]); 
% disp(['Total num lateral = ' num2str(total_lateral)]); 
% disp(['Sign pro vs anti instr = ' num2str(lateral_diff_pro_anti_instr)]);
% disp(['Proportion sign pro vs anti instr = ' num2str(proportion_lateral_instr)]);
% disp(['Sign pro vs anti sacc = ' num2str(lateral_diff_pro_anti_sacc)]);
% disp(['Proportion sign pro vs anti sacc = ' num2str(proportion_lateral_sacc)]);
% disp(['Sign diff base vs instr pro = ' num2str(lateral_diff_instr_base_pro)]);
% disp(['Sign diff base vs sacc pro = ' num2str(lateral_diff_sacc_base_pro)]);
% disp(['Sign diff base vs instr anti = ' num2str(lateral_diff_instr_base_anti)]);
% disp(['Sign diff base vs sacc anti = ' num2str(lateral_diff_sacc_base_anti)]);
% disp(['Sign diff base vs instrDir pro = ' num2str(lateral_diff_instrDir_base_pro)]);
% disp(['Sign diff base vs instrDir anti = ' num2str(lateral_diff_instrDir_base_anti)]);
% disp(['Sign diff instrDir vs goCue pro = ' num2str(lateral_diff_goCue_instrDir_pro)]);
% disp(['Sign diff instrDir vs goCue anti = ' num2str(lateral_diff_goCue_instrDir_anti)]);