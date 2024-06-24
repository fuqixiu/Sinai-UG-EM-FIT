% Convert to format that is more R-friendly

load("example_fits\WEIGHTED_FIT_LEAP_online_June_2024_4G_9M.mat")

M.modid     = { 'ms_UG0_f0f_adaptiveNorm', 'ms_UG0_f0f_adaptiveNorm_v2', 'ms_UG0_adaptiveNorm', ...
                'ms_UG0_FS_10','ms_UG0_fixedNorm', 'ms_UG0_adaptiveNorm_fixedalpha', ...
                'ms_UG0_Bayes_5', 'ms_UG0_Bayes_10', 'ms_UG0_Bayes_variable'}; 

expList = fieldnames(s);
for xp = 1:size(fieldnames(s))                
    cur_exp = expList{xp};

    for im = 1:numel(M.modid)
        
       IDs = s.(cur_exp).ID;  
       BIC = s.(cur_exp).em.(M.modid{im}).fit.bic;
       q = s.(cur_exp).em.(M.modid{im}).q;
        
       params = zeros(size(q));
       for i = 1:size(IDs)
            qi = q(i,:);
            params(i,:) = norm2par(M.modid{im},qi);
       end
    
       model_name_save = ['WEIGHTED_FIT_LEAP_online_' cur_exp '_t60_noFlat_' M.modid{im} '_June24.csv'];
   
       T= table(IDs,BIC,params);
    
       writetable(T,model_name_save )

    end
end