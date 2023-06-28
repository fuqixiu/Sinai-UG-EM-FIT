% Convert to format that is more R-friendly

load('WEIGHTED_FIT_ED_t30_May24_2023.mat')

M.modid     = {'ms_UG0_f0f_adaptiveNorm',...
                 'ms_UG1_etaf_f0f_adaptiveNorm',... 
                 'ms_UG2_etaf_f0f_adaptiveNorm',...
                 'ms_UG3_etaf_f0f_adaptiveNorm' };
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
    
       model_name_save = ['WEIGHTED_FIT_ED_Round1_' cur_exp '_t30_noFlat_' M.modid{im} '_May24.csv'];
   
       T= table(IDs,BIC,params);
    
       writetable(T,model_name_save )

    end
end