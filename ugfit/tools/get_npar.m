function [ npar ] = get_npar( modelID)
% Lookup table to get number of free parameters per model
% Created by MKW 2018
% Adapted by BRK Shevlin 2023


%%%%% main models
if       strcmp(modelID,'ms_UG0_adaptiveNorm'),            npar = 4;
elseif       strcmp(modelID,'ms_UG1_etaf_adaptiveNorm'),   npar = 5;
elseif       strcmp(modelID,'ms_UG2_etaf_adaptiveNorm'),   npar = 5;
elseif       strcmp(modelID,'ms_UG3_etaf_adaptiveNorm'),   npar = 5;
elseif       strcmp(modelID,'ms_UG0_f0f_adaptiveNorm'),   npar = 3;
elseif       strcmp(modelID,'ms_UG1_etaf_f0f_adaptiveNorm'),   npar = 4;
elseif       strcmp(modelID,'ms_UG2_etaf_f0f_adaptiveNorm'),   npar = 4;
elseif       strcmp(modelID,'ms_UG0_f0f_adaptiveNorm_blocked'),   npar = 4;
elseif       strcmp(modelID,'ms_UG1_etaf_f0f_adaptiveNorm_blocked'),   npar = 4;
elseif       strcmp(modelID,'ms_UG2_etaf_f0f_adaptiveNorm_blocked'),   npar = 4;
elseif       strcmp(modelID,'ms_UG3_etaf_f0f_adaptiveNorm_blocked'),   npar = 4;
elseif       strcmp(modelID,'ms_UG3_etaf_f0f_adaptiveNorm'),   npar = 4;
elseif       strcmp(modelID,'ms_UG0_fixedNorm'),           npar = 3;
elseif       strcmp(modelID,'ms_UG1_etaf_fixedNorm'),      npar = 4;
elseif       strcmp(modelID,'ms_UG2_etaf_fixedNorm'),      npar = 4;
elseif       strcmp(modelID,'ms_UG3_etaf_fixedNorm'),      npar = 4;

end




end

