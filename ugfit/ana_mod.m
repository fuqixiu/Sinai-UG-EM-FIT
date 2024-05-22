% Fits UG2 models and does model comparison
% Created by MK Wittmann, October 2018
% Adapted by BRK Shevlin, April 2023
%
%%
%== -I) Prepare workspace: ============================================================================================
clearvars
addpath('models'); 
addpath('tools');
setFigDefaults; 
%%

%%== 0) Load and organise data: ==========================================================================================
% load data:
%indir = '~example_data';
indir = 'C:\Users\fuq01\Documents\GitHub\Sinai-UG-EM-FIT\example_data';
load(fullfile(indir,'DATA_LEAP_online_baseline_2024.mat'));
% define data set(s) of interest:
%expids = {'hc','dep', 'anh', 'both'}; % fit by groups version
expids = {'all'}; % fit by all version
% how to fit RL:
M.dofit     = 1;                                                                                                     % whether to fit or not                                                           
M.doMC      = 1;                                                                                                     % whether to do model comparison or not  
M.quickfit  = 0;   % Shawn tells me that quickfit is too liberal to be trusted, but is fine for testing                                                                                                  % whether to lower the convergence criterion for model fitting (makes fitting quicker) (1=quick fit)
M.omitBMS   = 0;   % omit bayesian model comparison if you don't have SPM instaslled
M.modid     = { 'ms_UG0_f0f_adaptiveNorm', ...
     'ms_UG0_adaptiveNorm','ms_UG0_fixedNorm'}; 
              
                % list of main models to fit

%%
%== I) RUN MODELS: ======================================================================================================
trials = (1:60); %online 2 blocks, 60 trials
for iexp = 1:numel(expids)
   if M.dofit == 0,  break; end
   cur_exp = expids{iexp};                                                   
   
   %%% EM fit %%%
   for im = 1:numel(M.modid)
      dotry=1;
      while 1==dotry
         try
            close all;
            % If things are not working, run this line manually to set an
            % idea where the error in the code is
            s.(cur_exp).em = EMfit_changept_ms(s.(cur_exp),M.modid{im},M.quickfit,trials);
            %s.(cur_exp).em = EMfit_fmri_ms(s.(cur_exp),M.modid{im},M.quickfit,trials);

            dotry=0;
         catch
            dotry=1; disp('caught');
         end
      end
   end

   %%% calc BICint for EM fit
   for im = 1:numel(M.modid)
      s.(cur_exp).em.(M.modid{im}).fit.bicint =  cal_BICint_ms(s.(cur_exp).em, M.modid{im},trials);     
   end   

end


save('WEIGHTED_FIT_LEAP_online_May_2024.mat','s')

cur_exp = 'all';
%%
%== II) COMPARE MODELS: ================================================================================================

M.modid = {'ms_UG0_f0f_adaptiveNorm', ...
     'ms_UG0_adaptiveNorm','ms_UG0_fixedNorm'};

for iexp = 1:numel(expids)
   if M.doMC~=1, break; end
   cur_exp = expids{iexp};
   EMmc_ms_Hess_fix(s.(cur_exp),M.modid);  
end


%% Look at params of winning model
modelID =  'ms_UG0_fixedNorm';
winning_model = s.(cur_exp).em.(modelID);

params = winning_model.q;

est_params = zeros(size(params,1),size(params,2));

for i = 1:size(params,1)
    est_params(i,1) = norm2alpha(params(i,1));
    est_params(i,2) = norm2beta(params(i,2));
    est_params(i,3) = norm2alpha(params(i,3));
    %est_params(i,4) = norm2delta(params(i,4));
end

%% ONLY RELEVANT FOR PARAMETER RECOVERY EXERCISES
% This loads in param_tru
load('recover_nRv_f3_cap2_t20_etaf_no0_IC.mat')


% Temperature
corrplot([est_params(:,2) param_tru(:,1)    ])
% Envy
corrplot([est_params(:,1) param_tru(:,2)    ])
% Adaptation
corrplot([est_params(:,3) param_tru(:,4)    ])
% Delta
corrplot([est_params(:,4) param_tru(:,5)    ])

