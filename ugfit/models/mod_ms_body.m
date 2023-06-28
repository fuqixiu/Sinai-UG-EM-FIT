% Body of RL model that can run CL-trace, CS-trace and RT-trace. Compatible with all RL models.
% Ran directly from modms_RL_cl_cs_rt etc scripts
% MKW, Oct 2018
% Adapted BRK Shevlin April 2023
%%

% -------------------------------------------------------------------------------------
% 1 ) Define stuff, get variables
% -------------------------------------------------------------------------------------

% define behaviour on which model is fitted:
outcome        = behavData.choice;                                       
offers         = behavData.offer;                       
reward         = behavData.reward;



% prepare variables to collect information of interest
num_opt     = 2;
allvals     = nan(numel(outcome),num_opt);
ChoiceProb  = nan(numel(outcome),1);
o1_val      = nan(numel(outcome),1);  
o2_val      = nan(numel(outcome),1);
PE          = nan(numel(outcome),1);
perew       = nan(numel(outcome),1);
cloc        = nan(numel(outcome),1);
ch_trace    = nan(numel(outcome),3);
csbonus     = zeros(numel(outcome),2);
rew_trace   = nan(numel(outcome),1);

% define starting values:
EstState       = repmat(.5,1,num_opt);                                                    % order: [option1 option2 option3] % was 0 before sep18
EstChTrace     = [0 0 0];
allvals(1,:)   = EstState;
cloc(1)        = 0;
rew_trace(1)   = 0 ;

% -------------------------------------------------------------------------------------x
% 2 ) Learning model: 
% -------------------------------------------------------------------------------------


for it = 1:numel(outcome)
   
    % 1) get offer IDs and values:
    EstState      = allvals(it,:);
    o1_id         = offers(it,1);         o2_id    = offers(it,2);         
    o1_val(it)    = EstState(o1_id);  o2_val(it)   = EstState(o2_id); 
       
    % 2) Update basic RL part:
    norm = RW(f0,epsilon,offers);
    
    % 3) Update CL-trace
    cPE         = choice_loc_scale(it) - cloc(it);
    cloc(it+1)  = cloc(it)+lr_cl*cPE;   
    
    % 4) Update CS-trace 
    for itrace = 1:3
       ch_trace(it,itrace) = EstChTrace(itrace) * tau_cs;
    end
    csbonus(it,:)                = [ ch_trace(it,o1_id) ch_trace(it,o2_id)];
    EstChTrace                   = ch_trace(it,:);
    EstChTrace(choice_id(it))    = 1;                    

    % 5) Update reward trace
    perew(it)           = outcome(it) - rew_trace(it);
    rew_trace(it+1)     = rew_trace(it) + lr_rt*perew(it);
  
end
allvals     = allvals(1:size(outcome,1),:); % cut off last outcome
rew_trace   = rew_trace(1:size(outcome,1));
cloc        = cloc(1:numel(outcome));



% -------------------------------------------------------------------------------------
% 3 ) Observation model:  
% -------------------------------------------------------------------------------------

DVbasic = o1_val - o2_val - sblr_bias;                                        % coded in terms of evidence for option 1
DV      = DVbasic - wcl*cloc + wcs*(csbonus(:,1)-csbonus(:,2));

% calculate probability of chosing 1:                                                                                             
o1_prob     = 1 ./ ( 1 + exp(-beta.*DV));                          

% find when 1 was actually choosen:
pick1 = find(choice_id == offers(:,1));
pick2 = find(choice_id == offers(:,2));

ChoiceProb(pick1) = o1_prob(pick1);
ChoiceProb(pick2) = 1 - o1_prob(pick2);


% -------------------------------------------------------------------------------------
% 4 ) Calculate model fit:  
% -------------------------------------------------------------------------------------

nll =-nansum(log(ChoiceProb));                                                % the thing to minimize                      

if doprior == 0                                                               % NLL fit
   fval = nll;
elseif doprior == 1                                                           % EM-fit:   P(Choices | h) * P(h | O) should be maximised, therefore same as minimizing it with negative sign   
   fval = -(-nll + prior.logpdf(q));
end

if sum(isnan(ChoiceProb))>0, disp('ERROR'); keyboard; return; end             % error check                  











