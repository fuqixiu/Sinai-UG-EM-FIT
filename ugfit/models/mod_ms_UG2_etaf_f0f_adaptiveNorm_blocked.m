function [fval,fit] = mod_ms_UG2_etaf_f0f_adaptiveNorm_blocked(behavData,q, doprior,dofit,varargin)
% runs standard 0 Step UG Model
% Created MK Wittmann, Oct 2018
% Adapted BRK Shevlin, April 2023
%
% INPUT:    - behavData: behavioural input file
% OUTPUT:   - fval and fitted variables
% 
% 
%
%%
% -------------------------------------------------------------------------------------
% 1 ) Define free parameters
% -------------------------------------------------------------------------------------

if nargin > 4
    prior      = varargin{1};
end

[qt, bounds] = norm2par('ms_UG2_etaf_f0f_adaptiveNorm',q); % transform parameters from gaussian space to model space

% Define free parameters and set unused ones to zero
alpha = qt(1); % envy
if (alpha<min(bounds(:,1)) || alpha>max(bounds(:,1))), fval=10000000; return; end
beta = qt(2); % inverse temperature
if (beta<min(bounds(:,2)) || beta>max(bounds(:,2))), fval=10000000; return; end

epsilon = qt(3); % norm adaptation rate
if (epsilon<min(bounds(:,3)) || epsilon>max(bounds(:,3))), fval=10000000; return; end

delta = qt(4); % expected influence
if (delta<min(bounds(:,4)) || delta>max(bounds(:,4))), fval=10000000; return; end


fixed = [0.8, 10];
free = {alpha beta epsilon delta};
% -------------------------------------------------------------------------------------
% 2-4) Middle code is specific the the model
% -------------------------------------------------------------------------------------
if doprior == 1
    [fval,norm,V,ChoiceProb] = lik_UG2_etaf_f0f_adaptiveNorm_blocked(behavData.offer, behavData.choice,fixed,free,doprior,prior,q);
else
    [fval,norm,V,ChoiceProb] = lik_UG2_etaf_f0f_adaptiveNorm_blocked(behavData.offer, behavData.choice,fixed,free,doprior);

end
% -------------------------------------------------------------------------------------
% 5) Calculate additional Parameters and save: 
% -------------------------------------------------------------------------------------

if dofit ==1

   fit         = struct;
   fit.q = qt;
   fit.xnames  = {'alpha'; 'beta';'epsilon';'delta'};
   
   fit.ChoiceProb    = ChoiceProb;
   fit.norm  = norm;
   fit.V = V;
end




end

