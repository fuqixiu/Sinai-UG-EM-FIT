function [fval,fit] = mod_ms_UG0_Bayes_variable(behavData,q, doprior,dofit,varargin)
% runs standard 0 Step UG Model
% Created MK Wittmann, Oct 2018
% Adapted BRK Shevlin, April 2023
% Written AD July 2023
% Bayesian Observer variable initial Norm mu
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

[qt, bounds] = norm2par('ms_UG0_Bayes_variable',q); % transform parameters from gaussian space to model space

% Define free parameters and set unused ones to zero
alpha = qt(1); % envy
if (alpha<min(bounds(:,1)) || alpha>max(bounds(:,1))), fval=10000000; return; end
beta = qt(2); % inverse temperature
if (beta<min(bounds(:,2)) || beta>max(bounds(:,2))), fval=10000000; return; end
f0 = qt(3); % initial norm
if (f0<min(bounds(:,3)) || f0>max(bounds(:,3))), fval=10000000; return; end

delta = 0; % expected influence

fixed = 4; % Fixed k paramter
free = {alpha beta f0};
% -------------------------------------------------------------------------------------
% 2-4) Middle code is specific the the model
% -------------------------------------------------------------------------------------
if doprior == 1
    [fval,norm,V,ChoiceProb] = lik_UG0_Bayes_variable(behavData.offer, behavData.choice,fixed,free,doprior,prior,q);
else
    [fval,norm,V,ChoiceProb] = lik_UG0_Bayes_variable(behavData.offer, behavData.choice,fixed,free,doprior);

end
% -------------------------------------------------------------------------------------
% 5) Calculate additional Parameters and save: 
% -------------------------------------------------------------------------------------

if dofit ==1

   fit         = struct;
   fit.xnames  = {'alpha'; 'beta'; 'f0'};
   
   fit.mat    = [ChoiceProb norm V];
   fit.names  = {'ChoiceProb' ; 'norm'; 'V'};
end




end

